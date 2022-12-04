/*
compDist_merge.cu

comparing distributions organized into the rows of two matrices A and B.
This is part of the computation of the Third Lower Bound (TLB) of network Gromov-Wasserstain distance 
per the work of Chowdhury & Memoli (2019).

This is a parallelization of part of the function "tools/compareRealDistibutions.m" in the repo (https://github.com/samirchowdhury/GWnets).
Here we assume, however, the inputs are already sorted, i.e. A, B, and the corresponding cumulative probability measures cmA, cmB. 

Compare to an earlier version of this program "compDist_mf.cu" (9/6/2019), the present version use a more efficient merge before sorting. 

This is not designed to compare matrices greater than 2000x2000. 
More precisely, 2 + 4 * (NA + NB) * sizeof(float) should not exceed 48k.
I.e. NA + NB <= 3071.

== to compile in MATLAB
mexcuda compDist_merge.cu
----------------
created by Mengsen Zhang, mengsenzhang@gmail.com (9/18/2019)

*/
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "mex.h"

using namespace std;

// #define UNIQUETOL 1e-10 // uncomment if you want to add a "unique" procedure while computing the joint distribution; also uncomment corresponding lines in "mergesort".
#define N_REGS_PER_THREAD 52 // remember to update this if the program is change!

cudaError_t compdist(const double* A, const double* B, const double* cmA, const double* cmB, const unsigned int NA, const unsigned int NB, double* dist);
unsigned int nextPow2(unsigned int n);
__global__ void compdistKernel(const double* A, const double* B, const double* CMA, const double* CMB, const unsigned int NA, const unsigned int NB, double* dist);
__device__ void mergesort(const float* x, const float* y, float* xy, int* zeroflag, const int nx, const int ny);
__device__ int lastPow2(int n);
__device__ void findIdx(const float* cm, const float* CM, int* idx, const int N_cm, const int N_CM);
template<typename num> __device__ num sum(num* x, const int len);
template<typename num> __device__ void reset(num* x, const int len);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/*
	Interface with matlab, need 6 inputs, and 1 output.
	*/
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("mexFun:nlhs", "need 1 output: dist");
	}
	if (nrhs != 6) {
		mexErrMsgIdAndTxt("mexFun:nrhs", "need 6 inputs: sorted_A, sorted_B, sorted_cmA, sorted_cmB, NA, NB");
	}

	// read input
	double *A = mxGetPr(prhs[0]);
	double *B = mxGetPr(prhs[1]);
	double *cmA = mxGetPr(prhs[2]);
	double *cmB = mxGetPr(prhs[3]);
	int NA = mxGetScalar(prhs[4]);
	int NB = mxGetScalar(prhs[5]);

	// prep output
	plhs[0] = mxCreateDoubleMatrix(NB, NA, mxREAL);
	double *dist = mxGetPr(plhs[0]);

	// compute!
	cudaError_t cudaStatus = compdist(A, B, cmA, cmB, (unsigned int) NA, (unsigned int) NB, dist);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:compdist", "main function didn't run correctly!");
	}
}

// host function for launching the kernel
cudaError_t compdist(const double* A, const double* B, const double* cmA, const double* cmB, 
					 const unsigned int NA, const unsigned int NB, 
					 double* dist) {
	const unsigned int N_max = NA + NB + 1;
	// device data points
	double *dev_cmA, *dev_cmB, *dev_A, *dev_B, *dev_dist;

	// check device
	cudaError_t cudaStatus = cudaSuccess;
	cudaDeviceReset();
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	
	// declare variables
	int count;
	size_t requiredsharemem;
	size_t dyn_sharedmemsize;
	const unsigned int max_regs_threads = (prop.regsPerBlock-1)/N_REGS_PER_THREAD + 1;//max # threads per block limited by # registers
	const unsigned int N_optim_threads = min(N_max, prop.warpSize*prop.multiProcessorCount);//optimal threads per block for all MP to be working
	const unsigned int blocksize = min(N_optim_threads, min(max_regs_threads,prop.maxThreadsPerBlock));
	// mexPrintf("blocksize=%d\n", blocksize);
	dim3 griddims(NA, NB, 1);
	dim3 blockdims(blocksize, 1, 1);

	cudaStatus = cudaGetDeviceCount(&count);
	if (cudaStatus != cudaSuccess) {
		mexPrintf("Number of device: %d\n", count);
		mexErrMsgIdAndTxt("compDist:devicecount ", "Could not find device.!");
		goto Error;
	}

	// check shared memory size
	requiredsharemem = 2 + 4 * (NA + NB) * sizeof(float);
	
	cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
	dyn_sharedmemsize = prop.sharedMemPerBlock;
	if (requiredsharemem > dyn_sharedmemsize) {
		mexPrintf("Requested shared memory %d bytes > available %d bytes.\n", requiredsharemem, dyn_sharedmemsize);
		mexWarnMsgTxt("The shared memory required for your job may exceed the capacity of your GPU.");
	}

	// choose device
	if (count>0) {
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			mexErrMsgIdAndTxt("compdist:cudaSetDevice", "cannot set device 0!\n");
			goto Error;
		}
	}

	// allocate memory at device
	cudaStatus = cudaMalloc((void**)&dev_cmA, NA * NA * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Malloc dev_cmA!\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_cmB, NB * NB * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Malloc dev_cmB!\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_A, NA * NA * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Malloc dev_A!\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_B, NB * NB * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Malloc dev_B!\n");
		goto Error;
	} 
	cudaStatus = cudaMalloc((void**)&dev_dist, NA * NB * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Malloc dev_dist!\n");
		goto Error;
	} 
    cudaStatus = cudaMemset(dev_dist, 99.0, NA * NB * sizeof(double));
    if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Memset dev_dist!\n");
		goto Error;
	}
	// transfer data
	cudaStatus = cudaMemcpy(dev_cmA, cmA, NA * NA * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMemcpyHostToDevice", "Can't transfer cmA to device!\n");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_cmB, cmB, NB * NB * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMemcpyHostToDevice", "Can't transfer cmB to device!\n");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_A, A, NA * NA * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMemcpyHostToDevice", "Can't transfer A to device!\n");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_B, B, NB * NB * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMemcpyHostToDevice", "Can't transfer A to device!\n");
		goto Error;
	}

	// launch kernel	

	compdistKernel<<<griddims, blockdims, min(requiredsharemem,dyn_sharedmemsize)>>>(dev_A, dev_B, dev_cmA, dev_cmB, NA, NB, dev_dist);

	// check for errors
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		mexPrintf("cudaGetLastError returned error code %d after launching compdistKernel!\n", cudaStatus);
		mexErrMsgIdAndTxt("compdist:cudaGetLastError", "kernel didn't run correctly!\n");
		goto Error;
	}
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		mexPrintf("cudaDeviceSynchronize returned error code %d after launching compdistKernel!\n", cudaStatus);
		mexErrMsgIdAndTxt("compdist:cudaDeviceSynchronize", "device didn't sync!\n");
		goto Error;
	}

	// copy results to host
	cudaStatus = cudaMemcpy(dist, dev_dist, NA * NB * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMemcpyDeviceToHost", "couldn't get data back!\n");
		goto Error;
	}

	cudaDeviceReset();
Error:
	cudaDeviceReset();
	return cudaStatus;
}

unsigned int nextPow2(unsigned int n) {
	n--;
	n |= n >> 1;
	n |= n >> 2;
	n |= n >> 4;
	n |= n >> 8;
	n |= n >> 16;
	n++;
	return n;
}

__global__ void compdistKernel(const double* A, const double* B, 
								const double* CMA, const double* CMB,
								const unsigned int NA, const unsigned int NB,
								double* dist) {
	// basic info
	int idx_vA0 = blockIdx.x*NA; // start of the ii-th row of A
	int idx_vB0 = blockIdx.y*NB; // start of the jj-th row of B
	int N_merged = NA + NB + 1;

	// prepare shared memory
	extern __shared__ float ins[];
	float *cmA = ins;
	float *cmB = (float*)&cmA[NA]; // # in brack is the size of the previous variable!
	float *cmMerged = (float*)&cmB[NB];
	int *zeroflag = (int*)&cmMerged[N_merged];
	float *diffAB = (float*)&zeroflag[N_merged];

	// load cdf's from global memory
	int N_steps = (NA - 1) / blockDim.x + 1;// for large vectors, we load in steps
	for (int n = 0; n < N_steps; n++) {
		int idx_v = n*blockDim.x + threadIdx.x;
		int idx_A = idx_vA0 + idx_v;
		if (idx_v < NA) {
			cmA[idx_v] = __double2float_rn(CMA[idx_A]);
		}
	}
	__syncthreads();
	N_steps = (NB - 1) / blockDim.x + 1;// for large vectors, we load in steps
	for (int n = 0; n < N_steps; n++) {
		int idx_v = n*blockDim.x + threadIdx.x;
		int idx_B = idx_vB0 + idx_v;
		if (idx_v < NB) {
			cmB[idx_v] = __double2float_rn(CMB[idx_B]);
		}
	}
	__syncthreads();

	// merge and sort cdf's
	reset<int>(zeroflag,N_merged);
	reset<float>(cmMerged,N_merged);
	__syncthreads();
	if (NA>NB) {
		mergesort(cmA, cmB, cmMerged, zeroflag, NA, NB);// the algorithm assume the first argument is larger
	} else{
		mergesort(cmB, cmA, cmMerged, zeroflag, NB, NA);
	}
	
	__syncthreads();

	int N_dCM = N_merged - 1;
	reset<int>(zeroflag,N_merged);
	// finding indices for vA (reuse "zeroflag" to store indices)
	findIdx(cmA, cmMerged, zeroflag, NA, N_dCM);
	// load A by idx_vA
	N_steps = (N_dCM - 1) / blockDim.x + 1;
	for (int n = 0; n < N_steps; n++) {
		int idx_idx_vA = n*blockDim.x + threadIdx.x;
		if (idx_idx_vA < N_dCM) {
			diffAB[idx_idx_vA] = __double2float_rn(A[idx_vA0 + zeroflag[idx_idx_vA]]);// reuse cmA as vA
		}
	}
	__syncthreads();

	// finding indices for vB
	findIdx(cmB, cmMerged, zeroflag, NB, N_dCM);

	// load B by idx_vB and subtract A
	for (int n = 0; n < N_steps; n++) {
		int idx_idx_vB = n*blockDim.x + threadIdx.x;
		if (idx_idx_vB < N_dCM) {
			diffAB[idx_idx_vB] = powf(diffAB[idx_idx_vB] - __double2float_rn(B[idx_vB0 + zeroflag[idx_idx_vB]]), 2)*(cmMerged[idx_idx_vB + 1] - cmMerged[idx_idx_vB]);
		}
	}
	__syncthreads();

	// output final distance
	double dist_ij = (double) sum<float>(diffAB, N_dCM); 
	__syncthreads();
	if (threadIdx.x == 0) {
		dist[blockIdx.x*NB + blockIdx.y] = sqrtf(dist_ij);//cmMerged[(blockIdx.x*NB + blockIdx.y)%N_merged];//
	}
	__syncthreads();
}

__device__ void findIdx(const float* cm, const float* CM, int* idx, const int N_cm, const int N_CM) {
	reset<int>(idx, N_CM);
	int N_steps = (N_CM - 1) / blockDim.x + 1;
	// for each cm : compare to CM values in parallel
	for (int n = 0; n < N_steps; n++) {
		int idx_CM = n*blockDim.x + threadIdx.x;
		if ( idx_CM < N_CM ) {
			for (int m = 0; m < N_cm; m++) {// loop through values of cm
				idx[idx_CM] = min(idx[idx_CM] + (cm[m] <= CM[idx_CM]), N_cm - 1);//count the number of cm <=CM_i
				__syncthreads();
			}
		}
		__syncthreads();
	}	
	__syncthreads();
}

__device__ void mergesort(const float* x, const float* y,
	float* xy, int* zeroflag,
	const int nx, const int ny) {
	int N = nx + ny; //max N
	int Nsteps = (N - 1) / blockDim.x + 1;
	// -- merge: interleaving two vectors assuming nx>ny
	// distribute x over a NA+NB long vector
	for (int n = 0; n<Nsteps; n++) {
		int idx_x = n*blockDim.x + threadIdx.x;
		if (idx_x < nx) {
			xy[(int) (idx_x*N)/nx + 1] = x[idx_x];
		}
	}
	__syncthreads();
	// distribute y between x
	for (int n = 0; n<Nsteps; n++) {
		int idx_xy = n*blockDim.x + threadIdx.x;
		if ((idx_xy + 1 < N) && (xy[idx_xy + 2] == 0) ) {// if there is a gap to the right of idx_xy
			xy[idx_xy + 2] = y[idx_xy*ny/N];//insert y
		}
	}
	__syncthreads();

	// sort
	int N_comp = (N + 1) / 2;// max number of comparisons per loop
	Nsteps = (N_comp - 1) / blockDim.x + 1;
	int change;
	do {
		// use zeroflag to keep track of change
		reset<int>(zeroflag, N_comp);
		change = 0;
		__syncthreads();

		// comparing neighboring points
		for (int n = 0; n<Nsteps; n++) {
			int idxcomp = n*blockDim.x + threadIdx.x;			
			if (idxcomp * 2 + 1 < N + 1) {
				// compare : round 1
				float left = xy[idxcomp * 2];
				float right = xy[idxcomp * 2 + 1];
				// bool sim = (abs(left - right) < UNIQUETOL);
				xy[idxcomp * 2] = min(left, right);// * (float) (!sim);// zero if two numbers are similar
				xy[idxcomp * 2 + 1] = max(left, right);// * (float) (!(sim && (min(left, right) == 0)));
				zeroflag[idxcomp] += (left != xy[idxcomp * 2]);// || (right != xy[idxcomp * 2 + 1]);
			}
		}
		__syncthreads();
		for (int n = 0; n<Nsteps; n++) {
			int idxcomp = n*blockDim.x + threadIdx.x;
			if (idxcomp * 2 + 2 < N + 1) {
				// round 2
				float left = xy[idxcomp * 2 + 1];
				float right = xy[idxcomp * 2 + 2];
				// bool sim = (abs(left - right) < UNIQUETOL);
				xy[idxcomp * 2 + 1] = min(left, right);// * (float) (!sim);// zero if two numbers are similar
				xy[idxcomp * 2 + 2] = max(left, right);// * (float) (!(sim && (min(left, right) == 0)));
				zeroflag[idxcomp] += (left != xy[idxcomp * 2 + 1]);// || (right != xy[idxcomp * 2 + 2]);
			}
		}
		__syncthreads();
		change = sum<int>(zeroflag, N_comp);
		__syncthreads();
	} while (change);
}

// the largest power of 2 smaller than input
__device__ int lastPow2(int n) {
	// next power of 2
	n--;
	n |= n >> 1;
	n |= n >> 2;
	n |= n >> 4;
	n |= n >> 8;
	n |= n >> 16;
	n++;
	// last power of 2
	return n >> 1;
}

// sum elements along x within a block
template<typename num>
__device__ num sum(num* x, const int len) {
	int n_sums0 = lastPow2(len); // max number of cols in the initial sum step
	for (int n_sums = n_sums0; n_sums > 0; n_sums >>= 1) {
		int N_steps = (n_sums - 1) / blockDim.x + 1;
		for (int n = 0; n<N_steps; n++) {
			int idx = n*blockDim.x + threadIdx.x;
			if (idx < n_sums && (idx + n_sums) < len) {
				x[idx] += x[idx + n_sums];
			}
			__syncthreads();
		}
		__syncthreads();
	}
	__syncthreads();
	return x[0];
}

// reset all elements to zero
template<typename num>
__device__ void reset(num* x, const int len) {
	int N_steps = (len - 1) / blockDim.x + 1;
	for (int n = 0; n < N_steps; n++) {
		int idx = n*blockDim.x + threadIdx.x;
		if (idx < len) {
			x[idx] = 0;
		}
		__syncthreads();
	}
}
