/*
compDist_uni.cu

comparing distributions organized into the rows of two matrices A and B.
This is part of the computation of the Third Lower Bound (TLB) of network Gromov-Wasserstain distance 
per the work of Chowdhury & Memoli (2019).

Here we assume probability distributions mA and mB are uniform. This greatly simplify the computation. 


== to compile in MATLAB
mexcuda compDist_uni.cu

== warning
this is not designed to handle comparisons greater than 4000*4000;
----------------
created by Mengsen Zhang, mengsenzhang@gmail.com (9/8/2019).

*/
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "mex.h"

using namespace std;

#define N_THREADS_PER_BLOCK 1024

cudaError_t compdist(const double* A, const double* B, const double* idxA, const double* idxB, const double* dCM,
					 const unsigned int NA, const unsigned int NB, const unsigned int NdCM,
					 double* dist);
__global__ void compdistKernel(const double* A, const double* B, 
								const double* idxA, const double* idxB, const double* dCM,
								const unsigned int NA, const unsigned int NB, const unsigned int NdCM,
								double* dist);
__device__ float sum(float* x, const int len);
__device__ int lastPow2(int n);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/*
	Interface with matlab, need 6 inputs, and 1 output.
	*/
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("mexFun:nlhs", "need 1 output: dist");
	}
	if (nrhs != 8) {
		mexErrMsgIdAndTxt("mexFun:nrhs", "need 8 inputs: sorted_A, sorted_B, idxA, idxB, dCM, NA, NB, NdCM");
	}

	// read input
	double *A = mxGetPr(prhs[0]);
	double *B = mxGetPr(prhs[1]);
	double *idxA = mxGetPr(prhs[2]);
	double *idxB = mxGetPr(prhs[3]);
	double *dCM = mxGetPr(prhs[4]);
	int NA = mxGetScalar(prhs[5]);
	int NB = mxGetScalar(prhs[6]);
	int NdCM = mxGetScalar(prhs[7]);

	// prep output
	plhs[0] = mxCreateDoubleMatrix(NB, NA, mxREAL);
	double *dist = mxGetPr(plhs[0]);

	// compute!
	cudaError_t cudaStatus = compdist(A, B, idxA, idxB, dCM, (unsigned int) NA, (unsigned int) NB, (unsigned int) NdCM, dist);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:compdist", "main function didn't run correctly!");
	}
}

// host function for launching the kernel
cudaError_t compdist(const double* A, const double* B, const double* idxA, const double* idxB, const double* dCM,
					 const unsigned int NA, const unsigned int NB, const unsigned int NdCM,
					 double* dist) {
	
	// device data points
	double *dev_A, *dev_B, *dev_dist, *dev_idxA, *dev_idxB, *dev_dCM;

	// check device
	cudaError_t cudaStatus = cudaSuccess;
	cudaDeviceReset();
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	
	// declare variables
	int count;
	size_t requiredsharemem;
	size_t dyn_sharedmemsize;
	const unsigned int blocksize = min(NdCM, N_THREADS_PER_BLOCK);
	// size_t smsize = (size_t) min(requiredsharemem,dyn_sharedmemsize);
	dim3 griddims(NA, NB, 1);
	dim3 blockdims(blocksize, 1, 1);

	cudaStatus = cudaGetDeviceCount(&count);
	if (cudaStatus != cudaSuccess) {
		mexPrintf("Number of device: %d\n", count);
		mexErrMsgIdAndTxt("compDist:devicecount ", "Could not find device.!");
		goto Error;
	}

	// check shared memory size
	requiredsharemem = (NA + NB + NdCM) * sizeof(float);
	
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

	cudaStatus = cudaMalloc((void**)&dev_idxA, NdCM * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Malloc dev_cmA!\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_idxB, NdCM * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Malloc dev_cmB!\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_dCM, NdCM * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Malloc dev_cmB!\n");
		goto Error;
	}
	// initialize output
    cudaStatus = cudaMemset(dev_dist, 99.0, NA * NB * sizeof(double));
    if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMalloc", "Can't Memset dev_dist!\n");
		goto Error;
	}
	// transfer data
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
	cudaStatus = cudaMemcpy(dev_idxA, idxA, NdCM * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMemcpyHostToDevice", "Can't transfer cmA to device!\n");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_idxB, idxB, NdCM * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMemcpyHostToDevice", "Can't transfer cmB to device!\n");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_dCM, dCM, NdCM * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		mexErrMsgIdAndTxt("compdist:cudaMemcpyHostToDevice", "Can't transfer cmB to device!\n");
		goto Error;
	}


	// launch kernel	

	compdistKernel<<<griddims, blockdims, dyn_sharedmemsize>>>(dev_A, dev_B, dev_idxA, dev_idxB, dev_dCM, NA, NB, NdCM, dev_dist);

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

__global__ void compdistKernel(const double* A, const double* B, 
								const double* idxA, const double* idxB, const double* dCM,
								const unsigned int NA, const unsigned int NB, const unsigned int NdCM,
								double* dist) {

	// prepare shared memory
	extern __shared__ float ins[]; // these memory blocks will be reused for different purposes.
	float *vA = ins; 
	float *vB = (float*) &vA[NdCM];

	// load indices of A and B
	int N_steps = (NdCM - 1) / blockDim.x + 1;// for large vectors, we load in steps
	for (int n = 0; n < N_steps; n++) {
		int idx_v = n*blockDim.x + threadIdx.x;
		if (idx_v < NdCM) {
			vA[idx_v] = __double2float_rn(idxA[idx_v]) - 1;
		}
	}
	__syncthreads();
	for (int n = 0; n < N_steps; n++) {
		int idx_v = n*blockDim.x + threadIdx.x;
		if (idx_v < NdCM) {
			vB[idx_v] = __double2float_rn(idxB[idx_v]) - 1;
		}
	}
	__syncthreads();

	// load values of A and B
	for (int n = 0; n < N_steps; n++) {
		int idx_v = n*blockDim.x + threadIdx.x;
		if (idx_v < NdCM) {
			vA[idx_v] = __double2float_rn(A[(int) (blockIdx.x*NA + vA[idx_v])]);
		}
	}
	__syncthreads();
	for (int n = 0; n < N_steps; n++) {
		int idx_v = n*blockDim.x + threadIdx.x;
		if (idx_v < NdCM) {
			vB[idx_v] = __double2float_rn(B[(int) (blockIdx.y*NB + vB[idx_v])]);
			// compute difference
			vA[idx_v] = powf(vA[idx_v] - vB[idx_v],2);
		}
	}
	__syncthreads();

	// load dCM and calculate integral
	for (int n = 0; n < N_steps; n++) {
		int idx_v = n*blockDim.x + threadIdx.x;
		if (idx_v < NdCM) {
			vA[idx_v] = vA[idx_v] * __double2float_rn(dCM[idx_v]);
		}
	}
	__syncthreads();

	// output final distance
	double dist_ij = (double) sum(vA, NdCM); 
	__syncthreads();
	if (threadIdx.x == 0) {
		dist[blockIdx.x*NB + blockIdx.y] = sqrtf(dist_ij);
	}
	__syncthreads();
}

// sum elements along x within a block
__device__ float sum(float* x, const int len) {
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