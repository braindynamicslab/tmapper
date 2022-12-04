# GWTLBgpu
Adapted by Mengsen Zhang (Sep, 2019) from algorithm by Chowdhury, S. and Mémoli, F. (2018; https://github.com/samirchowdhury/GWnets.git). See also paper by Chowdhury, S. and Mémoli, F. The Gromov-Wasserstein distance between networks and stable network invariants. 2018. https://arxiv.org/abs/1808.04337

mexEMD code from https://github.com/gpeyre/2016-ICML-gromov-wasserstein

_GWTLBgpu is a Matlab/Octave + CUDA package for computing lower bounds on the Gromov-Wasserstein distance between networks._ Some functions requires a CUDA-enabled NVIDIA graphic card. The GPU codes were initially developed by Mengsen Zhang (Sep, 2019) with a Geforce GTX 1070 NVIDIA graphic card on a personal laptop.  

## different implementations can be used to accomandate different speed and scalability requirements
Five different implementations of the algorithm has been included in this package. Three of them run on CPU alone, which are slower but can take matrices of arbitrary size allowed by memory. 

Two of them require collaboration between CPU and GPU (i.e. heterogenous computing). They are much faster (~100x their CPU-only counterparts), but they do not admit very large matrices.

### `emd2RTB(A,B,mA,mB)`
This is the original version from Chowdhury, S. and Mémoli, F..

### `emd2RTB_simple(A,B,mA,mB)`
This is a vectorized version of the original to reduce the use of for-loop in MATLAB. It behaves identically to the original version.

### `emd2RTB_uni(A,B)`
This is a specialized version for networks whose node measures (mA, mB, arguments of the two version above) are uniformly distributed. It behaves identically to the original version when node measures are uniformly distributed.

### `emd2RTLB_hetero(A,B,mA,mB)`
This is a heterogeneous verion (CPU + GPU) of `emd2TLB_simple`. The GPU part use single-precision arithmetics, so that the results agree with the CPU version upto 8-th significant digits. It cannot handle the cases where the total number of nodes in the two networks exceed 3071. This may require compiling the CUDA code `compDist_merge.cu` for your computing platform.

### `emd2RTLB_unih(A,B)`
This is a heterogeneous version of `emd2TLB_uni`. The GPU part use single-precision arithmetics, so that the results agree with the CPU version upto 8-th significant digits. It is NOT designed to handle the matrices greater than 4000x4000. This may require compiling the CUDA code `compDist_uni.cu` for your computing platform.

## test and compare different implementations

run `test.m`

## compiling CUDA code
The CUDA codes `compDist_merge.cu` and `compDist_uni.cu` are actually mex functions. You can easily compile them in MATLAB

```
mexcuda compDist_merge.cu
mexcuda compDist_uni.cu
```
which requires you to have the correct version of [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) already installed. If you have not, MATLAB will give you an error when you use `mexcuda` and tell you the correct version it needs. Or, you can check the Toolkit version using `gpuDevice`. Note that different CUDA toolkit maybe required for different versions of MATLAB.

If you use it on Stanford Sherlock, you should first test it in an interactive GPU node by requesting:
```
srun -p gpu --gres gpu:1 --pty bash
```

Once you have the node to yourself, first load required modules, e.g.
```
ml load matlab
ml load cuda/10.0.130
```

Then start MATLAB
```
matlab -nodesktop
```
and use `mexcuda` as mentioned above. 

