# BHM-BCP
A Simple and Efficient GPU-Accelerated Spectral Scheme 
for the Block Copolymer Equation via the 
Biharmonic Modified Method

## Overview
We present the Biharmonic Modified Method (BHM) applied 
to the Block Copolymer (BCP) equation for the first time. 
The BHM scheme is semi-implicit, requiring only a single 
linear solve per time step with no auxiliary variables or 
system enlargement. The scheme is naturally GPU-compatible 
via FFT operations and scales efficiently to 3D simulations 
on consumer GPU hardware. All classical BCP morphologies 
including lamellae, spheres, and gyroids are reproduced 
in both 2D and 3D with full reproducibility.

## Requirements
- MATLAB with Parallel Computing Toolbox
- NVIDIA GPU (any CUDA-capable GPU)
