module load cudatoolkit/11.7
module unload darshan
./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$HOME/gkylsoft --lapack-inc=${HOME}/gkylsoft/OpenBLAS/include --lapack-lib=${HOME}/gkylsoft/OpenBLAS/lib/libopenblas.a --superlu-inc=${HOME}/gkylsoft/superlu/include --superlu-lib=${HOME}/gkylsoft/superlu/lib/libsuperlu.a --cudamath-libdir=/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/11.7/lib64;


