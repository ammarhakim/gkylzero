module load cudatoolkit/11.7
#module unload darshan
./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$HOME/g0g2-gpu/gkylsoft --lapack-inc=${HOME}/g0g2-gpu/gkylsoft/OpenBLAS/include --lapack-lib=${HOME}/g0g2-gpu/gkylsoft/OpenBLAS/lib/libopenblas.a --superlu-inc=${HOME}/g0g2-gpu/gkylsoft/superlu/include --superlu-lib=${HOME}/g0g2-gpu/gkylsoft/superlu/lib/libsuperlu.a --cudamath-libdir=/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/11.7/lib64;


