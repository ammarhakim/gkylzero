./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$HOME/gkylsoft-amd --lapack-inc=${HOME}/gkylsoft-amd/OpenBLAS/include --lapack-lib=${HOME}/gkylsoft-amd/OpenBLAS/lib/libopenblas.a --superlu-inc=${HOME}/gkylsoft-amd/superlu/include --superlu-lib=${HOME}/gkylsoft-amd/superlu/lib/libsuperlu.a;


