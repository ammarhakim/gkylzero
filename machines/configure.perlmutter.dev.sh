./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$HOME/gkylsoft --lapack-inc=${HOME}/gkylsoft/OpenBLAS/include --lapack-lib=${HOME}/gkylsoft/OpenBLAS/lib/libopenblas.a --superlu-inc=${HOME}/gkylsoft/superlu/include --superlu-lib=${HOME}/gkylsoft/superlu/lib/libsuperlu.a;


