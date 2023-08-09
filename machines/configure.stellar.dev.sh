export GKYLSOFT=$HOME/g0g2gpu/gkylsoft-amd
./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$HOME/g0g2gpu/gkylsoft-amd --lapack-inc=$GKYLSOFT/OpenBLAS/include --lapack-lib=$GKYLSOFT/OpenBLAS/lib/libopenblas.a --superlu-inc=$GKYLSOFT/superlu/include --superlu-lib=$GKYLSOFT/superlu/lib/libsuperlu.a;


