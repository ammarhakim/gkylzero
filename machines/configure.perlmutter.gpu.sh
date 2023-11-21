module load cudatoolkit/12.0
module load openmpi/5.0.0rc12
module unload darshan
: "${PREFIX:=$HOME/gkylsoft}"
./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --cudamath-libdir=/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/12.0/lib64;


