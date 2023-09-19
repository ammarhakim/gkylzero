module load cudatoolkit/12.0
module load openmpi/5.0.0rc12
module unload darshan
export GKYLSOFT=$HOME/perlmutter/gkeyll/code/gpu/gkylsoft
./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$GKYLSOFT --lapack-inc=$GKYLSOFT/OpenBLAS/include --lapack-lib=$GKYLSOFT/OpenBLAS/lib/libopenblas.a --superlu-inc=$GKYLSOFT/superlu/include --superlu-lib=$GKYLSOFT/superlu/lib/libsuperlu.a --cudamath-libdir=/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/12.0/lib64;


