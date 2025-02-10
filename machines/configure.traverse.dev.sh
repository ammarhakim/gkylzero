module load cudatoolkit/12.0
module load nccl/cuda-11.7/2.14.3
module load openmpi/cuda-11.7/nvhpc-22.5/4.1.4/64

# Comment out lines 7 and 29 of zero/gkyl_dflt.h because gkyl doesn't have the
# architecture intrinsics for IBM POWER9

: "${PREFIX:=$HOME/gkylsoft}"
: "${MPI_SOURCE:=/usr/local/openmpi/cuda-11.7/4.1.4/nvhpc225/ppc64le}"
: "${NCCL_HOME:=$HOROVOD_NCCL_HOME}"
./configure CC=nvcc  ARCH_FLAGS="-mcpu=native" --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --use-mpi=yes --mpi-inc=$MPI_SOURCE/include --mpi-lib=$MPI_SOURCE/lib64 --use-nccl=yes --nccl-inc=$NCCL_HOME/include --nccl-lib=$NCCL_HOME/lib64 --use-lua=yes;
