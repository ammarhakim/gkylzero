module load cudatoolkit/12.4
module load openmpi/cuda-11.1/gcc/4.1.1
: "${PREFIX:=/scratch/gpfs/manaurer/gkeyll/code/g0-amd-hack/gkylsoft}"
./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib64 --use-nccl=yes --nccl-inc=/usr/include --nccl-lib=/usr/lib64 --use-cudss=yes;


