#.MF 2024/03/07: At the time we got this to work I had the following modules loaded
#.  1) craype-x86-milan                        7) cpe/23.12                   13) craype-accel-nvidia80
#.  2) libfabric/1.15.2.0                      8) gpu/1.0                     14) cray-mpich/8.1.28     (mpi)
#.  3) craype-network-ofi                      9) craype/2.7.30       (c)     15) cudatoolkit/12.0      (g)
#.  4) xpmem/2.6.2-2.5_2.38__gd067c3f.shasta  10) cray-dsmml/0.2.2            16) nccl/2.18.3-cu12
#.  5) gcc-native/12.3                        11) cray-libsci/23.12.5 (math)
#.  6) perftools-base/23.12.0                 12) PrgEnv-gnu/8.5.0    (cpe)
#.Most of these are loaded by default, so we just load some extra/key ones here.
module load PrgEnv-gnu/8.5.0
module load craype-accel-nvidia80
module load cray-mpich/8.1.28
module load cudatoolkit/12.0
module load nccl/2.18.3-cu12

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --cudamath-libdir=/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/12.0/lib64 --use-mpi=yes --mpi-inc=$CRAY_MPICH_DIR/include --mpi-lib=$CRAY_MPICH_DIR/lib --use-nccl=yes --nccl-inc=$NCCL_DIR/include --nccl-lib=$NCCL_DIR/lib --use-lua=yes;


