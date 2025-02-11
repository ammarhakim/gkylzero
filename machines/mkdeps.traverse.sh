module load cudatoolkit/12.0
module load nccl/cuda-11.7/2.14.3
module load openmpi/cuda-11.7/nvhpc-22.5/4.1.4/64

cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$PREFIX --build-luajit=yes
