module load cudatoolkit/12.0
module load openmpi/cuda-11.1/gcc/4.1.1
cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-luajit=yes --prefix=$PREFIX
