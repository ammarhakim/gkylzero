module load intel/2021.1.2
module load openmpi/intel-2021.1/4.1.0 
cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-luajit=yes --prefix=$PREFIX
