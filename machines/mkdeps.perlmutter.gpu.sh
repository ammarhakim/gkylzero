module load cudatoolkit/11.7
module unload darshan
cd install-deps
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$HOME/gkylsoft/
