module load cudatoolkit/12.0
cd install-deps
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$HOME/gkylsoft-amd/
