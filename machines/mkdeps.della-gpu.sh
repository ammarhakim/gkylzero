module load cudatoolkit/12.0
cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --download-adas=yes --prefix=$PREFIX
