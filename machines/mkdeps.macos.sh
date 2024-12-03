cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=no --build-superlu=yes --prefix=$PREFIX --build-openmpi=no --use-adas=no
