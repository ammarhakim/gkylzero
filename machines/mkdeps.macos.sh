cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=no --build-superlu=no --use-adas=yes --prefix=$PREFIX
