cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=no --build-superlu=yes --build-luajit=yes --prefix=$PREFIX
