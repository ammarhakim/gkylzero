module unload darshan
cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --download-adas=yes  --prefix=$PREFIX
