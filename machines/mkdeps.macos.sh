cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
echo 'here 0'
echo $PREFIX
echo 'here 1'
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$PREFIX
