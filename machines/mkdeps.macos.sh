cd install-deps
: "${PREFIX:=/Users/mfrancis/Documents/gkeyll/code/g0comm/gkylsoft}"
./mkdeps.sh --build-openblas=no --build-superlu=yes --prefix=$PREFIX --build-openmpi=no --use-adas=no
