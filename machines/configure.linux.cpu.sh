: "${PREFIX:=$HOME/gkylsoft}"
: "${MPI_HOME:=$HOME/gkylsoft/openmpi}"
./configure CC=clang --prefix=$PREFIX --use-lua=yes
