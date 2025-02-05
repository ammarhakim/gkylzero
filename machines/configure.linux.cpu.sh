: "${PREFIX:=$HOME/gkylsoft}"
: "${MPI_HOME:=$HOME/gkylsoft/openmpi}"
./configure CC=clang --prefix=$PREFIX --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib --use-adas=yes --use-lua=yes
