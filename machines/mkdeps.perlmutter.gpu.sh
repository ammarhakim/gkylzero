module load cudatoolkit/12.0
module load openmpi/5.0.0rc12
module load nccl/2.18.3-cu12
module unload darshan
cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$PREFIX --build-openmpi=no MPICC=mpicc  MPICXX=mpicxx
