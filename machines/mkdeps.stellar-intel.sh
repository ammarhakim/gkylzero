module load intel/2021.1.2 
module load openmpi/intel-2021.1/4.1.2 
cd install-deps
: "${PREFIX:=/scratch/gpfs/as6384/gkapp-cpu/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$PREFIX
