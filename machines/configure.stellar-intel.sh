module load intel/2021.1.2 
module load openmpi/intel-2021.1/4.1.2 

: "${PREFIX:=$HOME/gkylsoft}"
./configure CC=cc --prefix=$PREFIX --use-adas=yes
