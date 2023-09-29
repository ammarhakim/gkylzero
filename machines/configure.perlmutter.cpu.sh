module unload darshan
export GKYLSOFT=$HOME/gkylsoft
./configure CC=cc --prefix=$GKYLSOFT --lapack-inc=$GKYLSOFT/OpenBLAS/include --lapack-lib=$GKYLSOFT/OpenBLAS/lib/libopenblas.a --superlu-inc=$GKYLSOFT/superlu/include --superlu-lib=$GKYLSOFT/superlu/lib/libsuperlu.a;


