: "${PREFIX:=$HOME/gkylsoft}"
./configure CC=nvcc  ARCH_FLAGS="-mcpu=native" --prefix=$PREFIX --use-lua=yes
