: "${PREFIX:=$HOME/gkylsoft}"
./configure CC=cc  ARCH_FLAGS="-mcpu=native" --prefix=$PREFIX --use-lua=yes
