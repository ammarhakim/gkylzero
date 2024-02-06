#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT
ADAS_DIR="$GKYLSOFT/gkylzero/share/adas"

python ../data/adas/download_adas.py
echo "Converting ADAS data to numpy .."
python ../data/adas/adas_to_numpy.py
rm *.dat
cp *.npy $ADAS_DIR/.
mv *.npy ../data/adas/.
echo "ADAS data downloaded to $ADAS_DIR"
