#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT
ADAS_DIR="$GKYLSOFT/greorg/share"

mkdir -p $ADAS_DIR/adas
#cd -
python ../gyrokinetic/data/adas/download_adas.py
echo "Converting ADAS data to numpy .."
python ../gyrokinetic/data/adas/adas_to_numpy.py
rm *.dat
cp *.npy $ADAS_DIR/adas/.
mv *.npy ../gyrokinetic/data/adas/.
echo "ADAS data downloaded to $ADAS_DIR/adas"
