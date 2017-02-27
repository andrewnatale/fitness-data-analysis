#!/bin/sh

#INPUT_FILE=$1
YSCALE=5
START=1

CWD=$(pwd)
OUTPUT_DIR=$CWD/pngs

mkdir -p $OUTPUT_DIR
#echo $INPUT_FILE
#echo $OUTPUT_DIR

for file in *.fasta; do
    OUTPUT_FILENAME=$file.png
    weblogo \
        --fin $file \
        --datatype fasta \
        --fout $OUTPUT_DIR/$OUTPUT_FILENAME \
        --format png \
        --sequence-type protein \
        --units bits \
        --first-index $START \
        --size large \
        --yaxis $YSCALE \
        --resolution 600 \
        --errorbars NO \
        --scale-width NO \
        --fineprint ''
done
