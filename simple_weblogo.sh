#!/bin/sh

FIRST_INDEX="$1"
YSCALE="$2"
INPUT_FILENAME="$3"
OUTPUT_FILENAME="$4"

CWD=$(pwd)
OUTPUT_DIR=$CWD/weblogo

mkdir -p $OUTPUT_DIR
#echo $INPUT_FILE
#echo $OUTPUT_DIR

weblogo \
  --fin $INPUT_FILENAME \
  --datatype fasta \
  --fout $OUTPUT_DIR/$OUTPUT_FILENAME \
  --format png \
  --sequence-type protein \
  --units bits \
  --first-index $FIRST_INDEX \
  --size large \
  --yaxis $YSCALE \
  --resolution 600 \
  --errorbars NO \
  --scale-width NO \
  --fineprint ''
