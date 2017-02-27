#!/bin/sh

# quick score test

ROSETTA_BIN=/Users/anatale/Rosetta/main/source/bin
ROSETTA_EXE=score_jd2.default.macosclangrelease
ROSETTA_DB=/Users/anatale/Rosetta/main/database

PDB_IN=$1
PARAMS_IN=$2

$ROSETTA_BIN/$ROSETTA_EXE \
  -database $ROSETTA_DB \
  -s $PDB_IN \
  -extra_res_fa $PARAMS_IN \
  -out:no_nstruct_label \
  -out:file:score_only \
  -out:file:scorefile withparams.sc

$ROSETTA_BIN/$ROSETTA_EXE \
  -database $ROSETTA_DB \
  -s $PDB_IN \
  -ignore_unrecognized_res \
  -out:no_nstruct_label \
  -out:file:score_only \
  -out:file:scorefile withoutparams.sc
