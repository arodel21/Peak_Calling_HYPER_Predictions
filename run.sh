#!/bin/bash

# ---------------------------------------------------------------------------- #
# Parameters
# ---------------------------------------------------------------------------- #

TRACK=<track_name>
SIGNAL=<path_to_bigwig>
OUT_DIR=<path_to_outdir>
BIN_SIZE=         #Leave empty for original resolution (128 bp)
CHROMS=           #Leave empty for all chromosomes
BACKGROUND=auto
MODEL=hybrid      #hybrid, bulk, or full
BULK_QUANTILE=auto
FOLD_CHANGE=5
THRESHOLD_MODE=density
DENSITY_THRESHOLD=1e-3
PVALUE=1e-2
LOC=auto

# ---------------------------------------------------------------------------- #

OUTPUT="${OUT_DIR}/${TRACK}_${MODEL}_lognorm_fc${FOLD_CHANGE}"

python call_peaks.py $SIGNAL \
    --output            "${OUTPUT}.bed" \
    --plot              "${OUTPUT}_fit.png" \
    --background        "$BACKGROUND" \
    --model             "$MODEL" \
    --bulk-quantile     "$BULK_QUANTILE" \
    --loc               "$LOC" \
    --fold-change       "$FOLD_CHANGE" \
    --threshold-mode    "$THRESHOLD_MODE" \
    --density-threshold "$DENSITY_THRESHOLD" \
    --pvalue            "$PVALUE" \
    --name-prefix       "$TRACK" \
    "${EXTRA[@]}"
