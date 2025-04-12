#!/bin/bash

# --- Configuration ---
n=$1
OUTNAME="output_$2"
FILENAME="$OUTNAME/output00000120.xml"
CONFIGNAME="config/PhysiCell_settings.xml.$2"
OUTDIR="$2"
# --- Configuration End ---

echo "-----------------------------------------------------"
echo "Starting script for $n runs."
echo "Will check for file: '$FILENAME' in the current directory."
echo "-----------------------------------------------------"

# --- Main Loop ---
for ((i=1; i<=n; i++)); do

  echo "[Run $i/$n] Starting run..."

  # --- Inner Retry Loop ---
  # Keep checking inside this while loop until the file exists
  while [[ ! -f "$FILENAME" ]]; do
    rm -Rf $OUTNAME/*
    ./pred_prey $CONFIGNAME
    echo "[Run $i/$n] File '$FILENAME' not found. retrying..."
    sleep 2
  done
  # --- End Inner Retry Loop ---

  # If the script exits the 'while' loop, the file was found.
  echo "[Run $i/$n] File '$FILENAME' found! Proceeding."
  echo "-----------------------------------------------------"
  python calculate_stats_mcds.py $OUTNAME
  mkdir -p "$OUTDIR/run$i"
  cp -a $OUTNAME "$OUTDIR/run$i"/.
  rm -Rf $OUTNAME/*

done
# --- End Main Loop ---

echo "Script finished after completing $n runs."
exit 0
