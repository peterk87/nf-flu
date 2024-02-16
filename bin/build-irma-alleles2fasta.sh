#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status.
set -euo pipefail

# Enable debugging to see the commands being executed, can be turned off in production for cleaner logs
# set -x

SOURCE_FILE="irma-alleles2fasta.v"
OUTPUT_FILE="irma-alleles2fasta"

# Check if V is installed and in PATH
if ! command -v v &> /dev/null; then
    echo "V compiler not found. Please ensure it is installed and in your PATH."
    exit 1
fi

# Compile the V program in production mode
echo "Compiling $SOURCE_FILE..."
v -prod -cg -cflags '--static' "$SOURCE_FILE"

# Check if strip is installed and in PATH
if command -v strip &> /dev/null; then
    echo "Stripping symbols from $OUTPUT_FILE..."
    strip -s "$OUTPUT_FILE"
else
    echo "Warning: 'strip' not found. Skipping symbol stripping."
fi

# Check if UPX is installed and in PATH
if command -v upx &> /dev/null; then
    echo "Compressing $OUTPUT_FILE with UPX..."
    upx --best "$OUTPUT_FILE"
else
    echo "Warning: 'upx' not found. Skipping compression."
fi

echo "Build and optimization complete."
