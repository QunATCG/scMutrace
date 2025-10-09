#!/bin/bash
# Script to check integrity of gzipped FASTQ files
# Usage: ./checkFQ.sh

for file in *.gz
do
    # Check if the file exists to prevent empty directory cases
    if [[ ! -f "$file" ]]; then
        echo "No .gz files found."
        break
    fi

    echo "Checking $file ..."

    # Check gzip integrity
    gunzip -t "$file"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: $file failed gunzip integrity test"
        continue
    fi

    # Count lines
    lines=$(zcat "$file" | wc -l)
    if [[ $((lines % 4)) -ne 0 ]]; then
        echo "ERROR: $file line count $lines not a multiple of 4"
    else
        echo "$file OK: $lines lines"
    fi
done
