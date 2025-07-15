#!/bin/bash
## mapfile format: barcode<TAB>celltype
bam="$1"
mapfile="$2"

if [ $# -ne 2 ]; then
    echo "Usage: $0 <bam_file> <barcode_celltype_file>"
    exit 1
fi

# extract cell type
types=$(cut -f2 $mapfile | sort | uniq)

for ct in $types; do
    # extract all barcodes for this cell type
    barcodes=$(awk -v type=$ct '$2==type {print $1}' $mapfile | paste -sd'|' -)

    # skip if no barcodes found
    [ -z "$barcodes" ] && continue

    echo "Splitting $ct ..."

    samtools view -h $bam | \
    awk -v bc="$barcodes" 'BEGIN{OFS="\t"} 
        $0 ~ /^@/ || $0 ~ "CB:Z:("bc")" {print $0}' | \
    samtools view -Sb - > ${ct}.bam

    samtools index ${ct}.bam
done