#!/bin/bash

# Improved script for calculating depth of coverage

# Usage information
usage() {
    echo "Usage: $0 <input.bam> [quality_threshold] [depth_threshold_gt5] [depth_threshold_gt8] [depth_threshold_gt10]"
    exit 1
}

# Check if input BAM file is provided
if [ -z "$1" ]; then
    usage
fi

# Assign input BAM file
bamInput=$1

# Set thresholds (use defaults if not provided)
quality_threshold=${2:-20}
depth_threshold_gt5=${3:-5}
depth_threshold_gt8=${4:-8}
depth_threshold_gt10=${5:-10}

# Check if input BAM file exists
if [ ! -f "$bamInput" ]; then
    echo "Error: File '$bamInput' not found!"
    exit 1
fi

# Check if required tools are installed
for tool in samtools bedtools; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: $tool is not installed."
        exit 1
    fi
done

echo "Calculating depth of coverage for $bamInput"
echo "Quality threshold: $quality_threshold"
echo "Depth thresholds: gt5=$depth_threshold_gt5, gt8=$depth_threshold_gt8, gt10=$depth_threshold_gt10"

# Create temporary files
depth_file="${bamInput}.depth"
depth_gt5_file=$(mktemp)
depth_gt8_file=$(mktemp)
depth_gt10_file=$(mktemp)
depth_gt5_bed_file=$(mktemp)
depth_gt8_bed_file=$(mktemp)
depth_gt10_bed_file=$(mktemp)

# Cleanup function to remove temporary files
cleanup() {
    rm -f "$depth_gt5_file" "$depth_gt8_file" "$depth_gt10_file"
    rm -f "$depth_gt5_bed_file" "$depth_gt8_bed_file" "$depth_gt10_bed_file"
}
trap cleanup EXIT

# Calculate depth of coverage
samtools depth -q "$quality_threshold" "$bamInput" > "$depth_file"
if [ $? -ne 0 ]; then
    echo "Error: samtools depth command failed."
    exit 1
fi

# Filter positions based on depth thresholds
awk -v threshold="$depth_threshold_gt5" '$3 >= threshold' "$depth_file" > "$depth_gt5_file"
awk -v threshold="$depth_threshold_gt8" '$3 >= threshold' "$depth_gt5_file" > "$depth_gt8_file"
awk -v threshold="$depth_threshold_gt10" '$3 >= threshold' "$depth_gt8_file" > "$depth_gt10_file"

# Convert filtered positions to BED format
awk -v OFS='\t' '{print $1, $2, $2}' "$depth_gt5_file" > "$depth_gt5_bed_file"
awk -v OFS='\t' '{print $1, $2, $2}' "$depth_gt8_file" > "$depth_gt8_bed_file"
awk -v OFS='\t' '{print $1, $2, $2}' "$depth_gt10_file" > "$depth_gt10_bed_file"

# Merge overlapping BED intervals
bedtools merge -i "$depth_gt5_bed_file" > "${bamInput}.depth.gt5.bed.merge"
if [ $? -ne 0 ]; then
    echo "Error: bedtools merge failed for depth > $depth_threshold_gt5."
    exit 1
fi

bedtools merge -i "$depth_gt8_bed_file" > "${bamInput}.depth.gt8.bed.merge"
if [ $? -ne 0 ]; then
    echo "Error: bedtools merge failed for depth > $depth_threshold_gt8."
    exit 1
fi

bedtools merge -i "$depth_gt10_bed_file" > "${bamInput}.depth.gt10.bed.merge"
if [ $? -ne 0 ]; then
    echo "Error: bedtools merge failed for depth > $depth_threshold_gt10."
    exit 1
fi

echo "Depth of coverage calculation completed successfully."
echo "Output files:"
echo "  $depth_file"
echo "  ${bamInput}.depth.gt5.bed.merge"
echo "  ${bamInput}.depth.gt8.bed.merge"
echo "  ${bamInput}.depth.gt10.bed.merge"
