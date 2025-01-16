#!/bin/bash

bamInput=$1
# Set a quality threshold for depth calculation
quality_threshold=20
# Set a depth threshold for filtering
depth_threshold=5

# Calculate depth of coverage with a specified mapping quality threshold
samtools depth -q ${quality_threshold} ${bamInput} > ${bamInput}.depth

# Filter positions where depth is greater than specified threshold
awk -v threshold=${depth_threshold} '$3 > threshold' ${bamInput}.depth > ${bamInput}.depth.gt

# Convert to BED format, keeping chromosome, start, end, and depth as fields
awk '{print $1, $2, $2+1}' ${bamInput}.depth.gt > ${bamInput}.depth.gt.bed

# Merge overlapping BED intervals
bedtools merge -i ${bamInput}.depth.gt.bed > ${bamInput}.depth.gt.bed.merge

rm ${bamInput}.depth.gt ${bamInput}.depth.gt.bed