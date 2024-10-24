import os
import argparse
import pandas as pd

# Function to generate gCellSNP matrix
def initialize_parser():
    parser = argparse.ArgumentParser(description='generate gCellSNP matrix')
    parser.add_argument('--outdir', type=str, default='', help='out file path', required = True)
    parser.add_argument('--vcfpath', type=str, default='', help='Vcf path', required = True)
    return (parser)

# generate gCellSNP matrix
def generate_gCellSNP_matrix(outdir, vcfpath):
    cellSNP = {}
    # Read vcf file
    with open(vcfpath, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            if line[0] == "chr":
                continue
            