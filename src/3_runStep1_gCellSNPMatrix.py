import os
import timeit
import argparse
import pandas as pd

# Function to generate gCellSNP matrix
def initialize_parser():
    parser = argparse.ArgumentParser(description='generate gCellSNP matrix')
    parser.add_argument('--outpath', type=str, default='', help='out file path', required = True)
    parser.add_argument('--vcfpath', type=str, default='', help='Vcf path', required = True)
    return (parser)

# get cell info
def get_cell_info(cellMutInfo):
    cellsInfo = cellMutInfo.split(";")[1]
    if cellMutInfo != "":
        cells = [cell.split("-")[0] for cell in cellsInfo.split(":")]
    else:
        cells = []   
    return cells

# generate gCellSNP matrix
def generate_gCellSNP_matrix(outpath, vcfpath):
    cellSNP = {}
    # Read vcf file
    with open(vcfpath, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            line_len = len(line)
            if line[0] == "chr":
                continue
            loctmp = line[0] + "_" + line[1]
            celltmp = []
            for i in range(14, line_len):
                cellMutInfo = line[i]
                cells = get_cell_info(cellMutInfo)
                celltmp.extend(cells)
            cellSNP[loctmp] = celltmp
    # Generate gCellSNP matrix
    cellSNP_matrix = pd.DataFrame(0, index=cellSNP.keys(), columns=list(set([item for sublist in cellSNP.values() for item in sublist])))
    for loc, cells in cellSNP.items():
        cellSNP_matrix.loc[loc, cells] = 1
    cellSNP_matrix.to_csv(outpath)
    
# Main function
def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outpath = args.outpath
    vcfpath = args.vcfpath
    generate_gCellSNP_matrix(outpath, vcfpath)

if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    end = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(end - start)) + ' seconds')