import os
import re
import time
import pysam
import random
import timeit
import argparse
from py_utils import *
import multiprocessing as mp

# get argparse
def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to obtain a list of base and cell counts in scRNA bam file')
    parser.add_argument('--bam', type=str, default=1, help='BAM file to be analysed', required = True)
    parser.add_argument('--ref', type=str, default=1, help='Path to reference genome version. *fai must be available in the same directory as the reference genome file', required = True)
    parser.add_argument('--contig', type=str, help='Config. Chr_Start_End', required = True)
    parser.add_argument('--minmq', default = 255, help='mapping quality', required = False, type = int)
    parser.add_argument('--cbtag', default= "CB", help='cell barcode tag', required = False)
    parser.add_argument('--outdir', default = '.', help='Out folder', required = False)
    return (parser)

def main():
    # arguments
    parser = initialize_parser()
    args = parser.parse_args()
    BamPath = args.bam
    referencePath = args.ref
    config_info = args.contig
    min_mq = args.minmq
    outDirPath = args.outdir
    cbtag = args.cbtag
    make_dir_check(outDirPath)
    target_extractsnp_fromBAM(BamPath, referencePath, outDirPath, config_info.split("_"), 0, min_mq, True, False, False, False, 2, 2, 2, 5, cbtag, 50, "test")

if __name__=='__main__':
    starttime = timeit.default_timer()
    main()
    stoptime = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stoptime - starttime)) + ' seconds')