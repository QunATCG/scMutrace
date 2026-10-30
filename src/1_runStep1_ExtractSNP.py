# import
import os
import re
import pysam
import timeit
import argparse
from py_utils import *
import multiprocessing as mp
import subprocess

# get argparse
def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to obtain a list of base and cell counts in scRNA bam file')
    parser.add_argument('--bam', type=str, default='', help='BAM file to be analysed', required=True)
    parser.add_argument('--ref', type=str, default='', help='Path to reference genome version. *fai must be available in the same directory as the reference genome file', required=True)
    parser.add_argument('--cbtag', type=str, default='CB', help='cell barcode tag', required=True)
    parser.add_argument('--contig', type=str, help='Contig. Chr_Start_End', required=True)
    parser.add_argument('--minmq', type=int, default=255, help='min mapping quality [Default: 255]', required=False)
    parser.add_argument('--minbq', type=int, default=20, help='min base quality [[Default: 20]]', required=False)
    parser.add_argument('--mindp', type=int, default=3, help='min sequencing depth [Default: 3]', required=False)
    parser.add_argument('--minaltdp', type=int, default=2, help='min alteration read number [Default: 2]', required=False)
    parser.add_argument('--mincell', type=int, default=3, help='min cell number [Default: 3]', required=False)
    parser.add_argument('--minaltcell', type=int, default=2, help='Minimum cell count with alternative allele. [Default: 2]', required=False)
    parser.add_argument('--nmcheck', action='store_true', help='A boolean flag for number of mismatches checking')
    parser.add_argument('--maxnm', type=int, default=3, help='max mismatch number [Default: 3]', required=False)
    parser.add_argument('--id', help='Prefix used to name output file', required=False)
    parser.add_argument('--maxcheck', action='store_true', help='A boolean flag for best performance checking, it will use full cpus')
    parser.add_argument('--nprocs', default=1, help='Number of processes [Default: 1]', required=False, type=int)
    parser.add_argument('--ratio', default='0.005', help='cutoff to filter noisy mutations', required=False)
    parser.add_argument('--outdir', default='.', help='Out folder', required=False)
    return parser

def err_call_back(err):
    print(f'apply_async error: {str(err)}')

def finished_callback(result):
    print(
        f"[DONE] {result['region']} | "
        f"total_sites={result['total_sites']:,} | "
        f"pass_sites={result['final_variants']:,} | "
        f"log={result['log']}",
        flush=True,
    )

def merge_VCF_filesShell(outPath, VCFtmpPath, IDSample, chrTag):
    """
    Merges VCF files for a specific chromosome tag and sample ID, and removes the temporary files.

    Args:
        outPath (str): The output directory where the merged VCF file will be saved.
        VCFtmpPath (str): The directory containing the temporary VCF files to be merged.
        IDSample (str): The sample ID to be used in the output file name.
        chrTag (str): The chromosome tag to filter the VCF files to be merged.

    Raises:
        subprocess.CalledProcessError: If an error occurs while executing the shell commands.
    """
    filesTotal = os.listdir(VCFtmpPath)
    VCFtmpTotal = [i for i in filesTotal if re.match("VCF_" + str(chrTag) + "_", i) is not None]
    VCFtmpTotalPaths = [os.path.join(VCFtmpPath, i) for i in VCFtmpTotal]
    VCFtmpTotalPathsSQL = " ".join(VCFtmpTotalPaths)
    cmd_merge_vcf = f"cat {VCFtmpTotalPathsSQL} >> {os.path.join(outPath, f'VCF_{IDSample}_{chrTag}.vcf')}"
    cmd_rm_vcf = f"rm -f {VCFtmpTotalPathsSQL}"
    try:
        subprocess.run(cmd_merge_vcf, shell=True, check=True)
        subprocess.run(cmd_rm_vcf, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.output}")

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    bam_Path = args.bam
    ref_Path = args.ref
    config_info = args.contig
    cbtag = args.cbtag
    nameSplitBamTag = args.id
    minmq = args.minmq
    minbq = args.minbq
    min_cov = args.mindp
    min_altcov = args.minaltdp
    min_cell = args.mincell
    min_altcell = args.minaltcell
    max_nm = args.maxnm
    nmcheck = args.nmcheck
    maxcheck = args.maxcheck
    nmCheckStates = nmcheck
    COREUSE = args.nprocs
    ratioUser = float(args.ratio)
    outDirPath = args.outdir
    logDir = os.path.join(outDirPath, "logs")
    outDirPathVCFTmp = os.path.join(outDirPath, "tmpVCF")
    
    # clean tmpVCF if it exists
    if os.path.exists(outDirPathVCFTmp):
        cmdclean = f"rm -f {outDirPathVCFTmp}/*"
        try:
            subprocess.run(cmdclean, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while executing command: {e.cmd}")
            print(f"Return code: {e.returncode}")
            print(f"Output: {e.output}")
    
    VcfPOSOUT = os.path.join(outDirPath, "VCFPOS")
    
    # make dir
    make_dir_check(outDirPath)
    make_dir_check(outDirPathVCFTmp)
    make_dir_check(logDir)
    make_dir_check(VcfPOSOUT)
    
    print("Processing " + bam_Path)
    
    # try to use max cores
    if maxcheck:
        if COREUSE != 1 and mp.cpu_count() > COREUSE:
            COREUSE = mp.cpu_count() - 1
            print("Try to use max cores number...")
        if COREUSE > mp.cpu_count():
            COREUSE = mp.cpu_count() - 1
            print("Cores is over computer max cpus...")
    else:
        if COREUSE > mp.cpu_count():
            COREUSE = mp.cpu_count() - 1
            print("Cores is over computer max cpus...")
    
    K = batch_config_window(config_info, groupNum=COREUSE)
    
    if K:
        if COREUSE > 1:
            print('Parent process %s.' % os.getpid())
            p = mp.get_context("fork").Pool(COREUSE)
            print("Your cores is " + str(COREUSE))
            for config_i in K:
                p.apply_async(count_mut_information, args=(bam_Path, ref_Path, outDirPathVCFTmp, config_i, minbq, minmq, True, True, True, nmCheckStates, min_cov, min_altcov, min_cell, min_altcell, max_nm, cbtag, nameSplitBamTag, ratioUser,logDir,), callback=finished_callback, error_callback=err_call_back)
            print('Waiting for all subprocesses done...')
            p.close()
            p.join()
            print('All subprocesses done.')
        else:
            print("Your cores is 1")
            for config_i in K:
                print(config_i)
                count_mut_information(BamPath=bam_Path, referencePath=ref_Path, outVcfPath=outDirPathVCFTmp, config_posList=config_i, min_bq=minbq, min_mq=minmq, dp_check=True, alt_check=True, cell_check=True, nm_check=nmCheckStates, min_cov=min_cov, min_altcov=min_altcov, min_cell=min_cell, min_altcell=min_altcell, max_nm=max_nm, cbtag=cbtag, nameSplitBamTag = nameSplitBamTag, ratioUser=ratioUser, logDir = logDir)
        
        print("Merging... " + bam_Path + " of Chromosome: " + str(config_info.split("_")[0]))
        merge_VCF_filesShell(VcfPOSOUT, outDirPathVCFTmp, nameSplitBamTag, str(config_info.split("_")[0]))
    else:
        print("No config info")
        print("Please check your config info")
        print("Exit")
        exit(0)

if __name__=='__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')