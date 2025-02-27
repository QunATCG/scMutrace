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
    parser.add_argument('--nmcheck', action='store_true', help='A boolean flag for number of mismatches checking')
    parser.add_argument('--maxnm', type=int, default=3, help='max mismatch number [Default: 3]', required=False)
    parser.add_argument('--id', help='Prefix used to name output file', required=False)
    parser.add_argument('--maxcheck', action='store_true', help='A boolean flag for best performance checking, it will use full cpus')
    parser.add_argument('--nprocs', default=1, help='Number of processes [Default: 1]', required=False, type=int)
    parser.add_argument('--outdir', default='.', help='Out folder', required=False)
    return parser

def err_call_back(err):
    print(f'apply_async error: {str(err)}')

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

def count_mut_information_getPosCNDepth(BamPath, referencePath, outVcfPath, config_posList, min_bq, min_mq, dp_check, nm_check, min_cov, max_nm, cbtag):
    """
    Count mutation information from BAM file and write to VCF file.

    Parameters:
    BamPath (str): Path to the BAM file.
    referencePath (str): Path to the reference genome file.
    outVcfPath (str): Output directory for the VCF file.
    config_posList (list): List containing chromosome, start, and end positions.
    min_bq (int): Minimum base quality.
    min_mq (int): Minimum mapping quality.
    dp_check (bool): Flag to check depth coverage.
    alt_check (bool): Flag to check alternative allele count.
    cell_check (bool): Flag to check cell count.
    nm_check (bool): Flag to check number of mismatches.
    min_cov (int): Minimum coverage.
    min_altcov (int): Minimum alternative allele count.
    min_cell (int): Minimum cell count.
    min_altcell (int): Minimum cell count with alternative allele.
    max_nm (int): Maximum number of mismatches.

    Returns:
    None
    """
    print('Run task ' + str(os.getpid()) + '...')
    pysamData = pysam.AlignmentFile(BamPath, 'rb')
    # chr
    chr = str(config_posList[0])
    # start
    start = int(config_posList[1])
    # end
    end = int(config_posList[2])
    # ref
    refFasta = pysam.FastaFile(referencePath)
    # pileups
    bam_pileups = pysamData.pileup(contig = str(chr), start = int(start), end = int(end), max_depth = 100000000, min_base_quality = min_bq, min_mapping_quality = min_mq)
    count_num_iter = 0
    # out vcf
    outPath = os.path.join(outVcfPath, f"VCF_{chr}_{start}_{end}.vcf")
    outVcfFile = open(outPath, "w+")
    for bam_pileup in bam_pileups:
        count_num_iter += 1
        if count_num_iter % 500000 == 0:
            print("Processed " + str(count_num_iter) + " reads")
        # PileupColumn
        # the position in the reference sequence (0-based)
        # -------------------------------------------------
        # bam_pileup > pysam.PileupColumn
        # -------------------------------------------------
        POS_PILEUP = bam_pileup.reference_pos
        # fetch sequences in a region
        REF_BASE = refFasta.fetch(chr, POS_PILEUP, POS_PILEUP + 1).upper()
        # return number of aligned bases at pileup column position. this method applies a base quality filter; return int
        DP = bam_pileup.get_num_aligned()
        # Check ref base
        if REF_BASE.upper() != "N":
            # Check DP
            if dp_check == True and DP < int(min_cov):
                continue
            # query mapping quality scores at pileup column position; return list
            # MAPPING_QUALITY_LIST = bam_pileup.get_mapping_qualities()
            # query/read names aligned at pileup column position; return list
            # READNAMES_LIST = bam_pileup.get_query_names()
            # positions in read at pileup column position; return list
            # POS_LIST = bam_pileup.get_query_positions()
            # query base quality scores at pileup column position; return list
            # query bases/sequences at pileup column position; return list
            BASE_PILEUP_LIST = bam_pileup.get_query_sequences(mark_matches=True, add_indels=True)
            NEW_BASE_PILEUP_LIST, COUNTNUMALT = count_base_typeNum(BASE_PILEUP_LIST, REF_BASE, False)
            # list of reads (pysam.PileupRead) aligned to this column
            # Representation of a read aligned to a particular position in the reference sequence.
            # -------------------------------------------------
            # READS_THIS_COLUMN > PileupColumn.pileups > pysam.PileupRead
            # -------------------------------------------------
            READS_THIS_COLUMN = bam_pileup.pileups
            # READS number
            READS_NUM = len(READS_THIS_COLUMN)
            # Setting output
            BASE_COUNTS_DIC = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, "D": 0, "I": 0}
            #BASE_QUALITIES_DIC = {"A": [], "T": [], "C": [], "G": [], "N": [], "D": [], "I": []}
            CELL_COUNTS_DIC = {"A": [], "T": [], "C": [], "G": [], "N": [], "D": [], "I": []}
            # Check each reads in this column
            for i in range(0, READS_NUM):
                # a pysam.AlignedSegment object of the aligned read
                # -------------------------------------------------
                # READ_ALIGNMENT > pysam.AlignedSegment
                # -------------------------------------------------
                READ_ALIGNMENT = READS_THIS_COLUMN[i].alignment
                # Check read quality
                # print(READ_ALIGNMENT.is_duplicate, READ_ALIGNMENT.is_secondary, READ_ALIGNMENT.is_supplementary, READ_ALIGNMENT.is_mapped)
                if READ_ALIGNMENT.is_duplicate == False and  READ_ALIGNMENT.is_secondary == False and READ_ALIGNMENT.is_supplementary == False and READ_ALIGNMENT.is_mapped == True:
                    # CELLS
                    tagCheck = READ_ALIGNMENT.has_tag(cbtag)
                    if tagCheck:
                        CELLBARCODE = READ_ALIGNMENT.opt(cbtag).split("-")[0]
                    else:
                        continue
                    # num_mismatch
                    if nm_check:
                        tagCheck_NM = READ_ALIGNMENT.has_tag("NM")
                        if tagCheck_NM:
                            NM_count = READ_ALIGNMENT.opt("NM")
                            if NM_count > max_nm:
                                continue
                        else:
                            tagCheck_NM = READ_ALIGNMENT.has_tag("nM")
                            if tagCheck_NM == False:
                                print("No NM/nM tag in your bam, please set nm_check as False")
                            if tagCheck_NM:
                                NM_count = READ_ALIGNMENT.opt("nM")
                                if NM_count > max_nm:
                                    continue
                    # BASE
                    BASE_READ_ALIGNMENT = NEW_BASE_PILEUP_LIST[i]
                    if BASE_READ_ALIGNMENT == "NA" or BASE_READ_ALIGNMENT == "O":
                        continue
                    # BASE COUNT DIC
                    BASE_COUNTS_DIC[BASE_READ_ALIGNMENT] += 1
                    # CELL COUNT DIC
                    CELL_COUNTS_DIC[BASE_READ_ALIGNMENT].append(CELLBARCODE)
                    #CELLS.append(CELLBARCODE)
                    # ALTERNATIVE BASE
            # get real DP of cells with barcodes
            DP_TOTAL = count_base_fromDic(BASE_COUNTS_DIC, ["N"])
            if dp_check:
                if int(DP_TOTAL) < int(min_cov):
                    continue
            # CELL NUMBER
            # CELL MUTINFO Total
            COUNTINFOCELL, COUNTINFOALT, COUNTOTALCELL, COUNTOTALALT, COUNTOTALCELLALT = get_mut_total(CELL_COUNTS_DIC, REF_BASE)
            SQL_0 = f"{chr}\t{POS_PILEUP + 1}\t{REF_BASE}\t{DP_TOTAL}\t{COUNTINFOALT}\t{COUNTOTALALT}\t{COUNTOTALCELL}\t{COUNTINFOCELL}\t{COUNTOTALCELLALT}"
            outVcfFile.writelines(SQL_0 + "\n")
            outVcfFile.flush()
    refFasta.close()
    outVcfFile.close()
    print(f'count_mut_information ###---### sub task was done: {chr}_{start}_{end}')

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
    max_nm = args.maxnm
    nmcheck = args.nmcheck
    maxcheck = args.maxcheck
    nmCheckStates = nmcheck
    COREUSE = args.nprocs
    outDirPath = args.outdir
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
                p.apply_async(count_mut_information_getPosCNDepth, args=(bam_Path, ref_Path, outDirPathVCFTmp, config_i, minbq, minmq, True, nmCheckStates, min_cov, max_nm, cbtag), error_callback=err_call_back)
            print('Waiting for all subprocesses done...')
            p.close()
            p.join()
            print('All subprocesses done.')
        else:
            print("Your cores is 1")
            for config_i in K:
                print(config_i)
                count_mut_information_getPosCNDepth(BamPath=bam_Path, referencePath=ref_Path, outVcfPath=outDirPathVCFTmp, config_posList=config_i, min_bq=minbq, min_mq=minmq, dp_check=True, nm_check=nmCheckStates, min_cov=min_cov, max_nm=max_nm, cbtag=cbtag)
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