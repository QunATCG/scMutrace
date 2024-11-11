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


def target_extractsnp_fromBAM(BamPath, referencePath, outVcfPath, config_posList, min_bq, min_mq, dp_check, alt_check, cell_check, nm_check, min_cov, min_alt, min_cell, max_nm, cbtag):
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
    min_alt (int): Minimum alternative allele count.
    min_cell (int): Minimum cell count.
    max_nm (int): Maximum number of mismatches.

    Returns:
    None
    """
    print('Run task ' + str(os.getpid()) + '...')
    pysamData = pysam.AlignmentFile(BamPath, 'rb')
    # chr
    chr = str(config_posList[0])
    # start
    start = int(config_posList[1]) - 50
    # end
    end = int(config_posList[2]) + 50
    # ref
    refFasta = pysam.FastaFile(referencePath)
    # pileups
    bam_pileups = pysamData.pileup(contig = str(chr), start = int(start), end = int(end), max_depth = 100000000, min_base_quality = min_bq, min_mapping_quality = min_mq)
    count_num_iter = 0
    # out vcf
    outPath = os.path.join(outVcfPath, f"VCF_{chr}_{config_posList[1]}_{config_posList[1]}_random.vcf")
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
        if POS_PILEUP < start or POS_PILEUP > end:
            continue
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
            BASE_QUALITY_LIST = bam_pileup.get_query_qualities()
            # query bases/sequences at pileup column position; return list
            BASE_PILEUP_LIST = bam_pileup.get_query_sequences(mark_matches=True, add_indels=True)
            NEW_BASE_PILEUP_LIST, COUNTNUMALT = count_base_typeNum(BASE_PILEUP_LIST, REF_BASE, False)
            # Check ALT number
            if alt_check == True and COUNTNUMALT < int(min_alt):
                continue
            if DP >= 1000 and COUNTNUMALT <= 0.005 * DP:
                continue
            # list of reads (pysam.PileupRead) aligned to this column
            # Representation of a read aligned to a particular position in the reference sequence.
            # -------------------------------------------------
            # READS_THIS_COLUMN > PileupColumn.pileups > pysam.PileupRead
            # -------------------------------------------------
            READS_THIS_COLUMN = bam_pileup.pileups
            # READS number
            READS_NUM = len(READS_THIS_COLUMN)
            # Setting output
            CELLS = {}
            BASE_COUNTS_DIC = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, "D": 0, "I": 0}
            BASE_COUNTS_F_DIC = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, "D": 0, "I": 0}
            BASE_COUNTS_R_DIC = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, "D": 0, "I": 0}
            #BASE_QUALITIES_DIC = {"A": [], "T": [], "C": [], "G": [], "N": [], "D": [], "I": []}
            CELL_COUNTS_DIC = {"A": [], "T": [], "C": [], "G": [], "N": [], "D": [], "I": []}
            CELL_COUNTS_F_DIC = {"A": [], "T": [], "C": [], "G": [], "N": [], "D": [], "I": []}
            CELL_COUNTS_R_DIC = {"A": [], "T": [], "C": [], "G": [], "N": [], "D": [], "I": []}
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
                    # BASE QUALITY
                    BASE_PILEUP = BASE_QUALITY_LIST[i]
                    if CELLBARCODE not in CELLS.keys():
                        CELLS[CELLBARCODE] = {"A": [], "T": [], "C": [], "G": [], "N": [], "D": [], "I": []}
                    CELLS[CELLBARCODE][BASE_READ_ALIGNMENT].append(BASE_PILEUP)
                    # BASE COUNT DIC
                    BASE_COUNTS_DIC[BASE_READ_ALIGNMENT] += 1
                    # check strand
                    if READ_ALIGNMENT.is_reverse == True:
                        BASE_COUNTS_R_DIC[BASE_READ_ALIGNMENT] += 1
                        CELL_COUNTS_R_DIC[BASE_READ_ALIGNMENT].append(CELLBARCODE)
                    else:
                        BASE_COUNTS_F_DIC[BASE_READ_ALIGNMENT] += 1
                        CELL_COUNTS_F_DIC[BASE_READ_ALIGNMENT].append(CELLBARCODE)
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
            #CELL_NUM = len(set(CELLS))
            if cell_check == True and COUNTOTALCELL < int(min_cell):
                continue
            CELLUNIQ, CELLMUTCOUNTINFO = get_mut_bycell(CELL_COUNTS_DIC)
            CELLUNIQ_F, CELLMUTCOUNTINFO_F = get_mut_bycell(CELL_COUNTS_F_DIC)
            CELLUNIQ_R, CELLMUTCOUNTINFO_R = get_mut_bycell(CELL_COUNTS_R_DIC)
            CELLUNIQ_STR = ";".join(CELLUNIQ)
            CELLMUTCOUNTINFO_STR = ";".join(CELLMUTCOUNTINFO)
            if len(CELLUNIQ_F) != 0:
                CELLUNIQ_F_INDEX = ";".join([str(CELLUNIQ.index(i)) for i in CELLUNIQ_F])
                CELLMUTCOUNTINFO_F_STR = ";".join(CELLMUTCOUNTINFO_F)
            else:
                CELLUNIQ_F_INDEX = "None"
                CELLMUTCOUNTINFO_F_STR = "None"
            if len(CELLUNIQ_R) != 0:
                CELLUNIQ_R_INDEX = ";".join([str(CELLUNIQ.index(i)) for i in CELLUNIQ_R])
                CELLMUTCOUNTINFO_R_STR = ";".join(CELLMUTCOUNTINFO_R)
            else:
                CELLUNIQ_R_INDEX = "None"
                CELLMUTCOUNTINFO_R_STR = "None"
            CELL_BASE_QUALITES = [integrate_data_fromQualityDic(CELLS[i]) for i in CELLUNIQ]
            BASE_QUALITIES_INFO = ";".join(CELL_BASE_QUALITES)
            SQL_0 = f"{chr}\t{POS_PILEUP + 1}\t{REF_BASE}\t{DP_TOTAL}\t{COUNTINFOALT}\t{COUNTOTALALT}\t{COUNTOTALCELL}\t{COUNTINFOCELL}\t{COUNTOTALCELLALT}"
            SQL_1 = f"{CELLUNIQ_STR}\t{CELLMUTCOUNTINFO_STR}\t{BASE_QUALITIES_INFO}\t{CELLUNIQ_F_INDEX}\t{CELLMUTCOUNTINFO_F_STR}\t{CELLUNIQ_R_INDEX}\t{CELLMUTCOUNTINFO_R_STR}"
            SQL = f"{SQL_0}\t{SQL_1}"
            outVcfFile.writelines(SQL + "\n")
            outVcfFile.flush()
    refFasta.close()
    outVcfFile.close()
    print(f'count_mut_information ###---### sub task was done: {chr}_{start}_{end}')

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
    target_extractsnp_fromBAM(BamPath, referencePath, outDirPath, config_info.split("_"), 10, min_mq, True, False, False, False, 1, 1, 1, 5, cbtag)

if __name__=='__main__':
    starttime = timeit.default_timer()
    main()
    stoptime = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stoptime - starttime)) + ' seconds')