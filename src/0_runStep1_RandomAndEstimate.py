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
    parser.add_argument('--nprocs', default = 1, help='Number of processes [Default: 1]', required=False, type = int)
    parser.add_argument('--maxcheck', action='store_true', help='A boolean flag for best performance checking, it will use full cpus')
    parser.add_argument('--outdir', default = '.', help='Out folder', required = False)
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

def random_extractsnp_fromBAM(BamPath, referencePath, outVcfFile, config_posList, min_bq, min_mq, dp_check, min_cov, maxloc):
    #start = time.time()
    #print('Run task ' + str(os.getpid()) + '...')
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
    #count_num_iter = 0
    # out vcf
    #outPath = outVcfPath + "VCF_" + str(chr) + "_" + str(start) + "_" + str(end) + ".vcf"
    #outVcfFile = open(outPath, "w+")
    # out Pos index
    #outPosIndexPath = outVcfPath + "Pos_" + str(chr) + "_" + str(start) + "_" + str(end) + ".Index"
    #outPosIndexFile = open(outPosIndexPath, "w+")
    for bam_pileup in bam_pileups:
        # PileupColumn
        # the position in the reference sequence (0-based)
        # -------------------------------------------------
        # bam_pileup > pysam.PileupColumn
        # -------------------------------------------------
        POS_PILEUP = bam_pileup.reference_pos
        if POS_PILEUP != start:
            continue
        if POS_PILEUP >= maxloc:
            print("over max selected pos")
            break
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
                    tagCheck = READ_ALIGNMENT.has_tag("CB")
                    if tagCheck:
                        CELLBARCODE = READ_ALIGNMENT.opt("CB").split("-")[0]
                    else:
                        continue
                    # BASE
                    BASE_READ_ALIGNMENT = NEW_BASE_PILEUP_LIST[i]
                    if BASE_READ_ALIGNMENT == "NA" or BASE_READ_ALIGNMENT == "O":
                        continue
                    # BASE COUNT DIC
                    BASE_COUNTS_DIC[BASE_READ_ALIGNMENT] += 1
                    CELL_COUNTS_DIC[BASE_READ_ALIGNMENT].append(CELLBARCODE)
                    # ALTERNATIVE BASE
            # get real DP of cells with barcodes
            DP_TOTAL = count_base_fromDic(BASE_COUNTS_DIC, ["N"])
            if dp_check == True and DP_TOTAL < int(min_cov):
                continue
            # CELL NUMBER
            # CELL MUTINFO Total
            COUNTINFOCELL, COUNTINFOALT, COUNTOTALCELL, COUNTOTALALT, COUNTOTALCELLALT = get_mut_total(CELL_COUNTS_DIC, REF_BASE)
            SQL_0 = str(chr) + "\t" + str(POS_PILEUP + 1) + "\t" + str(REF_BASE) + "\t" + str(DP_TOTAL) + "\t" + COUNTINFOALT + "\t" + str(COUNTOTALALT) + "\t" + str(COUNTOTALCELL) + "\t" + COUNTINFOCELL + "\t" + str(COUNTOTALCELLALT)
            #SQL_1 = CELLUNIQ_STR + "\t" + CELLMUTCOUNTINFO_STR  + "\t"+ str(BASE_QUALITIES_INFO) + "\t" + CELLUNIQ_F_INDEX + "\t" + CELLMUTCOUNTINFO_F_STR + "\t" + CELLUNIQ_R_INDEX + "\t" + CELLMUTCOUNTINFO_R_STR
            SQL = SQL_0
            #POS_INDEX = str(chr) + "_" + str(POS_PILEUP + 1) + "_" + str(REF_BASE)
            #print(SQL)
            outVcfFile.writelines(SQL + "\n")
            outVcfFile.flush()
            #outPosIndexFile.writelines(POS_INDEX + "\n")
    refFasta.close()
    pysamData.close()
    #outVcfFile.close()
    #outPosIndexFile.close()
    end = time.time()
    #print('count_mut_information ###---### sub task was done.')

def iter_random_process(locList, BamPath, referencePath, outVcfPath, chrContig, min_bq, min_mq, dp_check, min_cov):
    outPath = outVcfPath + "VCF_" + str(chrContig) + "_" + str(locList[0]) + "_random.vcf"
    NumCheck = int(len(locList)/10)
    maxloc = max(locList)
    count_num_iter = 0
    outVcfFile = open(outPath, "w")
    for loc in locList:
        count_num_iter += 1
        if count_num_iter % NumCheck == 0:
            print("Processed " + str(count_num_iter) + " locs")
        config_posList = [chrContig, loc, loc + 1]
        random_extractsnp_fromBAM(BamPath, referencePath, outVcfFile, config_posList, min_bq, min_mq, dp_check, min_cov, maxloc)

def main():
    # arguments
    parser = initialize_parser()
    args = parser.parse_args()
    bam_Path = args.bam
    ref_Path = args.ref
    config_info = args.contig
    COREUSE = args.nprocs
    mqmin = args.minmq
    outDirPath = args.outdir
    make_dir_check(outDirPath)
    numSitesUsedForEst = 1000000
    chr = config_info.split("_")[0]
    #start = config_info.split("_")[1]
    end = config_info.split("_")[2]
    maxcheck = args.maxcheck
    numRandomList = random.sample(range(1, int(end)), numSitesUsedForEst)
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
    groupLength = int(numSitesUsedForEst/COREUSE)
    difNum = numSitesUsedForEst - COREUSE * groupLength
    numRandomListSub = []
    tmpCount = 0
    for i in range(0, COREUSE):
        tmpCount += 1
        if tmpCount < COREUSE - 1:
            subList = numRandomList[(i*groupLength):(i+1) * groupLength]
        else:
            subList = numRandomList[(i*groupLength):(i+1) * groupLength + difNum]
        numRandomListSub.append(subList)
    if COREUSE > 1:
        print('Parent process %s.' % os.getpid())
        print("Your cores is " + str(COREUSE))
        p = mp.get_context("fork").Pool(COREUSE)
        #p = Pool(COREUSE)
        print("Your cores is " + str(COREUSE))
        print(len(numRandomListSub))
        #countSubList = 0
        for subList in numRandomListSub:
            #countSubList += 1
            #outVcfFile = outDirPath + "VCF_" + str(chr) + "_" + str(countSubList) + "_random.vcf"
            p.apply_async(iter_random_process, args = (subList, bam_Path, ref_Path, outDirPath, chr, 20, mqmin, True, 3), error_callback = err_call_back)
        #    print(subList)
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')
    else:
        #countSubList = 0
        for subList in numRandomListSub:
            #outVcfFile = outDirPath + "VCF_" + str(chr) + "_" + str(countSubList) + "_random.vcf"
            iter_random_process(subList, bam_Path, ref_Path, outDirPath, chr, 20, 255, True, 3)
    randomTmpfiles = [i for i in os.listdir(outDirPath) if re.search("VCF_" + chr, i ) != None]
    randomTmpfilePaths = [outDirPath + i for i in randomTmpfiles]
    fileSql = " ".join(randomTmpfilePaths)
    cmd0 = "cat " + fileSql + " > " + outDirPath + "random_vcf_" + str(chr) + ".vcf"
    cmd1 = "rm " + fileSql
    os.system(cmd0)
    os.system(cmd1)

if __name__=='__main__':
    start = timeit.default_timer()
    #19_1_59128983
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')

# python 0_runStep1_RandomAndEstimate.py --bam /Users/liqun/BrainRNA/SNP/testData/scRNA/bam/P6_cSCC_scRNA.chr19.bam --ref /Users/liqun/reference/hg19/GRCh37/genome.fa --contig 19_1_59128983  --nprocs 10 --outdir /Users/liqun/BrainRNA/SNP/testData/scRNA/tmp/subam/tmpRandomVCF/