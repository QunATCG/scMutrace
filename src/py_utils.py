####################
####################
####################
# def functions
####################
####################
####################
import os
import re
import math
import pysam
import random
import linecache
import subprocess
import collections
import numpy as np
from scipy import stats
from random import shuffle
from scipy.stats import betabinom
from scipy.optimize import minimize
from itertools import takewhile, repeat
from statsmodels.stats.proportion import proportions_ztest
####################
# make dir
def make_dir_check(dirPath):
    """
    Ensures that a directory exists at the specified path. If the directory does not exist, it is created.

    Args:
        dirPath (str): The path of the directory to check or create.

    Returns:
        None

    Side Effects:
        - Creates the directory at the specified path if it does not exist.
        - Prints a message indicating whether the directory was created or already existed.
    """
    dirPath = dirPath.strip()
    dirPath = dirPath.rstrip("\\")
    isExist = os.path.exists(dirPath)
    if not isExist:
        os.makedirs(dirPath, exist_ok=True)
    else:
        print(f"{dirPath} already exists")
    print(f"Directory created: {dirPath}")

# split file
def split_file_by_lines(inputfile, tag, lines_per_file, pathdir):
    """
    Splits a file into multiple smaller files, each containing a specified number of lines.

    Args:
        inputfile (str): The path to the input file that needs to be split.
        tag (str): A tag to be appended to the output file names.
        lines_per_file (int): The number of lines each split file should contain.
        pathdir (str): The directory where the split files will be saved.

    Returns:
        None
    """
    cmd = ["split", "-l", str(lines_per_file), inputfile, os.path.join(pathdir, "splitVcf_" + tag)]
    subprocess.run(cmd, check=True)

# split muti file
def split_multi_files(CORES, pathdir):
    """
    Splits multiple VCF files in a directory into smaller files for parallel processing.

    Args:
        CORES (int): The number of cores to use for parallel processing.
        pathdir (str): The directory containing the VCF files to be split.

    Returns:
        None

    Side Effects:
        - Removes previously split files in the directory.
        - Splits the VCF files into smaller files based on the number of cores.
    """
    splitfileCheck = [i for i in os.listdir(pathdir) if re.search("splitVcf", i) is not None]
    if splitfileCheck:
        print("Removing pre-split files...")
        splitfileCheckPaths = [os.path.join(pathdir, i) for i in splitfileCheck]
        for filePath in splitfileCheckPaths:
            os.remove(filePath)

    totalfilenames = [i for i in os.listdir(pathdir) if re.search("vcf", i) is not None]
    totalfilenamePaths = [os.path.join(pathdir, i) for i in totalfilenames]
    fileposnum = [get_Vcf_rownum(i) for i in totalfilenamePaths]
    totalpos = sum(fileposnum)
    binsize = int(totalpos / CORES)

    for filename, posnum in zip(totalfilenames, fileposnum):
        tmpfilePath = os.path.join(pathdir, filename)
        # if tmpfilePath size is zero then skip
        if posnum == 0:
            continue
        if posnum > binsize:
            split_file_by_lines(tmpfilePath, filename, binsize, pathdir)
        else:
            split_file_by_lines(tmpfilePath, filename, posnum, pathdir)

# batch window
# input example: 1_1_249250621
def batch_config_window(configItem, groupNum = None):
    """
    Splits a configuration item into a list of windows based on the specified group number.
    Args:
        configItem (str): A string in the format "chr_start_end" where 'chr' is the chromosome,
                            'start' is the starting position, and 'end' is the ending position.
        groupNum (int, optional): The number of groups to split the window into. If not provided,
                                    defaults to 5.
    Returns:
        list: A list of lists, where each sublist contains the chromosome, start, and end positions
                of each window.
    Example:
        >>> batch_config_window("chr1_0_1000", 5)
        [['chr1', 1, 200], ['chr1', 200, 400], ['chr1', 400, 600], ['chr1', 600, 800], ['chr1', 800, 1000]]
    """
    res = []
    chr = configItem.split("_")[0]
    end = int(configItem.split("_")[2])
    # Set default groupNum if not provided
    if groupNum is None:
        groupNum = 5
    # Determine window size based on groupNum
    winSize = int(end / groupNum)
    # Create list of window start points, ensuring the last point matches 'end'
    configWindowList = list(range(1, end, winSize)) + [end]
    # Generate windows based on configWindowList
    for i in range(1, len(configWindowList)):
        restmp = [chr, configWindowList[i-1], configWindowList[i]]
        res.append(restmp)
    return res

# check re_exp
def check_string_re(re_exp, string):
    """
    Checks if a given string matches a specified regular expression pattern.

    Args:
        re_exp (str): The regular expression pattern to match.
        string (str): The string to be checked against the pattern.

    Returns:
        bool: True if the string matches the pattern, False otherwise.
    """
    return re.search(re_exp, string) is not None

# extract CB from pysam data
def extract_cbMeta_data(pysamData, outFilePath):
    cellBarcodesCheck = []
    file = open(outFilePath, 'w+')
    count_num = 0
    for read in pysamData.fetch():
        count_num += 1
        tagCheck = read.has_tag("CB")
        if count_num % 5000000 == 0:
            print("Processed " + str(count_num) + " reads")
        if tagCheck:
            cellbarcode_read = read.opt('CB').split('-')[0]
            if cellbarcode_read not in cellBarcodesCheck:
                cellBarcodesCheck.append(cellbarcode_read)
                file.writelines(cellbarcode_read + "\n")
        else:
            pass

# split bam to small bam based on CB
# "nM" tag (number of mismatches); default: None
# "NH" tag (number of hits); default: None
def split_bam_data_obo(pysamData, barcodeList, subam_path, seqTag, num_Mis = None, num_Hit = None):
    count_num = 0
    DICT_files = {}
    # creat and open small pysam files
    print("###########################################")
    print("Split " + seqTag)
    print("Creating tmp bam files in "+ subam_path)
    print("###########################################")
    for barcode in barcodeList:
        pathTmp = subam_path + barcode + ".tmp."+ seqTag +".bam"
        Outfile = pysam.AlignmentFile(pathTmp, "wb", template=pysamData)
        DICT_files[barcode] = Outfile
    # split bam based on barcode
    for read in pysamData.fetch():
        count_num += 1
        tagCheck = read.has_tag("CB")
        if count_num % 500000 == 0:
            print("Processed " + str(count_num) + " reads")
        if tagCheck:
            cb_read = read.opt('CB').split('-')[0]
            if cb_read in barcodeList:
                DICT_files[cb_read].write(read)
        else:
            pass
    # pysamData.close()
    # close bam files
    print("Close pysam files...")
    for key in DICT_files.keys():
        DICT_files[key].close()
        pysam.index(DICT_files[key])
    print("Split bam ###---### Done!")

# shuffle cell barcodes
def shuffle_cell_barcodes(barcodeList, groupNum):
    n_total = len(barcodeList)
    shuffle(barcodeList)
    m = int(n_total/groupNum)
    res = []
    for i in range(0, n_total, m):
        res.append(barcodeList[i:i+m])
    print("Shuffle barcodes ###---### Done!")
    return res

# count mutations
"""
#Information on match, mismatch, indel, strand, mapping
#quality and start and end of a read are all encoded at the
#read base column. At this column, a dot stands for a match
#to the reference base on the forward strand, a comma for a
#match on the reverse strand, a '>' or '<' for a reference
#skip, `ACGTN' for a mismatch on the forward strand and
#`acgtn' for a mismatch on the reverse strand. A pattern
#`\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion
#between this reference position and the next reference
#position. The length of the insertion is given by the
#integer in the pattern, followed by the inserted
#sequence. Similarly, a pattern `-[0-9]+[ACGTNacgtn]+'
#represents a deletion from the reference. The deleted bases
#will be presented as `*' in the following lines. Also at
#the read base column, a symbol `^' marks the start of a
#read. The ASCII of the character following `^' minus 33
#gives the mapping quality. A symbol `$' marks the end of a
#read segment
"""
# count alt
def count_alt_num(pileupList, refBase):
    #count_num = 0
    CountList_Total = {"Mis": 0, "Insert": 0, "Del": 0, "Skip": 0, "DR": 0, "Other": 0}
    CountList_Forward = {"Mis": 0, "Insert": 0, "Del": 0, "Skip": 0}
    CountList_REVERSE = {"Mis": 0, "Insert": 0, "Del": 0, "Skip": 0}
    BASIC_BASE = ["A", "T", "C", "G", "N"]
    for pileupItem in pileupList:
        len_pileupItem = len(pileupItem)
        if len_pileupItem == 1:
            if pileupItem.upper() != refBase.upper():
                # check Forward Mis
                if pileupItem.isupper() and pileupItem.upper() in BASIC_BASE:
                    CountList_Forward["Mis"] += 1
                # check REVERSE Mis
                elif pileupItem.islower() and pileupItem.upper() in BASIC_BASE:
                    CountList_REVERSE['Mis'] += 1
                # check Forward Skip
                elif pileupItem == ">":
                    CountList_Forward["Skip"] += 1
                # check REVERSE Skip
                elif pileupItem == "<":
                    CountList_REVERSE["Skip"] += 1
                elif pileupItem == "*":
                    CountList_Total["DR"] += 1
            else:
                pass
        elif len_pileupItem > 1:
            if check_string_re('\+[0-9]+[ACGTN]+', pileupItem):
                CountList_Forward["Insert"] += 1
            elif check_string_re('\+[0-9]+[acgtn]+', pileupItem):
                CountList_REVERSE["Insert"] += 1
            elif check_string_re('\-[0-9]+[ACGTN]+', pileupItem):
                CountList_Forward["Del"] += 1
            elif check_string_re('\-[0-9]+[acgtn]+', pileupItem):
                CountList_REVERSE["Del"] += 1
            else:
                CountList_Total["Other"] += 1
        else:
            CountList_Total["Other"] += 1     
    CountList_Total["Mis"] = CountList_Forward["Mis"] + CountList_REVERSE["Mis"]
    CountList_Total["Insert"] = CountList_Forward["Insert"] + CountList_REVERSE["Insert"]
    CountList_Total["Del"] = CountList_Forward["Del"] + CountList_REVERSE["Del"]
    CountList_Total["Skip"] = CountList_Forward["Skip"] + CountList_REVERSE["Skip"]
    ALT_NUM = CountList_Total["Mis"] + CountList_Total["Insert"] + CountList_Total["Del"]
    return CountList_Total, CountList_Forward, CountList_REVERSE, ALT_NUM

# count base number exclude other base
def count_base_fromDic(baseDic, excludeBaseList):
    """
    Counts the number of bases in a dictionary, excluding specified bases.
    Args:
        baseDic (dict): Dictionary containing base counts.
        excludeBaseList (list): List of bases to exclude from the count.
    Returns:
        int: The count of bases excluding the specified bases.
    """
    check_bases = list(set(baseDic.keys()).difference(set(excludeBaseList)))
    ALT_NUM = sum(baseDic[base] for base in check_bases)
    return ALT_NUM

# count alt base types
# K = ['C', 'T', 'T', 'T', 'a', '<', '<', 't', '>', 'n', '+2AG', '-2ct', "*"]
def count_base_typeNum(pileupList, refBase, SkipCheck):
    """
    Counts the occurrences of different base types in a given pileup list and categorizes them.
    Args:
        pileupList (list): A list of bases and other characters representing the pileup.
        refBase (str): The reference base to compare against.
        SkipCheck (bool): A flag to determine whether to skip certain checks.
    Returns:
        tuple: A tuple containing:
            - NEW_PILE_LIST (list): A list of categorized bases and characters.
            - ALT_NUM (int): The count of alternative bases excluding the reference base and certain other categories.
    """
    BASIC_BASE = ["A", "T", "C", "G", "N"]
    NEW_PILE_LIST = []
    NEW_BASE_DIC = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, "D": 0, "I": 0, "FS": 0, "RS": 0, "O": 0, "NA": 0}
    for pileupItem in pileupList:
        if pileupItem.upper() in BASIC_BASE:
            NEW_BASE_DIC[pileupItem.upper()] += 1
            NEW_PILE_LIST.append(pileupItem.upper())
        else:
            if SkipCheck:
                if pileupItem.upper() == ">":
                    NEW_BASE_DIC["FS"] += 1
                    NEW_PILE_LIST.append("FS")
                elif pileupItem.upper() == "<":
                    NEW_BASE_DIC["RS"] += 1
                    NEW_PILE_LIST.append("RS")
                elif check_string_re(r'\+[0-9]+[ACGTNacgtn]+', pileupItem):
                    NEW_BASE_DIC["I"] += 1
                    NEW_PILE_LIST.append("I")
                elif check_string_re(r'\-[0-9]+[ACGTNacgtn]+', pileupItem):
                    NEW_BASE_DIC["D"] += 1
                    NEW_PILE_LIST.append("D")
                elif pileupItem.upper() == "*":
                    NEW_BASE_DIC["O"] += 1
                    NEW_PILE_LIST.append("O")
                else:
                    NEW_BASE_DIC["NA"] += 1
                    NEW_PILE_LIST.append("NA")
            else:
                if check_string_re(r'\+[0-9]+[ACGTNacgtn]+', pileupItem):
                    NEW_BASE_DIC["I"] += 1
                    NEW_PILE_LIST.append("I")
                elif check_string_re(r'\-[0-9]+[ACGTNacgtn]+', pileupItem):
                    NEW_BASE_DIC["D"] += 1
                    NEW_PILE_LIST.append("D")
                elif pileupItem.upper() == "*":
                    NEW_BASE_DIC["O"] += 1
                    NEW_PILE_LIST.append("O")
                else:
                    NEW_BASE_DIC["NA"] += 1
                    NEW_PILE_LIST.append("NA")
    ALT_NUM = count_base_fromDic(NEW_BASE_DIC, [refBase, "FS", "RS", "O", "NA", "N"])
    return NEW_PILE_LIST, ALT_NUM

# return cell mutations info
def get_mut_bycellExcludeRef(CELL_COUNTS_DIC, refBase):
    basic_base = ["A", "T", "C", "G", "D", "I"]
    check_base = list(set(basic_base).difference(set([refBase])))
    totalCells = list(CELL_COUNTS_DIC[check_base[0]])
    check_base_len = int(len(check_base))
    for i in range(1, check_base_len):
        totalCells = totalCells + list(CELL_COUNTS_DIC[check_base[i]])
    totalCells = set(totalCells)
    cellOut = []
    countOut = []
    for cell in totalCells:
        num_a = CELL_COUNTS_DIC['A'].count(cell)
        num_t = CELL_COUNTS_DIC['T'].count(cell)
        num_c = CELL_COUNTS_DIC['C'].count(cell)
        num_g = CELL_COUNTS_DIC['G'].count(cell)
        num_d = CELL_COUNTS_DIC['D'].count(cell)
        num_i = CELL_COUNTS_DIC['I'].count(cell)
        count_info = "|".join([str(num_a), str(num_t), str(num_c), str(num_g), str(num_d), str(num_i)])
        cellOut.append(cell)
        countOut.append(count_info)
    return cellOut, countOut

# return cell mutations info
def get_mut_bycell(CELL_COUNTS_DIC):
    """
    Returns the mutation information for each cell based on the provided dictionary of cell counts.

    Args:
        CELL_COUNTS_DIC (dict): Dictionary containing cell counts for each base type.

    Returns:
        tuple: A tuple containing:
            - cellOut (list): List of unique cell identifiers.
            - countOut (list): List of mutation counts for each cell in the format "A|T|C|G|D|I".
    """
    totalCells = set(CELL_COUNTS_DIC["A"] + CELL_COUNTS_DIC["T"] + CELL_COUNTS_DIC["C"] + CELL_COUNTS_DIC["G"] + CELL_COUNTS_DIC["D"] + CELL_COUNTS_DIC["I"])
    cellOut = []
    countOut = []
    ca = collections.Counter(CELL_COUNTS_DIC["A"])
    ct = collections.Counter(CELL_COUNTS_DIC["T"])
    cc = collections.Counter(CELL_COUNTS_DIC["C"])
    cg = collections.Counter(CELL_COUNTS_DIC["G"])
    cd = collections.Counter(CELL_COUNTS_DIC["D"])
    ci = collections.Counter(CELL_COUNTS_DIC["I"])
    for cell in totalCells:
        num_a = ca.get(cell, 0)
        num_t = ct.get(cell, 0)
        num_c = cc.get(cell, 0)
        num_g = cg.get(cell, 0)
        num_d = cd.get(cell, 0)
        num_i = ci.get(cell, 0)
        count_info = "|".join([str(num_a), str(num_t), str(num_c), str(num_g), str(num_d), str(num_i)])
        cellOut.append(cell)
        countOut.append(count_info)
    return cellOut, countOut

def get_mut_total(CELL_COUNTS_DIC, refBase):
    """
    Calculate the total number of mutations and cell counts for each base type.
    Args:
        CELL_COUNTS_DIC (dict): Dictionary containing cell counts for each base type.
        refBase (str): Reference base to exclude from the count.
    Returns:
        tuple: A tuple containing:
            - countCellNumSTR (str): String representation of cell counts for each base type.
            - countAltNumSTR (str): String representation of alternative base counts for each base type.
            - countCellNum (int): Total number of cells.
            - countAltNum (int): Total number of alternative bases.
            - countAltCellNum (int): Total number of cells with alternative bases.
    """
    basic_base = ["A", "T", "C", "G", "D", "I"]
    CELL_COUNTS_DIC_NEW_0 = {x: len(set(CELL_COUNTS_DIC[x])) for x in basic_base}
    CELL_COUNTS_DIC_NEW_1 = {x: len(CELL_COUNTS_DIC[x]) for x in basic_base}
    countCellNumSTR = ':'.join([str(CELL_COUNTS_DIC_NEW_0[x]) for x in basic_base])
    countAltNumSTR = ':'.join([str(CELL_COUNTS_DIC_NEW_1[x]) for x in basic_base])
    countCellNum = sum(CELL_COUNTS_DIC_NEW_0.values())
    countAltNum = count_base_fromDic(CELL_COUNTS_DIC_NEW_1, [refBase])
    countAltCellNum = count_base_fromDic(CELL_COUNTS_DIC_NEW_0, [refBase])
    return countCellNumSTR, countAltNumSTR, countCellNum, countAltNum, countAltCellNum

def integrate_data_fromQualityDic(baseQualityDic):
    """
    Integrates base quality data from a dictionary into a formatted string.
    Args:
        baseQualityDic (dict): Dictionary containing base quality data for each base type.
    Returns:
        str: Formatted string representing the integrated base quality data.
    """
    infoA = ":".join(str(i) for i in baseQualityDic["A"])
    infoT = ":".join(str(i) for i in baseQualityDic["T"])
    infoC = ":".join(str(i) for i in baseQualityDic["C"])
    infoG = ":".join(str(i) for i in baseQualityDic["G"])
    infoD = ":".join(str(i) for i in baseQualityDic["D"])
    infoI = ":".join(str(i) for i in baseQualityDic["I"])
    res = "|".join([infoA, infoT, infoC, infoG, infoD, infoI])
    return res

def cal_pError_fromBaseQuality(baseQuality):
    """
    Calculate the probability of error from base quality score.
    Args:
        baseQuality (int): The base quality score.
    Returns:
        float: The probability of error.
    """
    pError = 10 ** (-int(baseQuality) / 10)
    return pError

def count_mut_information(BamPath, referencePath, outVcfPath, config_posList, min_bq, min_mq, dp_check, alt_check, cell_check, nm_check, min_cov, min_alt, min_cell, max_nm, cbtag):
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
            DP_ALT_REAL = count_base_fromDic(BASE_COUNTS_DIC, [REF_BASE, "N"])
            if alt_check == True and int(DP_ALT_REAL) < int(min_alt):
                continue
            if DP_TOTAL >= 1000 and DP_ALT_REAL <= 0.005 * DP_TOTAL:
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

# target extract snp from bam
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
    pysamData = pysam.AlignmentFile(BamPath, 'rb')
    # chr
    chr = str(config_posList[0])
    # start
    start = int(config_posList[1]) - 100
    # end
    end = int(config_posList[2]) + 100
    # ref
    refFasta = pysam.FastaFile(referencePath)
    # pileups
    bam_pileups = pysamData.pileup(contig = str(chr), start = int(start), end = int(end), max_depth = 100000000, min_base_quality = min_bq, min_mapping_quality = min_mq)
    count_num_iter = 0
    # out vcf
    outPath = os.path.join(outVcfPath, f"VCF_{chr}_{config_posList[1]}_{config_posList[1]}_random.vcf")
    outVcfFile = open(outPath, "w+")
    for bam_pileup in bam_pileups:
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

# generate matrix data
def generate_matrix_fromvcf(vcfilePath, barcodePath, posIndexPath):
    # barcode index
    barcodeList = []
    barcodeListFile = open(barcodePath, "r")
    for line in barcodeListFile:
        if line.strip() not in barcodeList:
            barcodeList.append(line.strip())
    # position index
    posList = []
    posFile = open(posIndexPath, "r")
    for line in posFile:
        if line.strip() not in posList:
            posList.append(line.strip())
    # vcf loading
    tmpVcf = open(vcfilePath, "r")
    for line in tmpVcf:
        tagSets = line.split("\t")
        #print(line.strip())
        posTag = tagSets[0] + "_" + tagSets[1] + "_" + tagSets[2]
        cellList = tagSets[9].split(";")
        cellAltNumList = tagSets[10].split(";")
        cellIndex = [barcodeList.index(i) for i in cellList]
        posIndex = posList.index(posTag)
        #sql = str(posIndex) + "\t" + str(cellIndex) + "\t" + str(cellAltNumList)
        cellIndex_len = int(len(cellIndex))
        for i in range(0, cellIndex_len):
            sql = str(posIndex) + "\t" + str(cellIndex[i]) + "\t" + str(cellAltNumList[i])
            print(posTag + "\t" + sql)


############################################################
############################################################
# calculate probability of genotype
# 2024-03-11
############################################################
############################################################
# basic function
# Q test r22 for n >= 14
# https://pubs.acs.org/doi/pdf/10.1021/ac00002a010

Q90 = {6: 0.983, 7: 0.881, 9: 0.737, 10: 0.682, 11: 0.637, 12: 0.600, 13: 0.570,
       14: 0.546, 15: 0.525, 16: 0.507, 17: 0.490, 18: 0.475, 19: 0.462, 20: 0.450,
       21: 0.440, 22: 0.430}

Q95 = {6: 0.990, 7: 0.909, 8: 0.846, 9: 0.787, 10: 0.734, 11: 0.688, 12: 0.648, 13: 0.616,
       14: 0.590, 15: 0.568, 16: 0.548, 17: 0.531, 18: 0.516, 19: 0.503, 20: 0.491,
       21: 0.480, 22: 0.470}

Q99 = {6: 0.998, 7: 0.970, 8: 0.922, 9: 0.873, 10: 0.826, 11: 0.781, 12: 0.740, 13: 0.705,
       14: 0.674, 15: 0.647, 16: 0.624, 17: 0.605, 18: 0.589, 19: 0.575, 20: 0.562,
       21: 0.551, 22: 0.541}

def dixon_test_pGenotype(pGenotypeDic, left=True, right=True, q_dict=Q95):
    """
    Perform Dixon's Q test to identify outliers in genotype probabilities.

    Args:
        pGenotypeDic (dict): Dictionary of genotype probabilities.
        left (bool): Perform Q-test on the minimum value if True.
        right (bool): Perform Q-test on the maximum value if True.
        q_dict (dict): Dictionary of critical Q values for given confidence levels.

    Returns:
        list: List of two values for the outliers, or None.
    """
    data = list(pGenotypeDic.values())
    data_len = len(data)
    sdata = sorted(data)
    Q_mindiff, Q_maxdiff = (0, 0), (0, 0)
    if left:
        Q_min = (sdata[1] - sdata[0])
        try:
            Q_min /= (sdata[-1] - sdata[0])
        except ZeroDivisionError:
            pass
        Q_mindiff = (Q_min - q_dict.get(data_len, 0), sdata[0])   
    if right:
        Q_max = abs((sdata[-2] - sdata[-1]))
        try:
            Q_max /= abs((sdata[0] - sdata[-1]))
        except ZeroDivisionError:
            pass
        Q_maxdiff = (Q_max - q_dict.get(data_len, 0), sdata[-1])
    if not Q_mindiff[0] > 0 and not Q_maxdiff[0] > 0:
        outliers = [None, None]
        # Try z-score
        outliers_zscore = []
        mean = np.mean(sdata)
        std = np.std(sdata)
        # 99.7% threshold
        threshold = 3
        for x in sdata:
            z_score = (x - mean) / std
            if abs(z_score) > threshold:
                outliers_zscore.append(x)
        if outliers_zscore:
            outliers = [None, max(outliers_zscore)]
    elif Q_mindiff[0] == Q_maxdiff[0]:
        outliers = [Q_mindiff[1], Q_maxdiff[1]]
    elif Q_mindiff[0] > Q_maxdiff[0]:
        outliers = [Q_mindiff[1], None]
    else:
        outliers = [None, Q_maxdiff[1]]
    gt0 = [k for k, v in pGenotypeDic.items() if v == outliers[0]] if outliers[0] is not None else None
    gt1 = [k for k, v in pGenotypeDic.items() if v == outliers[1]] if outliers[1] is not None else None
    return [gt0, gt1]

def define_accuracy_calling(baseQuality):
    """
    Calculate the accuracy of base calling from base quality score.
    Args:
        baseQuality (int): The base quality score.
    Returns:
        float: The accuracy of base calling.
    """
    # Default 10x single cell mapping quality is 255. Also, should keep accuracy != 0
    # in our analysis, we filtered reads with low quality of mapping
    callingErrorP = cal_pError_fromBaseQuality(baseQuality)
    accuracySequencingP = 1 - callingErrorP
    if accuracySequencingP == 1:
        accuracySequencingP = 0.999999
    return accuracySequencingP

def shuffle_reads_features(refBase, baseQualityInfoOnesite):
    """
    Extracts features from base quality information for shuffling reads.
    Args:
        refBase (str): The reference base.
        baseQualityInfoOnesite (str): Base quality information for one site.
    Returns:
        tuple: A tuple containing:
            - readcountsNumNormal (list): Normalized read counts.
            - baseWeight (list): Weights for each base.
            - tmpDic (dict): Dictionary of base qualities.
            - tmpDicP (dict): Dictionary of base quality probabilities.
            - readcountsWeight (list): Weights for read counts.
    """
    basicBase = ["A", "T", "C", "G", "D", "I"]
    cellNum = []
    tmpDic = {"A": [], "T": [], "C": [], "G": [], "D": [], "I": []}
    tmpDicP = {"A": [], "T": [], "C": [], "G": [], "D": [], "I": []}
    num = len(baseQualityInfoOnesite.split(";"))
    for i in range(num):
        baseQuaOneCell = baseQualityInfoOnesite.split(";")[i].split("|")
        readcount = 0
        for index_base, j in enumerate(baseQuaOneCell):
            if j:
                readcount += len(j.split(":"))
                tmpDic[basicBase[index_base]].extend(j.split(":"))
                tmpDicP[basicBase[index_base]].extend(int(m) for m in j.split(":"))
        cellNum.append(readcount)
    baseWeight = [len(value) for value in tmpDic.values()]
    baseNumTotal = sum(baseWeight)
    # Check baseWeight, when there is no refBase, should consider the refBase exists
    if baseWeight[basicBase.index(refBase)] == 0:
        refBaseWeight = cal_pError_fromBaseQuality(min(tmpDic["A"] + tmpDic["T"] + tmpDic["C"] + tmpDic["G"] + tmpDic["D"] + tmpDic["I"]))
        baseWeight[basicBase.index(refBase)] = refBaseWeight * baseNumTotal
        tmpDicP[refBase] = [1]
    readcountsNum = list(set(cellNum))
    readcountsWeight = [cellNum.count(i) for i in readcountsNum]
    readcountsNumNormal = readcountsNum
    if max(readcountsNum) >= 10:
        readcountsNumNormal = [math.ceil((10 * i) / max(readcountsNum)) for i in readcountsNum]
    return readcountsNumNormal, baseWeight, tmpDic, tmpDicP, readcountsWeight

def shuffle_reads_baseonFeatures(readcountsNum, baseWeight, tmpDic, tmpDicP, readcountsWeight):
    """
    Shuffle reads based on features extracted from base quality information.
    Args:
        readcountsNum (list): List of normalized read counts.
        baseWeight (list): Weights for each base.
        tmpDic (dict): Dictionary of base qualities.
        tmpDicP (dict): Dictionary of base quality probabilities.
        readcountsWeight (list): Weights for read counts.
    Returns:
        str: Integrated base quality information for the shuffled reads.
    """
    basicBase = ["A", "T", "C", "G", "D", "I"]
    readsNumSelect = random.choices(readcountsNum, weights=readcountsWeight, k=1)
    baseSelect = random.choices(basicBase, weights=baseWeight, k=readsNumSelect[0])
    baseQua = {base: [] for base in basicBase}  
    for base in basicBase:
        if tmpDic[base]:
            baseQua[base] = random.choices(tmpDic[base], weights=tmpDicP[base], k=baseSelect.count(base))
    integrate_baseQuality = "|".join([":".join(baseQua[base]) for base in basicBase])
    return integrate_baseQuality

def shuffle_reads_resample(refBase,baseQualityInfoOnesite):
    # Example input: 32||||20|;27:32||||20|;32:37|||||;22:41:32:32:37:27|||||;41|||||;|||27|20|
    basicBase = ["A", "T", "C", "G", "D", "I"]
    cellNum = []
    tmpDic = {"A": [], "T": [], "C":[], "G": [], "D": [], "I":[]}
    tmpDicP = {"A": [], "T": [], "C":[], "G": [], "D": [], "I":[]}
    num = len(baseQualityInfoOnesite.split(";"))
    for i in range(0, num):
        baseQuaOneCell = baseQualityInfoOnesite.split(";")[i].split("|")
        readcount = 0
        for index_base, j in enumerate(baseQuaOneCell):
            if j != "":
                readcount += len(j.split(":"))
                tmpDic[basicBase[index_base]] += [m for m in j.split(":")]
                tmpDicP[basicBase[index_base]] += [define_accuracy_calling(int(m)) for m in j.split(":")]
        cellNum.append(readcount)
    baseWeight = [len(value) for value in tmpDic.values()]
    baseNumTotal = sum(baseWeight)
    # check baseWeight, when there is no refbase, should consider the refbase exist
    if baseWeight[basicBase.index(refBase)] == 0:
        refBaseWeight = cal_pError_fromBaseQuality(min(tmpDic["A"] + tmpDic["T"] + tmpDic["C"] + tmpDic["G"] + tmpDic["D"] + tmpDic["I"]))
        baseWeight[basicBase.index(refBase)] = refBaseWeight * baseNumTotal
        tmpDicP[refBase] = [1]
    readcountsNum = list(set(cellNum))
    readcountsWeight = [cellNum.count(i) for i in readcountsNum]
    readsNumSelect = random.choices(readcountsNum, weights = readcountsWeight, k = 1)
    baseSelect = random.choices(basicBase, weights = baseWeight, k = readsNumSelect[0])
    baseQuaA = []
    baseQuaT = []
    baseQuaC = []
    baseQuaG = []
    baseQuaD = []
    baseQuaI = []
    if len(tmpDic["A"]) != 0:
        baseQuaA = random.choices(tmpDic["A"], weights = tmpDicP["A"], k = baseSelect.count("A"))
    if len(tmpDic["T"]) != 0:
        baseQuaT = random.choices(tmpDic["T"], weights = tmpDicP["T"], k = baseSelect.count("T"))
    if len(tmpDic["C"]) != 0:
        baseQuaC = random.choices(tmpDic["C"], weights = tmpDicP["C"], k = baseSelect.count("C"))
    if len(tmpDic["G"]) != 0:
        baseQuaG = random.choices(tmpDic["G"], weights = tmpDicP["G"], k = baseSelect.count("G"))
    if len(tmpDic["D"]) != 0:
        baseQuaD = random.choices(tmpDic["D"], weights = tmpDicP["D"], k = baseSelect.count("D"))
    if len(tmpDic["I"]) != 0:
        baseQuaI = random.choices(tmpDic["I"], weights = tmpDicP["I"], k = baseSelect.count("I"))
    integrate_baseQuality = "|".join([":".join(baseQuaA), ":".join(baseQuaT), ":".join(baseQuaC), ":".join(baseQuaG), ":".join(baseQuaD), ":".join(baseQuaI)])
    return integrate_baseQuality

def split_oneCellReadCount_num(refBase, oneCellReadCountBase):
    # 0|3|0|0|0|0
    # A|T|C|G|D|I
    basicBase = ["A", "T", "C", "G", "D", "I"]
    baseCountinfoList = oneCellReadCountBase.split("|")
    baseCountinfoList_len = int(len(baseCountinfoList))
    #base = []
    base = [basicBase[i] for i in range(0, baseCountinfoList_len) if baseCountinfoList[i] != 0]
    refNum = baseCountinfoList[basicBase.index(refBase)]
    altBase = list(set(base).difference(set([refBase])))
    altBaseNum = sum([int(baseCountinfoList[basicBase.index(i)]) for i in altBase])
    return refNum, altBaseNum

def split_baseQua_num(refBase, oneCellBaseQua):
    # Example: 32:32||||20|
    basicBase = ["A", "T", "C", "G", "D", "I"]
    baseQuainfoList = oneCellBaseQua.split("|")
    base = []
    base = [basicBase[i] for i in range(0, 6) if baseQuainfoList[i] != ""]
    #for i in range(0, 6):
    #    if baseQuainfoList[i] != "":
    #        base.append(basicBase[i])
    refNum = len(baseQuainfoList[basicBase.index(refBase)].split(":"))
    altBase = list(set(base).difference(set([refBase])))
    altBaseNum = sum([int(len(baseQuainfoList[basicBase.index(i)].split(":"))) for i in altBase])
    return refNum, altBaseNum

def count_CellNum_fromQuainfo(refBase, quainfoOneSite):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    altBase = list(set(basicBase).difference(set([refBase])))
    altIndex = [basicBase.index(i) for i in altBase]
    cellsQua = quainfoOneSite.split(";")
    altNum = 0
    for cellQua in cellsQua:
        baseQua = cellQua.split("|")
        baseQuaList = [baseQua[i] for i in altIndex]
        #print(baseQuaList)
        #tmp = []
        tmp = [1 for i in baseQuaList if i != ""]
        #for i in baseQuaList:
        #    if i != "":
        #        tmp.append(1)
        numCheck = sum(tmp)
        #print(numCheck)
        if numCheck != 0:
            altNum += 1
    refNum = len(cellsQua) - altNum
    return refNum, altNum

def split_mutnotmutCell_num(refBase, oneSiteBaseQua):
    # Example: 32||||20|;27:32||||20|;32:37|||||
    basicBase = ["A", "T", "C", "G", "D", "I"]
    altBase = list(set(basicBase).difference(set([refBase])))
    altBaseIndex = [basicBase.index(i) for i in altBase]
    cellBaseList = oneSiteBaseQua.split(";")
    totalCellNum = len(cellBaseList)
    tmpCount = 0
    for cellbase in cellBaseList:
        tmplist = 0
        baselist = cellbase.split("|")
        for base_index, basequa in enumerate(baselist):
            if basequa != "" and base_index in altBaseIndex:
                tmplist += 1
        if tmplist != 0:
            tmpCount += 1
    mutCellNum = tmpCount
    notmutCellNum = totalCellNum - mutCellNum
    return mutCellNum, notmutCellNum

def integrate_baseQuality_onesite(oneSiteBaseQua):
    # Example: 32||||20|;27:32||||20|;32:37|||||
    basicBase = ["A", "T", "C", "G", "D", "I"]
    baseQuaDic = {"A": [], "T": [], "C": [], "G": [], "D": [], "I": []}
    cellsBaseQuaList = oneSiteBaseQua.split(";")
    for cell in cellsBaseQuaList:
        cellBaseList = cell.split("|")
        for index_base, cellBase in enumerate(cellBaseList):
            if cellBase != "":
                indexBase = basicBase[index_base]
                baseQuaDic[indexBase].append(cellBase)
    integrate_baseQuality = "|".join([":".join(baseQuaDic["A"]), ":".join(baseQuaDic["T"]), ":".join(baseQuaDic["C"]), ":".join(baseQuaDic["G"]), ":".join(baseQuaDic["D"]), ":".join(baseQuaDic["I"])])
    return integrate_baseQuality

def count_cell_reads(baseinfoCell):
    t = baseinfoCell.split("|")
    m = [len(i.split(":")) for i in t]
    return(sum(m))

####################
# Model 1
####################
def calculate_probility_oneCellBinomial(refBase, baseQualityInfoOneCell):
    # g(0): aa
    # g(1): ra
    # g(2): rr
    basicBaseDic = {"A": [], "T": [], "C": [], "G": [], "D": [], "I": []}
    basicBase = ["A", "T", "C", "G", "D", "I"]
    altBase = list(set(basicBase).difference(set([refBase])))
    baseQuainfoList = baseQualityInfoOneCell.split("|")
    baseQuainfoList_len = int(len(baseQuainfoList))
    for i in range(0, baseQuainfoList_len):
        if baseQuainfoList[i] != "":
            basicBaseDic[basicBase[i]] += [cal_pError_fromBaseQuality(i) for i in list(baseQuainfoList[i].split(":"))]
    refReadsQua = basicBaseDic[refBase]
    refReadsNum = len(refReadsQua)
    altReadsQua = []
    for i in altBase:
        altReadsQua += basicBaseDic[i]
    altReadsNum = len(altReadsQua)
    totalReadsNum = refReadsNum + altReadsNum
    p_0 = 1/(2**totalReadsNum)
    # likelihood of g(0)
    p_tmp_0_0 = 1
    p_tmp_0_1 = 1
    # likelihood of g(2)
    p_tmp_2_0 = 1
    p_tmp_2_1 = 1
    for i in range(0, refReadsNum):
        p_tmp_0_0 = p_tmp_0_0 * 2 * refReadsQua[i]
        p_tmp_2_0 = p_tmp_2_0 * 2 * (1 - refReadsQua[i])
    for i in range(0, altReadsNum):
        p_tmp_0_1 = p_tmp_0_1 * 2 * (1 - altReadsQua[i])
        p_tmp_2_1 = p_tmp_2_1 * 2 * altReadsQua[i]
    l0 = p_0 * p_tmp_0_0 * p_tmp_0_1
    l1 = p_0
    l2 = p_0 * p_tmp_2_0 * p_tmp_2_1
    genotype = {"aa": l0, "ra": l1, "rr": l2}
    outGenotype = max(genotype, key = lambda x: genotype[x]) 
    return l0, l1, l2, outGenotype

####################
# Model 2
####################
def neg_log_likelihood(params, successes, trials):
    shape1, shape2 = params
    kk = -np.sum(betabinom.logpmf(successes, trials, shape1, shape2))
    return kk

def estimate_shape_parameters_mle(successes, trials):
    # successes = np.array([10, 20, 30, 40])
    # trials = np.array([20, 30, 40, 50])
    # Initial guess for parameters
    initial_guess = [1, 1]
    # Minimize negative log-likelihood to estimate parameters
    result = minimize(neg_log_likelihood, initial_guess, args=(successes, trials), method='L-BFGS-B')
    shape1_est, shape2_est = result.x
    return shape1_est, shape2_est

def get_Vcf_rownum(vcfPath):
    buffer = 1024 * 1024
    with open(vcfPath) as f:
        buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
        return sum(buf.count('\n') for buf in buf_gen)

def random_extractsnp_fromVCF(vcfPath, tmpOutPath, numSitesUsedForEst, tag = False):
    # vcfPath
    fileRowNum = get_Vcf_rownum(vcfPath)
    tmpOutFile = open(tmpOutPath, "w")
    if fileRowNum < numSitesUsedForEst:
        print("#####Sites number is lower than number setting######")
        numSitesUsedForEst = int(fileRowNum * 0.4)
    numRandomList = random.sample(range(1, fileRowNum), numSitesUsedForEst)
    for i in numRandomList:
        lineTmp = linecache.getline(vcfPath, i)
        lineSets = lineTmp.split("\t")
        #print(i, "\t", lineTmp)
        chr = lineSets[0]
        pos = lineSets[1]
        refBase = lineSets[2]
        # ignore per cell info
        if tag:
            Cells = lineSets[10].split(";")
            CellIndex = random.randint(0, int(len(Cells) - 1))
            baseCountSite = Cells[CellIndex]
            refBaseNumOneCell, altBaseNumOneCell = split_oneCellReadCount_num(refBase, baseCountSite)
        else:
            refBaseNumOneCell = 0
            altBaseNumOneCell = 0
        refCellNum, altCellNum = int(lineSets[6].strip()) - int(lineSets[8].strip()), lineSets[8].strip()
        refReadNum, altReadNum = int(lineSets[3].strip()) - int(lineSets[5].strip()), lineSets[5].strip()
        sql = str(chr) + "\t" + str(pos) + "\t" + str(refBase) + "\t" + str(refBaseNumOneCell) + "\t" + str(altBaseNumOneCell) + "\t" + str(refCellNum) + "\t" + str(altCellNum) + "\t" + str(refReadNum) + "\t" + str(altReadNum)
        tmpOutFile.writelines(sql + "\n")

def generate_success_trailArray(randomExtractSNPPath):
    successesSingleCellList = []
    trialsSingleCellList = []
    successesCellNumList = []
    trialsCellNumList = []
    successesSiteList = []
    trialsSiteList = []
    with open(randomExtractSNPPath, "r") as f:
        for line in f:
            lineSets = line.split("\t")
            successesSingleCellTmp = int(lineSets[4])
            trialsSingleCellTmp = int(lineSets[3]) + int(lineSets[4])
            successesCellNumTmp = int(lineSets[6])
            trialsCellNumTmp = int(lineSets[5]) + int(lineSets[6])
            if trialsCellNumTmp < 3:
                continue
            successesSiteTmp = int(lineSets[8])
            trialsSiteTmp = int(lineSets[7]) + int(lineSets[8])
            successesSingleCellList.append(successesSingleCellTmp)
            trialsSingleCellList.append(trialsSingleCellTmp)
            successesCellNumList.append(successesCellNumTmp)
            trialsCellNumList.append(trialsCellNumTmp)
            successesSiteList.append(successesSiteTmp)
            trialsSiteList.append(trialsSiteTmp)
    f.close()
    successesCellArray = np.array(successesSingleCellList)
    trialsCellArray = np.array(trialsSingleCellList)
    successesCellNArray = np.array(successesCellNumList)
    trialsCellNArray = np.array(trialsCellNumList)
    successesSiteArray = np.array(successesSiteList)
    trialsSiteArray = np.array(trialsSiteList)
    return successesCellArray, trialsCellArray, successesCellNArray, trialsCellNArray, successesSiteArray, trialsSiteArray

def bayes_posterior_genotypeOneCell(AltReadNum, TotalReadsNumber, alpha, beta):
    p =  betabinom.sf(AltReadNum, TotalReadsNumber, alpha, beta)
    return p

####################
# Model 3
####################
def define_error_samplePreparing(genotype, altBase, perror):
    """
    Define the error probability for a given genotype and alternative base.

    Args:
        genotype (str): The genotype (e.g., "AA", "AT", etc.).
        altBase (str): The alternative base (e.g., "A", "T", "C", etc.).
        perror (float): The base error probability.

    Returns:
        float: The error probability for the given genotype and alternative base.
    """
    basicGenotype = ["AA", "AT", "AC", "AG", "AD", "AI", "TT", "TC", "TG", "TD", "TI", "CC", "CG", "CD", "CI", "GG", "GD", "GI", "DD", "DI", "II"]
    if genotype not in basicGenotype:
        raise ValueError("Invalid genotype")
    errorDic = {}
    if genotype == "AA":
        errorDic = {"A": 1 - 5 * perror, "T": perror, "C": perror, "G": perror, "D": perror, "I": perror}
    elif genotype == "AT":
        errorDic = {"A": (1 - 4 * perror) / 2, "T": (1 - 4 * perror) / 2, "C": perror, "G": perror, "D": perror, "I": perror}
    elif genotype == "AC":
        errorDic = {"A": (1 - 4 * perror) / 2, "T": perror, "C": (1 - 4 * perror) / 2, "G": perror, "D": perror, "I": perror}
    elif genotype == "AG":
        errorDic = {"A": (1 - 4 * perror) / 2, "T": perror, "C": perror, "G": (1 - 4 * perror) / 2, "D": perror, "I": perror}
    elif genotype == "AD":
        errorDic = {"A": (1 - 4 * perror) / 2, "T": perror, "C": perror, "G": perror, "D": (1 - 4 * perror) / 2, "I": perror}
    elif genotype == "AI":
        errorDic = {"A": (1 - 4 * perror) / 2, "T": perror, "C": perror, "G": perror, "D": perror, "I": (1 - 4 * perror) / 2}
    elif genotype == "TT":
        errorDic = {"A": perror, "T": 1 - 5 * perror, "C": perror, "G": perror, "D": perror, "I": perror}
    elif genotype == "TC":
        errorDic = {"A": perror, "T": (1 - 4 * perror) / 2, "C": (1 - 4 * perror) / 2, "G": perror, "D": perror, "I": perror}
    elif genotype == "TG":
        errorDic = {"A": perror, "T": (1 - 4 * perror) / 2, "C": perror, "G": (1 - 4 * perror) / 2, "D": perror, "I": perror}
    elif genotype == "TD":
        errorDic = {"A": perror, "T": (1 - 4 * perror) / 2, "C": perror, "G": perror, "D": (1 - 4 * perror) / 2, "I": perror}
    elif genotype == "TI":
        errorDic = {"A": perror, "T": (1 - 4 * perror) / 2, "C": perror, "G": perror, "D": perror, "I": (1 - 4 * perror) / 2}
    elif genotype == "CC":
        errorDic = {"A": perror, "T": perror, "C": 1 - 5 * perror, "G": perror, "D": perror, "I": perror}
    elif genotype == "CG":
        errorDic = {"A": perror, "T": perror, "C": (1 - 4 * perror) / 2, "G": (1 - 4 * perror) / 2, "D": perror, "I": perror}
    elif genotype == "CD":
        errorDic = {"A": perror, "T": perror, "C": (1 - 4 * perror) / 2, "G": perror, "D": (1 - 4 * perror) / 2, "I": perror}
    elif genotype == "CI":
        errorDic = {"A": perror, "T": perror, "C": (1 - 4 * perror) / 2, "G": perror, "D": perror, "I": (1 - 4 * perror) / 2}
    elif genotype == "GG":
        errorDic = {"A": perror, "T": perror, "C": perror, "G": 1 - 5 * perror, "D": perror, "I": perror}
    elif genotype == "GD":
        errorDic = {"A": perror, "T": perror, "C": perror, "G": (1 - 4 * perror) / 2, "D": (1 - 4 * perror) / 2, "I": perror}
    elif genotype == "GI":
        errorDic = {"A": perror, "T": perror, "C": perror, "G": (1 - 4 * perror) / 2, "D": perror, "I": (1 - 4 * perror) / 2}
    elif genotype == "DD":
        errorDic = {"A": perror, "T": perror, "C": perror, "G": perror, "D": 1 - 5 * perror, "I": perror}
    elif genotype == "DI":
        errorDic = {"A": perror, "T": perror, "C": perror, "G": perror, "D": (1 - 4 * perror) / 2, "I": (1 - 4 * perror) / 2}
    elif genotype == "II":
        errorDic = {"A": perror, "T": perror, "C": perror, "G": perror, "D": perror, "I": 1 - 5 * perror}

    return errorDic[altBase]

def calculate_probility_oneBase(k, obeservedBase, obeservedBaseQuality, genotype, perror):
    """
    Calculate the probability of observing a base given the genotype and error probability.
    Args:
        k (str): The base being considered (e.g., "A", "T", "C", etc.).
        obeservedBase (str): The observed base.
        obeservedBaseQuality (int): The quality score of the observed base.
        genotype (str): The genotype (e.g., "AA", "AT", etc.).
        perror (float): The base error probability.
    Returns:
        float: The probability of observing the base given the genotype and error probability.
    """
    perror_oneBase = define_error_samplePreparing(genotype=genotype, altBase=k, perror=perror)
    if obeservedBase == k:
        w = define_accuracy_calling(obeservedBaseQuality)
    else:
        w = (1 - define_accuracy_calling(obeservedBaseQuality)) / 5
    probilityGT_onebase = perror_oneBase * w
    return probilityGT_onebase

def calculate_probility_oneRead(observedBase, observedBaseQuality, genotype, perror):
    """
    Calculate the probability of observing a base given the genotype and error probability.
    Args:
        observedBase (str): The observed base.
        observedBaseQuality (int): The quality score of the observed base.
        genotype (str): The genotype (e.g., "AA", "AT", etc.).
        perror (float): The base error probability.
    Returns:
        float: The probability of observing the base given the genotype and error probability.
    """
    basicBase = ["A", "T", "C", "G", "D", "I"]
    probabilityGT_oneread = sum(
        calculate_probility_oneBase(k=base, obeservedBase=observedBase, obeservedBaseQuality=observedBaseQuality, genotype=genotype, perror=perror)
        for base in basicBase
    )
    return probabilityGT_oneread

def create_model3_probabilityDic(errorp):
    """
    Create a dictionary of probabilities for model 3 based on error probability.
    Args:
        errorp (float): The base error probability.
    Returns:
        dict: A nested dictionary containing probabilities for each genotype, observed base, and quality score.
    """
    outDic = {}
    genotypes = ["AA", "AT", "AC", "AG", "AD", "AI", "TT", "TC", "TG", "TD", "TI", "CC", "CG", "CD", "CI", "GG", "GD", "GI", "DD", "DI", "II"]
    observed_bases = ["A", "T", "C", "G", "D", "I"]
    quainfos = list(range(20, 60, 1))
    # Initialize outDic with nested dictionaries
    for genotype in genotypes:
        outDic[genotype] = {}
        for observed_base in observed_bases:
            outDic[genotype][observed_base] = {}
    # Now, fill outDic with values
    for genotype in genotypes:
        for observed_base in observed_bases:
            for quainfo in quainfos:
                k = calculate_probility_oneRead(observed_base, quainfo, genotype, errorp)
                outDic[genotype][observed_base][str(quainfo)] = k
    return outDic

def calculate_probility_oneGenotype(baseQuaDic, genotype, model3ProbabilityDic, quainfos, perror):
    probility_oneGenotype = 1
    for key, value in baseQuaDic.items():
        if value != "":
            num = len(value.split(':'))
            if num == 1:
                if value in quainfos:
                    tmpProbability = model3ProbabilityDic[genotype][key][value]
                else:
                    tmpProbability = calculate_probility_oneRead(key, value, genotype, perror)
                probility_oneGenotype = probility_oneGenotype * tmpProbability
            else:
                for i in range(0, num):
                    if value in quainfos:
                        tmpProbability = model3ProbabilityDic[genotype][key][value.split(':')[i]]
                    else:
                        tmpProbability = calculate_probility_oneRead(key, value.split(':')[i], genotype, perror)
                    probility_oneGenotype = probility_oneGenotype * tmpProbability
    return probility_oneGenotype

def calculate_posterior_genotypeOneCell(baseQualityInfo, model3ProbabilityDic, quainfos, errorP):
    #readsInfo:   Infomation of readcounts of each single cell. (Format: A|T|C|G|D|I; Example: 2|0|0|0|0|0 )
    #baseQuality: (Example: 27:27|||||)
    #Generally, only ten genotypes were considered
    basicGenotypesP = {"AA": 0, "AT": 0, "AC": 0, "AG": 0, "AD": 0, "AI": 0, "TT": 0, "TC": 0, "TG": 0, "TD": 0, "TI": 0, 
                       "CC": 0, "CG": 0, "CD": 0, "CI": 0, "GG": 0, "GD": 0, "GI": 0, "DD": 0, "DI": 0, "II": 0}

    baseQuaInfoList = baseQualityInfo.split("|")

    baseQuaDic = {"A": baseQuaInfoList[0], "T": baseQuaInfoList[1], "C": baseQuaInfoList[2], 
                  "G": baseQuaInfoList[3], "D": baseQuaInfoList[4], "I": baseQuaInfoList[5]}
    
    for genotype in basicGenotypesP.keys():
        basicGenotypesP[genotype] = calculate_probility_oneGenotype(baseQuaDic = baseQuaDic, genotype=genotype, model3ProbabilityDic = model3ProbabilityDic, quainfos = quainfos, perror=errorP)
    maxGenotype = max(basicGenotypesP, key = lambda x: basicGenotypesP[x])
    return maxGenotype, basicGenotypesP

def calculate_posterior_genotypeOneSite(baseQualityInfo, model3ProbabilityDic, quainfos, errorP):
    # 19	288738	T
    #site basequality information: ||41:32:41|||;||37:41:41|||;||37|||;||41|||;||41:37:41|||;||22:22|||;||37:37:37:41:37|||;||41:41:37|||;||41:27|||;||37:41:32|||;||41:41:41:41|||
    basicBase = ["A", "T", "C", "G", "D", "I"]
    tmpDic = {"A": [], "T": [], "C":[], "G": [], "D": [], "I":[]}
    num = len(baseQualityInfo.split(";"))
    tmpMinBaseQua = []
    for i in range(0, num):
        baseQuaOneCell = baseQualityInfo.split(";")[i].split("|")
        for index_base, j in enumerate(baseQuaOneCell):
            if j != "":
                tmpDic[basicBase[index_base]].append(j)
                tmpMinBaseQua += [int(k) for k in j.split(":")]
    integrate_baseQuality = "|".join([":".join(tmpDic["A"]), ":".join(tmpDic["T"]), ":".join(tmpDic["C"]), ":".join(tmpDic["G"]), ":".join(tmpDic["D"]), ":".join(tmpDic["I"])])
    p_site = calculate_posterior_genotypeOneCell(integrate_baseQuality,model3ProbabilityDic, quainfos, errorP)
    return p_site, min(tmpMinBaseQua)

####################
# main calculation
####################
def good_of_fitTest(refBase, baseQualityInfoOnesite, alphaAltCountEst, betaRefCountEst, alphaMutCellCountEst, betaRefCellCountEst, resampleCheck, model3ProbabilityDic, quainfos, resampleNum, iternumK, errorP, backP, m1check, m2check, m3check):
    # Example: ||41:32:41|||;||37:41:41|||;||37|||;||41|||;||41:37:41|||
    # extract base information
    #baseQuaInfoCellsList = baseQualityInfoOnesite.split(";")
    #cellNum = len(baseQuaInfoCellsList)
    integrate_baseQua = integrate_baseQuality_onesite(baseQualityInfoOnesite)
    totalRef, totalAlt = split_baseQua_num(refBase,integrate_baseQua)
    notmutCellNum, mutCellNum = count_CellNum_fromQuainfo(refBase, baseQualityInfoOnesite)
    mean_cell_reads = (totalRef + totalAlt) / (notmutCellNum + mutCellNum)
    if not resampleCheck:
        site_Genotype_m1 = {"rr": 0, "ra": 0, "aa": 0}
        site_Genotype_m2 = 1
        cell_num_m2 = 1
        site_Genotype_m3 = [None, None]
        if m1check:
            # Model 1
            site_Genotype_tmp = calculate_probility_oneCellBinomial(refBase, integrate_baseQua)[3]
            site_Genotype_m1[site_Genotype_tmp] += 1
        if m2check:
            # Model 2
            site_Genotype_m2 = bayes_posterior_genotypeOneCell(totalAlt, int(totalAlt + totalRef), alphaAltCountEst, betaRefCountEst)
            cell_num_m2 = bayes_posterior_genotypeOneCell(mutCellNum, int(mutCellNum + notmutCellNum), alphaMutCellCountEst, betaRefCellCountEst)
        if m3check:
            # Model 3
            site_Genotype_m3 = dixon_test_pGenotype(calculate_posterior_genotypeOneCell(integrate_baseQua, model3ProbabilityDic, quainfos, errorP)[1])
            #print(calculate_posterior_genotypeOneCell(integrate_baseQua, model3ProbabilityDic, quainfos, errorP)[1])
        result = (site_Genotype_m1, site_Genotype_m2, cell_num_m2, site_Genotype_m3)
    elif resampleCheck and totalAlt >= errorP * (totalAlt + totalRef):
        samples_random_Genotype_m1 = {"rr": 0, "ra": 0, "aa": 0}
        samples_random_Genotype_m2 = {"total": 0, "mut": 0}
        samples_random_Genotype_m3 = {"AA": 0, "AT": 0, "AC": 0, "AG": 0, "AD": 0, "AI": 0, "TT": 0, "TC": 0, "TG": 0, "TD": 0, "TI": 0, 
                                  "CC": 0, "CG": 0, "CD": 0, "CI": 0, "GG": 0, "GD": 0, "GI": 0, "DD": 0, "DI": 0, "II": 0}
        site_Genotype_m2 = 1
        cell_num_m2 = 1
        if m2check:
            # Model 2
            site_Genotype_m2 = bayes_posterior_genotypeOneCell(totalAlt, int(totalAlt + totalRef), alphaAltCountEst, betaRefCountEst)
            cell_num_m2 = bayes_posterior_genotypeOneCell(mutCellNum, int(mutCellNum + notmutCellNum), alphaMutCellCountEst, betaRefCellCountEst)
        # generate random samples
        samples_random_baseQuaInfo = []
        readcountsNum, baseWeight, tmpDic, tmpDicP, readcountsWeight = shuffle_reads_features(refBase,baseQualityInfoOnesite)
        for i in range(0, resampleNum):
            sample_random_baseQuaInfo_tmp = shuffle_reads_baseonFeatures(readcountsNum, baseWeight, tmpDic, tmpDicP, readcountsWeight)
            samples_random_baseQuaInfo.append(sample_random_baseQuaInfo_tmp)
        # setting iternum
        if mean_cell_reads >= 100:
            if notmutCellNum + mutCellNum < 50:
                iternum = iternumK
                sampleCellNum = 1
            else:
                iternum = iternumK
                sampleCellNum = 1
        else:
            if notmutCellNum + mutCellNum < 50:
                iternum = iternumK
                sampleCellNum = 3
            else:
                iternum = iternumK
                sampleCellNum = 5
        #counti = 0
        for i in range(0, iternum):
            newsiteQuainfo = random.sample(samples_random_baseQuaInfo, sampleCellNum)
            new_integrate_baseQua = integrate_baseQuality_onesite(";".join(newsiteQuainfo))
            totalRefK, totalAltK= split_baseQua_num(refBase,new_integrate_baseQua)
            if totalAltK >= 2:
                countCheck = count_cell_reads(new_integrate_baseQua)
                if countCheck >= 5000:
                    continue
                #counti += 1
                if m1check:
                    # Model 1
                    sample_random_probility_tmp_m1_1 = calculate_probility_oneCellBinomial(refBase, new_integrate_baseQua)[3]
                    samples_random_Genotype_m1[sample_random_probility_tmp_m1_1] += 1
                if m2check:
                    RefTmp, AltTmp = split_baseQua_num(refBase, new_integrate_baseQua)
                    m2_p_tmp = bayes_posterior_genotypeOneCell(AltTmp, int(AltTmp + RefTmp), alphaAltCountEst, betaRefCountEst)
                    if AltTmp > 0 and m2_p_tmp < 0.01:
                        samples_random_Genotype_m2["mut"] += 1
                if m3check:
                    # Model 3
                    maxGTTmp = calculate_posterior_genotypeOneCell(new_integrate_baseQua, model3ProbabilityDic, quainfos, errorP)
                    check_m3 = dixon_test_pGenotype(maxGTTmp[1])
                    if check_m3[0] != None or check_m3[1] != None:
                        samples_random_Genotype_m3[maxGTTmp[0]] += 1
            else:
                if m1check:
                    samples_random_Genotype_m1["rr"] += 1
                if m3check:
                    samples_random_Genotype_m3[refBase*2] += 1
        samples_random_Genotype_m2['total'] = iternum
        statpval_m1 = (1, 1)
        statpval_m2 = (1, 1)
        statpval_m3 = (1, 1)
        if m1check:
            statpval_m1 = proportions_ztest([iternum - samples_random_Genotype_m1["rr"], int(backP * iternum)], [iternum,iternum], alternative="larger")
        if m2check:
            statpval_m2 = proportions_ztest([samples_random_Genotype_m2["mut"], int(backP * iternum)], [iternum, iternum], alternative="larger")
        if m3check:
            statpval_m3 = proportions_ztest([iternum - samples_random_Genotype_m3[refBase*2], int(backP * iternum)], [iternum,iternum], alternative="larger")
        if m3check:
            iternum_alt_m3 = iternum - samples_random_Genotype_m3[refBase*2]
        if not m3check:
            iternum_alt_m3 = 0
        result = (statpval_m1[1], site_Genotype_m2, cell_num_m2, statpval_m2[1], statpval_m3[1], samples_random_Genotype_m1["ra"] + samples_random_Genotype_m1["aa"], samples_random_Genotype_m2["mut"], iternum_alt_m3)
    elif resampleCheck and totalAlt < errorP * (totalAlt + totalRef):
        site_Genotype_m2 = 0.05
        cell_num_m2 = 0.05
        result = (1, site_Genotype_m2, cell_num_m2, 1, 1, 10, 10, 10)
    return result

def assign_site_genotype(refBase, baseQualityInfoOnesite, alphaAltCountEst, betaRefCountEst, alphaMutCellCountEst, betaRefCellCountEst, resampleCheck, model3ProbabilityDic, quainfos, resampleNum, iternumK, errorP, backP, m1check, m2check, m3check):
    model_check = {"m1": 0, "m2": 0, "m3": 0, "p1": 1, "p2": 1, "p3": 1, "c1": 0, "c2": 0, "c3": 0, "ic1": 0, "ic2": 0, "ic3": 0}
    if resampleCheck:
        result = good_of_fitTest(refBase, baseQualityInfoOnesite, alphaAltCountEst, betaRefCountEst, alphaMutCellCountEst, betaRefCellCountEst, resampleCheck, model3ProbabilityDic, quainfos, resampleNum, iternumK, errorP, backP, m1check, m2check, m3check)
        # Model 1
        if result[0] <= 0.01:
            model_check["m1"] += 1
        # Model 2
        if result[1] <= 0.01 and result[2] <= 0.01 and result[3] <= 0.01:
            model_check["m2"] += 1
        # Model 3
        if result[4] <= 0.01:
            model_check["m3"] += 1
        model_check["p1"] = result[0]
        model_check["p2"] = result[3]
        model_check["p3"] = result[4]
        if result[5] > backP * resampleNum:
            model_check["c1"] += 1
        if result[6] > backP * resampleNum:
            model_check["c2"] += 1
        if result[7] > backP * resampleNum:
            model_check["c3"] += 1
        model_check["ic1"] = result[5]
        model_check["ic2"] = result[6]
        model_check["ic3"] = result[7]
    if not resampleCheck:
        result = good_of_fitTest(refBase, baseQualityInfoOnesite, alphaAltCountEst, betaRefCountEst, alphaMutCellCountEst, betaRefCellCountEst, resampleCheck, model3ProbabilityDic, quainfos, resampleNum, iternumK, errorP, backP, m1check, m2check, m3check)
        if m1check:
            # Model 1
            if max(result[0], key = lambda x: result[0][x]) != "rr":
                model_check["m1"] += 1
                model_check["p1"] = 0.001
                model_check["c1"] = 1
        if m2check:
            # Model 2
            if result[1] <= 0.01 and result[2] <= 0.01:
                model_check["m2"] += 1
                model_check["p2"] = 0.001
                model_check["c2"] = 1
        if m3check:
            # Model 3
            if result[3][0] != None or result[3][1] != None:
                if result[3][0] == None and refBase *2 not in result[3][1]:
                    model_check["m3"] += 1
                    model_check["p3"] = 0.001
                    model_check["c3"] = 1
                if result[3][1] == None and refBase *2 not in result[3][0]:
                    model_check["m3"] += 1
                    model_check["p3"] = 0.001
                    model_check["c3"] = 1
    return model_check