import os
import numpy as np
from scipy.stats import fisher_exact

BASIC_BASE = ["A", "T", "C", "G", "D", "I"]

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

# get the number of reads supporting the reference and alternate alleles
def getMutRefReadsNum(numInfo, refBase):
    # A|T|C|G|D|I
    altBase = list(set(BASIC_BASE).difference(set([refBase])))
    altBaseIndex = [BASIC_BASE.index(i) for i in altBase]
    numList = numInfo.split(";")
    mutReadsNum = 0
    refReadsNum = 0
    mutFreqList = []
    mutCellNum = 0
    for numInfotmp in numList:
        countReads = numInfotmp.split("|")
        mutNum = sum([int(countReads[i]) for i in altBaseIndex])
        totalNum = sum([int(countReads[i]) for i in range(6)])
        mutFreqList.append(mutNum / totalNum)
        if mutNum > 0:
            mutCellNum += 1
        for base_index, basequa in enumerate(countReads):
            if basequa != 0 and base_index in altBaseIndex:
                mutReadsNum += int(basequa)
            elif basequa != 0 and base_index == BASIC_BASE.index(refBase):
                refReadsNum += int(basequa)
    return mutReadsNum, refReadsNum, mutFreqList, mutCellNum, len(numList)

# strand fisher test
def strand_fisher_test(alt_f, ref_f, alt_r, ref_r):
    """
    Performs a Fisher's exact test to determine if the strand bias for a variant is significant.

    Args:
        alt_f: Count of alt on forward strand
        ref_f: Count of ref on forward strand
        alt_r: Count of alt on reverse strand
        ref_r: Count of ref on reverse strand

    Returns:
        float: The p-value for the Fisher's exact test.
    """
    contingency_table = np.array([[alt_f, ref_f], [alt_r, ref_r]])
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')
    return p_value

# shuffle data for mitochondrial data
def filter_mito_VCF(vcfPathIn, vcfPathOut):
    """
    Filters a VCF file to include only variants on the mitochondrial chromosome.
    """
    vcfIn = open(vcfPathIn, "r")
    vcfOut = open(vcfPathOut, "w+")
    # write header
    vcfOut.write("chr\tpos\tref\ttotal_reads\ttotal_alt_reads\tmutFreq_total\tmutFreq_forward\tmutFreq_reverse\ttotal_cells\tmutCellnum_total\tmutCellnum_forward\tmutCellnum_reverse\ttotalCellNum_total\ttotalCellNum_forward\ttotalCellNum_reverse\tmutfrqlist_total\tcellFreq_total\tcellFreq_forward\tcellFreq_reverse\tp_value_strandFT\n")
    for line in vcfIn:
        lineSets = line.split("\t")
        # chr
        chr = lineSets[0]
        # pos
        pos = lineSets[1]
        # ref
        ref = lineSets[2]
        # total reads
        total_reads = lineSets[3]
        # total alt reads
        total_alt_reads = lineSets[5]
        # total cells
        total_cells = lineSets[9]
        # total reads info
        total_reads_info = lineSets[10]
        # forward reads info
        forward_reads_info = lineSets[13]
        # reverse reads info
        reverse_reads_info = lineSets[15]
        alt_total_reads, ref_total_reads, mutfrqlist_total, mutCellnum_total, totalCellNum_total = getMutRefReadsNum(total_reads_info, ref)
        alt_forward_reads, ref_forward_reads, mutfrqlist_forward, mutCellnum_forward, totalCellNum_forward = getMutRefReadsNum(forward_reads_info, ref)
        alt_reverse_reads, ref_reverse_reads, mutfrqlist_reverse, mutCellnum_reverse, totalCellNum_reverse = getMutRefReadsNum(reverse_reads_info, ref)
        # strand fisher test
        p_value_strandFT = strand_fisher_test(alt_forward_reads, ref_forward_reads, alt_reverse_reads, ref_reverse_reads)
        # variants alt freq
        mutFreq_total = alt_total_reads / (alt_total_reads + ref_total_reads)
        mutFreq_forward = alt_forward_reads / (alt_forward_reads + ref_forward_reads)
        mutFreq_reverse = alt_reverse_reads / (alt_reverse_reads + ref_reverse_reads)
        # mutfrqlist_total
        mutfrqlist_total_str = ";".join([str(i) for i in mutfrqlist_total])
        # mutfrqlist_forward
        #mutfrqlist_forward_str = ";".join([str(i) for i in mutfrqlist_forward])
        # mutfrqlist_reverse
        #mutfrqlist_reverse_str = ";".join([str(i) for i in mutfrqlist_reverse])
        # cell alt freq
        cellFreq_total = mutCellnum_total / totalCellNum_total
        cellFreq_forward = mutCellnum_forward / totalCellNum_forward
        cellFreq_reverse = mutCellnum_reverse / totalCellNum_reverse
        # write to file with header
        vcfOut.write(f"{chr}\t{pos}\t{ref}\t{total_reads}\t{total_alt_reads}\t{mutFreq_total}\t{mutFreq_forward}\t{mutFreq_reverse}\t{total_cells}\t{mutCellnum_total}\t{mutCellnum_forward}\t{mutCellnum_reverse}\t{totalCellNum_total}\t{totalCellNum_forward}\t{totalCellNum_reverse}\t{mutfrqlist_total_str}\t{cellFreq_total}\t{cellFreq_forward}\t{cellFreq_reverse}\t{p_value_strandFT}\n")

def main():
    vcfPathIn = "/Users/liqun/scMutraceProject/Data/GSE142742/data/VCF/VCF_CD34_800_Day08_MT.vcf"
    vcfPathOut = "/Users/liqun/scMutraceProject/Data/GSE142742/data/VCF/VCF_CD34_800_Day08_MT.vcf.out"
    filter_mito_VCF(vcfPathIn, vcfPathOut)

if __name__ == "__main__":
    main()