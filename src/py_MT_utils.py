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
    for numInfotmp in numList:
        countReads = numInfotmp.split("|")
        for base_index, basequa in enumerate(countReads):
            if basequa != 0 and base_index in altBaseIndex:
                mutReadsNum += int(basequa)
            elif basequa != 0 and base_index == BASIC_BASE.index(refBase):
                refReadsNum += int(basequa)
    return mutReadsNum, refReadsNum

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
        # forward reads info
        forward_reads_info = lineSets[13]
        # reverse reads info
        reverse_reads_info = lineSets[15]
        alt_forward_reads, ref_forward_reads = getMutRefReadsNum(forward_reads_info, ref)
        alt_reverse_reads, ref_reverse_reads = getMutRefReadsNum(reverse_reads_info, ref)
        # strand fisher test
        p_value_strandFT = strand_fisher_test(alt_forward_reads, ref_forward_reads, alt_reverse_reads, ref_reverse_reads)
        print(p_value_strandFT)

def main():
    vcfPathIn = "/Users/liqun/scMutraceProject/Data/GSE142742/data/VCF/VCF_CD34_800_Day08_MT.test.vcf"
    vcfPathOut = "/Users/liqun/scMutraceProject/Data/GSE142742/data/VCF/VCF_CD34_800_Day08_MT.test.vcf.out"
    filter_mito_VCF(vcfPathIn, vcfPathOut)

if __name__ == "__main__":
    main()