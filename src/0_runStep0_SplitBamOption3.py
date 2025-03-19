# import
import os
import pysam
import timeit
import argparse
from py_utils import make_dir_check
import multiprocessing as mp

'''
Option 1:
https://github.com/10XGenomics/subset-bam

Option 2:
There is also another method to extract bam with desired barcodes.
# filter.txt
CB:Z:GCCAAATTCACATACG-1
CB:Z:GGACATTGTGATGATA-1
CB:Z:TCAGGTACATTAGGCT-1

export BAM_FILE='possorted_genome_bam.bam'

# Save the header lines
samtools view -H $BAM_FILE > SAM_header

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $BAM_FILE | LC_ALL=C grep -F -f filter.txt > filtered_SAM_body

# Combine header and body
cat SAM_header filtered_SAM_body > filtered.sam

# Convert filtered.sam to BAM format
samtools view -b filtered.sam > filtered.bam

Option 3:
this script pysam
'''

# get argparse
def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to exclude noisy reads based on barcode')
    parser.add_argument('--bam', type=str, help='BAM file to be analysed', required = True)
    parser.add_argument('--ref', type=str, help='Path to reference genome version. *fai must be available in the same directory as the reference genome file', required = True)
    parser.add_argument('--barcode', type=str, help='Barcode path', required = True)
    parser.add_argument('--barcodetag', type=str, help='Barcode tag', required = True)
    parser.add_argument('--contigFile', type=str, help='Contig information in tab-separated format: Chr\tstart\tend', required = True)
    parser.add_argument('--splitCheck', action='store_true', help='A boolean flag')
    parser.add_argument('--id', help='Prefix used to name output file.', required = True)
    parser.add_argument('--ncut', default = 10, help='Number of cut reads bases [Default: 10]',required=False,type = int)
    parser.add_argument('--nprocs', default = 1, help='Number of processes [Default: 1]',required=False,type = int)
    parser.add_argument('--outdir', default = '.', help='Out folder', required = False)
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

# get config data from fasta
def generate_config_range(refasta, contigSelfInPath):
    #basicChr = list(str(i) for i in range(1, 40)) + ["chr" + str(i) for i in range(1, 40)] + ['X', 'Y', 'chrX', 'chrY', 'M', 'chrM', 'm', 'chrm', 'MT', 'chrMT', 'mt', 'chrmt','x', 'y', 'chrx', 'chry']
    chrListSelfDefined = []
    with open(contigSelfInPath, "r") as f:
        for line in f:
            chrListSelfDefined.append(line.split("\t")[0].strip())
    refChrTotal = refasta.references
    configDic = {}
    chrList = []
    for refChr in refChrTotal:
        if refChr in chrListSelfDefined:
            sqlConfig = f"{refChr}_1_{refasta.get_reference_length(refChr)}"
            chrList.append(refChr)
            configDic[refChr] = sqlConfig
    return configDic, chrList

# split bam based on barcode
def split_bam_basedList(bam_Path, barcodeSubList, ChrList, configValue, subam_path, seqTag, nameOut, ncut, barcodetag):
    bam = pysam.AlignmentFile(bam_Path, 'rb')
    count_num = 0
    chrtmp = str(configValue.split("_")[0])
    startmp = str(configValue.split("_")[1])
    endtmp = str(configValue.split("_")[2])
    # create and open small pysam files
    outBamPath = os.path.join(subam_path, f"{nameOut}_{chrtmp}_{startmp}_{endtmp}.tmp.{seqTag}.bam")
    Outfile = pysam.AlignmentFile(outBamPath, "wb", template=bam)
    print("###########################################")
    print("Split " + seqTag)
    print("Creating tmp bam files in " + outBamPath)
    print("###########################################")
    for read in bam.fetch(contig=chrtmp, start=int(startmp), stop=int(endtmp)):
        count_num += 1
        if count_num % 500000 == 0:
            print("Processed " + str(count_num) + " reads")
        if read.has_tag(barcodetag):
            cb_read = read.get_tag(barcodetag).split('-')[0]
            chr_read = read.reference_name
            if cb_read in barcodeSubList and chr_read in ChrList:
                # cut each read's first left and right bases
                CIGAR_read = read.cigartuples
                # check if CIGAR exists
                if CIGAR_read and len(CIGAR_read) > 1:
                    cut_l = CIGAR_read[0][1] + ncut if CIGAR_read[0][0] == 4 else ncut
                    cut_r = CIGAR_read[-1][1] + ncut if CIGAR_read[-1][0] == 4 else ncut
                else:
                    cut_l = ncut
                    cut_r = ncut
                base_quality = read.query_qualities
                index_cut_bq = list(range(0, cut_l)) + list(range(-1, -cut_r - 1, -1))
                if base_quality is not None:
                    for index_q in index_cut_bq:
                        base_quality[index_q] = 0
                    read.query_qualities = base_quality
                Outfile.write(read)
    Outfile.close()
    print("Indexing bam --- " + outBamPath)
    pysam.index(outBamPath)
    print("Split bam ###---### Done!")
    bam.close()

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    bam_Path = args.bam
    ref_Path = args.ref
    barcodePath = args.barcode
    barcodetag = args.barcodetag
    contigSelfInPath = args.contigFile
    split_check = args.splitCheck
    nameSplitBamTag = args.id
    seqTag = "SMT"
    ncut = args.ncut
    COREUSE = args.nprocs
    outDirPath = args.outdir
    # make dir
    make_dir_check(outDirPath)
    # extract barcode
    barcodeListCell = []
    barcodeListFile = open(outDirPath + "barcodeList.txt", "w+")
    totalBarcode = open(barcodePath, "r")
    for line in totalBarcode:
        barcodeTmp = line.split("\t")[0].strip().split("-")[0]
        if barcodeTmp not in barcodeListCell:
            barcodeListFile.writelines(barcodeTmp + "\n")
            barcodeListCell.append(barcodeTmp)
    barcodeListFile.close()
    # generate_config_range
    configDic, chrList = generate_config_range(pysam.FastaFile(ref_Path), contigSelfInPath)
    # extract bam by cellbarcode
    print(split_check)
    if split_check:
        if COREUSE > 1:
            print('Parent process %s.' % os.getpid())
            p = mp.get_context("fork").Pool(COREUSE)
            print("Your cores is " + str(COREUSE))
            results = []
            for value in configDic.values():
                p.apply_async(split_bam_basedList, args = (bam_Path, barcodeListCell, chrList, value, outDirPath, seqTag, nameSplitBamTag, ncut, barcodetag,), error_callback = err_call_back)
            print('Waiting for all subprocesses done...')
            p.close()
            p.join()
            for result in results:
                print(result.get())
            print('All subprocesses done.')
        else:
            for value in configDic.values():
                split_bam_basedList(bam_Path = bam_Path, barcodeSubList = barcodeListCell, ChrList = chrList, configValue = value, subam_path = outDirPath, seqTag = seqTag, nameOut = nameSplitBamTag, ncut=ncut, barcodetag=barcodetag)

if __name__=='__main__':
    start = timeit.default_timer()
    #19_1_59128983
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')