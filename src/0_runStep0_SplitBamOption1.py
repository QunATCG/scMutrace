import os
import timeit
import argparse

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
this script based on pysam
'''

# get argparse
def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to obtain a list of base and cell counts in scRNA bam file')
    parser.add_argument('--subbampath', type=str, default='', help='subam script path', required = True)
    parser.add_argument('--bampath', type=str, default='', help='BAM file to be analysed', required = True)
    parser.add_argument('--bamtag', type=str, default="CB", help='BAM tag used to be subset', required = False)
    parser.add_argument('--barcodepath', type=str, help='Barcode path', required = True)
    parser.add_argument('--outpath', default = '.', help='Out folder', required = True)
    parser.add_argument('--nprocs', default = 1, help='Number of processes [Default: 1]',required=False,type = int)
    return (parser)

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    subpath = args.subbampath
    inbam = args.bampath
    barcode = args.barcodepath
    outbam = args.outpath
    bamtag = args.bamtag
    CORES = args.nprocs 
    if os.path.exists(inbam + ".bai"):
        pass
    else:
        print("samtools index " + inbam)
        cmdIndex = "samtools index -@ " + str(CORES) + " " + inbam
        os.system(cmdIndex)
    print("extract bam based on barcodes")
    cmd = subpath + " --bam " + inbam + " --cell-barcodes " + barcode + " --cores " + str(CORES) + " --out-bam " + outbam + " --bam-tag " + bamtag 
    print(cmd)
    os.system(cmd)
    

if __name__=='__main__':
    start = timeit.default_timer()
    #19_1_59128983
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')

# python 0_runStep0_SplitBamOption1.py  --subbampath /Users/liqun/QunScripts/subset_bam --bampath /Users/liqun/BrainRNA/SNP/testData/scRNA/bam/P4_cSCC_scRNA.chr19.bam --barcodepath /Users/liqun/BrainRNA/SNP/testData/scRNA/meta/P4_Tumor_meta_sub --outpath /Users/liqun/BrainRNA/SNP/testData/scRNA/tmp/subam/P4_cSCC_exCludeBG.tmp_1.RNA.bam  --nprocs 10