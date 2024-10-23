#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

# help info
helpFunction(){
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m            Step4 filter noisy sites               \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "Filter noisy sites"
    echo ""
    echo "Usage:"
    echo "$0 -v VCFDIR[str] -s SAMPLEID[str] -p nprocs[int] -r excludeitems[str] -i includeitems[str] -o OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -v \t[vcf path] vcf path"
    echo -e "\t -o \t[outdir] output directory"
    echo -e "\t -s \t[id] sample ID"
    echo -e "\t -c \t[barcode]  barcode file [tab separated file with header: Barcode CellType]"
    echo -e "\t -l \t[binsize] base number between two snps (MultiAlleles)"
    echo -e "\t -L \t[clustersize] number of snps in one cluster"
    echo -e "\t -d \t[minsequencingdepth] min depth"
    echo -e "\t -D \t[min altdpsize] min alternative bases depth"
    echo -e "\t -n \t[mincellnumber] min cell number"
    echo -e "\t -N \t[minaltcellnumber] min alteration cell number"
    echo -e "\t -B \t[lowBaseCutoff] cutoff of low base quality (default 20)"
    echo -e "\t -h \t[help] display this help message"
    echo -e "\t Author: QunLi (liqun95@163.com)"
    exit 1
}

