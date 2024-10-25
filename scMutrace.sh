#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

# Help function
helpFunction(){
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m                OneStep Process                    \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "Exclude background cells and trim sequencing reads"
    echo ""
    echo "Usage:"
    echo "$0 -b BAMFILE[str] -f REF[str] -c CELLBARCODE[str] -t ncut[int] -s SAMPLEID[str] -p nprocs[int] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -b \t[bam] bam file path to exclude background cells"
    echo -e "\t -f \t[ref] reference file path"
    echo -e "\t -c \t[barcode] cell barcode to be included"
    echo -e "\t -t \t[ncut] number of reads bases will be cut (default: 10)"
    echo -e "\t -s \t[id] sample id"
    echo -e "\t -p \t[nprocs] number of processes (default: 1)"
    echo -e "\t -O \t[outdir] output directory"
    echo -e "\t -h \t[help] display this help message"
    echo -e "\t Author: QunLi (liqun95@163.com) "
    exit 1
}
