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
    echo "$0 -v VCFPath[str] -s SAMPLEID[str] -c CELLBARCODEFILE[str] -l BINSIZE[int] -L CLUSTERSIZE[int] -d MINDP[int] -D MINALTDP[int] -n MINCELLNUM[int] -N MINALTCELLNUM[int] -B LOWBQ[int] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -v \t[vcf path] vcf path"
    echo -e "\t -O \t[outdir] output directory"
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

# check if no arguments are passed
if [ $# -eq 0 ]; then
    helpFunction
fi

# default parameters
binsize=10
clustersize=3
minsequencingdepth=5
minaltdpsize=2
mincellnumber=3
minaltcellnumber=2
lowBaseCutoff=20

# get parameters
while getopts "v:O:s:c:l:L:d:D:n:N:B:h" opt
do
   case "$opt" in
        v ) vcfpath="$OPTARG" ;;
        O ) outdir="$OPTARG" ;;
        s ) id="$OPTARG" ;;
        c ) cellbarcode="$OPTARG" ;;
        l ) binsize="$OPTARG" ;;
        L ) clustersize="$OPTARG" ;;
        d ) minsequencingdepth="$OPTARG" ;;
        D ) minaltdpsize="$OPTARG" ;;
        n ) mincellnumber="$OPTARG" ;;
        N ) minaltcellnumber="$OPTARG" ;;
        B ) lowBaseCutoff="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$vcfpath" ] || [ -z "$outdir" ] || [ -z "$cellbarcode" ] || [ -z "$id" ]
then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the parameters are empty"
    [ -z "$vcfpath" ] && echo "BAMFILE: ${vcfpath}"
    [ -z "$cellbarcode" ] && echo "REF: ${cellbarcode}"
    [ -z "$id" ] && echo "SAMPLEID: ${id}"
    [ -z "$outdir" ] && echo "OUTPUT: ${outdir}"
    helpFunction
fi

# echo all parameters
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
echo -e "\033[32m            Step4 filter noisy sites               \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "VCFPath: ${vcfpath}"
echo "SAMPLEID: ${id}"
echo "CELLBARCODEFILE: ${cellbarcode}"
echo "BINSIZE: ${binsize}"
echo "CLUSTERSIZE: ${clustersize}"
echo "MINDP: ${minsequencingdepth}"
echo "MINALTDP: ${minaltdpsize}"
echo "MINCELLNUM: ${mincellnumber}"
echo "MINALTCELLNUM: ${minaltcellnumber}"
echo "LOWBQ: ${lowBaseCutoff}"
echo "OUTPUT: ${outdir}"
echo ""
echo ""

# main steps
python -u ${scMutrace_src_dir}/src/3_runStep0_noisy.py --vcfpath ${vcfpath} --id ${id} --barcode ${cellbarcode} --binsize ${binsize} --clustersize ${clustersize} --dpsize ${minsequencingdepth} --altdpsize ${minaltdpsize} --cellsize ${mincellnumber} --altcellsize ${minaltcellnumber} --lowBaseCutoff ${lowBaseCutoff} --outdir ${outdir}