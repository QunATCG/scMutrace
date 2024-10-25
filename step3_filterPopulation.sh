#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

# help info
helpFunction(){
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m          Step3 filter public databases            \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "Filter SNPs based on public databases"
    echo ""
    echo "Usage:"
    echo "$0 -V VCFDIR[str] -s SAMPLEID[str] -p nprocs[int] -r excludeitems[str] -i includeitems[str] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -O \t[outdir] output directory"
    echo -e "\t -V \t[casevcfdir] case vcf directory"
    echo -e "\t -s \t[id] sample ID"
    echo -e "\t -r \t[removeitems] files to be excluded, tab seperated format: path\tid"
    echo -e "\t -i \t[includeitems] files to be included, tab seperated format: path\tid"
    echo -e "\t -p \t[nprocs] number of processes"
    echo -e "\t -h \t[help] display this help message"
    echo -e "\t Author: QunLi (liqun95@163.com)"
    exit 1
}

# check if no arguments are passed
if [ $# -eq 0 ]; then
    helpFunction
fi

# default parameters
nprocs=1
removeitem="NotAvailable"
includeitem="NotAvailable"

# get parameters
while getopts "O:V:s:r:i:p:h" opt
do
   case "$opt" in
        O ) outdir="$OPTARG" ;;
        V ) casevcfdir="$OPTARG" ;;
        s ) id="$OPTARG" ;;
        r ) removeitem="$OPTARG" ;;
        i ) includeitem="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$outdir" ] || [ -z "$casevcfdir" ] || [ -z "$id" ] 
then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the parameters are empty"
    [ -z "$outdir" ] && echo "Output directory is missing"
    [ -z "$casevcfdir" ] && echo "VCF directory is missing"
    [ -z "$id" ] && echo "SampleID is missing"
    helpFunction
fi

# echo all parameters
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
echo -e "\033[32m          Step3 filter public databases            \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "OUTPUT: ${outdir}"
echo "VCF: ${casevcfdir}"
echo "SAMPLEID: ${id}"
echo "ExcludeItems: ${removeitem}"
echo "IncludeItems: ${includeitem}"
echo "nprocs: ${nprocs}"
echo ""
echo ""

# main steps
python -u ${scMutrace_src_dir}/src/2_runStep1_FilterSNPHarder.py --outdir ${outdir} --casevcfdir ${casevcfdir} --caseid ${id} --exclude ${removeitem} --include ${includeitem} --nprocs ${nprocs}