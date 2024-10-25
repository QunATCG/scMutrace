#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

# help info
helpFunction(){
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m         Step5 generate Cell-SNP Matrix            \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "generate Cell-SNP Matrix"
    echo ""
    echo "Usage:"
    echo "$0 -v VCFPATH[str] -o OUTPUTPath[str]"
    echo ""
    echo "Options:"
    echo -e "\t -v \t[vcf path] vcf path"
    echo -e "\t -o \t[outputpath] output path"
    echo -e "\t -h \t[help] display this help message"
    echo -e "\t Author: QunLi (liqun95@163.com)"
    exit 1
}

# check if no arguments are passed
if [ $# -eq 0 ]; then
    helpFunction
fi

# get parameters
while getopts "v:o:h" opt
do
   case "$opt" in
        v ) vcfpath="$OPTARG" ;;
        o ) outputpath="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$vcfpath" ] || [ -z "$outputpath" ]
then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the parameters are empty"
    [ -z "$vcfpath" ] && echo "VCFDIR: ${vcfpath}"
    [ -z "$outputpath" ] && echo "OUTPUT: ${outputpath}"
    helpFunction
fi

# echo all parameters
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
echo -e "\033[32m         Step5 generate Cell-SNP Matrix            \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "VCFPATH: ${vcfpath}"
echo "OUTPUT: ${outputpath}"
echo ""
echo ""

# generate Cell-SNP Matrix
python3 ${scMutrace_src_dir}/3_runStep1_gCellSNPMatrix.py --vcfpath ${vcfpath} --outpath ${outputpath}
