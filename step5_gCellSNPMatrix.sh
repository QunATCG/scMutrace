#!/bin/bash
set -euo pipefail
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
    echo "$0 -g configfile[str] -O outputdir[str] -f reference[str] -q minBQ[int] -Q minMQ[int] -C cellbarcodetag[str] -p nprocs[int]"
    echo ""
    echo "Options:"
    echo -e "\t -g \t[config file path ] tab separate file, with the format [ID\tVCFPATH\tBAMPATH\tTAG]"
    echo -e "\t -O \t[output dir] output directory"
    echo -e "\t -f \t[reference] reference genome"
    echo -e "\t -q \t[min base quality] min base quality, default 0"
    echo -e "\t -Q \t[min mapping quality] min mapping quality, default 255"
    echo -e "\t -C \t[cell barcode tag] cell barcode tag, default CB"
    echo -e "\t -p \t[nprocs] number of processors, default 1"
    echo -e "\t -h \t[help] display this help message"
    echo -e "\t Author: QunLi (liqun95@163.com)"
    exit 1
}

# check if no arguments are passed
if [ $# -eq 0 ]; then
    helpFunction
fi

# set default values
min_bq=0
min_mq=255
cbtag="CB"
nprocs=1

# get parameters
while getopts "g:O:f:q:Q:C:p:h" opt
do
   case "$opt" in
        g ) contigfile="$OPTARG" ;;
        O ) outputdir="$OPTARG" ;;
        f ) ref="$OPTARG" ;;
        q ) min_bq="$OPTARG" ;;
        Q ) min_mq="$OPTARG" ;;
        C ) cbtag="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$contigfile" ] || [ -z "$outputdir" ] || [ -z "$ref" ]
then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the parameters are empty"
    [ -z "$contigfile" ] && echo "VCFDIR: ${contigfile}"
    [ -z "$outputdir" ] && echo "OUTPUT: ${outputdir}"
    [ -z "$ref" ] && echo "REF: ${ref}"
    helpFunction
fi

# echo all parameters
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
echo -e "\033[32m         Step5 generate Cell-SNP Matrix            \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "CONTIG: ${contigfile}"
echo "OUTPUT: ${outputdir}"
echo "REF: ${ref}"
echo "MIN BQ: ${min_bq}"
echo "MIN MQ: ${min_mq}"
echo "CBTAG: ${cbtag}"
echo "NPROCS: ${nprocs}"
echo ""
echo ""

# check input files
while IFS=$'\t' read -r sample vcf bam type
do
  if [[ ! -f "$vcf" ]]; then
    echo "[ERROR] VCF file not found for ${sample}: $vcf, please check the config file."
    exit 1
  fi

  if [[ ! -f "$bam" ]]; then
    echo "[ERROR] BAM file not found for ${sample}: $bam, please check the config file."
    exit 1
  fi
done < "${contigfile}"

# generate Cell-SNP Matrix
outputdir="${outputdir%/}"
python3 ${scMutrace_src_dir}/src/3_runStep1_gCellSNPMatrix.py --outdir ${outputdir} --jobconfig ${contigfile} --reference ${ref} --min_bq ${min_bq} --min_mq ${min_mq} --cbtag ${cbtag} --nprocs ${nprocs}
