#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

# help info
helpFunction(){
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m               Step1 extract SNPs                  \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "Extract SNPs from bam file"
    echo ""
    echo "Usage:"
    echo "$0 -b BAMFILE[str] -f REF[str] -g CONTIG[str] -s SAMPLEID[str] -p nprocs[int] -q minbq[int] -Q minmq[int] -n mincellnumber[int] -N minaltcellnumber[int] -m maxmismatch[int] -d minsequencingdepth[int] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -b \t[bam] bam file path to extract SNPs"
    echo -e "\t -f \t[ref] reference file path"
    echo -e "\t -C \t[cellbarcodetag] cell barcode tag in bam file (default: CB)"
    echo -e "\t -g \t[contig] contig information, tab separated format for each line: Chr\tstart\tend"
    echo -e "\t -s \t[id] sample id"
    echo -e "\t -p \t[nprocs] number of processes"
    echo -e "\t -R \t[ratio] ratio of cell number to total cell number"
    echo -e "\t -q \t[minbq] min base quality"
    echo -e "\t -Q \t[minmq] min mapping quality"
    echo -e "\t -n \t[mincellnumber] min cell number"
    echo -e "\t -N \t[minaltcellnumber] min alteration cell number"
    echo -e "\t -m \t[maxmismatch] number max mismatch number"
    echo -e "\t -d \t[minsequencingdepth] min sequecing depth"
    echo -e "\t -O \t[outdir] output directory"
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
minbq=20
minmq=255
mincellnumber=3
minaltcellnumber=2
maxmismatch=4
minsequencingdepth=3
cellbarcodetag="CB"
ratioUser=0.005

# get parameters
while getopts "b:f:C:g:s:p:R:q:Q:n:N:m:d:O:h" opt
do
   case "$opt" in
        b ) bam="$OPTARG" ;;
        f ) ref="$OPTARG" ;;
        C ) cellbarcodetag="$OPTARG" ;;
        g ) contig="$OPTARG" ;;
        s ) id="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        R ) ratioUser="$OPTARG" ;;
        q ) minbq="$OPTARG" ;;
        Q ) minmq="$OPTARG" ;;
        n ) mincellnumber="$OPTARG" ;;
        N ) minaltcellnumber="$OPTARG" ;;
        m ) maxmismatch="$OPTARG" ;;
        d ) minsequencingdepth="$OPTARG" ;;
        O ) outdir="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$bam" ] || [ -z "$ref" ] || [ -z "$contig" ] || [ -z "$id" ] || [ -z "$outdir" ]
then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the parameters are empty"
    [ -z "$bam" ] && echo "BAMFILE: ${bam}"
    [ -z "$ref" ] && echo "REF: ${ref}"
    [ -z "$contig" ] && echo "CONTIG: ${contig}"
    [ -z "$id" ] && echo "SAMPLEID: ${id}"
    [ -z "$outdir" ] && echo "OUTPUT: ${outdir}"
    helpFunction
fi

# echo all parameters
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
echo -e "\033[32m               Step1 extract SNPs                  \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "BAMFILE: ${bam}"
echo "REF: ${ref}"
echo "CELLBARCODETAG: ${cellbarcodetag}"
echo "CONTIG: ${contig}"
echo "SAMPLEID: ${id}"
echo "nprocs: ${nprocs}"
echo "ratioUser: ${ratioUser}"
echo "minbq: ${minbq}"
echo "minmq: ${minmq}"
echo "mincellnumber: ${mincellnumber}"
echo "minaltcellnumber: ${minaltcellnumber}"
echo "maxmismatch: ${maxmismatch}"
echo "minsequencingdepth: ${minsequencingdepth}"
echo "OUTPUT: ${outdir}"
echo ""
echo ""

# main steps
while read chrtmp starttmp endtmp
do
    tmpConfig=${chrtmp}_${starttmp}_${endtmp}
    echo "Extracting SNPs for ${tmpConfig}"
    python -u ${scMutrace_src_dir}/src/1_runStep1_ExtractSNP.py --bam ${bam} --ref ${ref} --cbtag ${cellbarcodetag} --contig ${tmpConfig} --nprocs ${nprocs} --nmcheck --outdir ${outdir} --id ${id} --minmq ${minmq} --minbq ${minbq} --mincell ${mincellnumber} \
        --minalt ${minaltcellnumber} --maxnm ${maxmismatch} --mindp ${minsequencingdepth} --ratio ${ratioUser}
done < ${contig}