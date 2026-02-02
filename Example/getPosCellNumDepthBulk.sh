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
    echo "get cell number and sequencing depth for each position"
    echo ""
    echo "Usage:"
    echo "$0 -b BAMFILE[str] -f REF[str] -g CONTIG[str] -s SAMPLEID[str] -p nprocs[int] -q minbq[int] -Q minmq[int] -m maxmismatch[int] -d minsequencingdepth[int] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -b \t[bam] bam file path to extract SNPs"
    echo -e "\t -f \t[ref] reference file path"
    echo -e "\t -g \t[contig] contig information, tab separated format for each line: Chr\tstart\tend"
    echo -e "\t -s \t[id] sample id"
    echo -e "\t -p \t[nprocs] number of processes"
    echo -e "\t -q \t[minbq] min base quality"
    echo -e "\t -Q \t[minmq] min mapping quality"
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
maxmismatch=4
minsequencingdepth=3


# get parameters
while getopts "b:f:g:s:p:q:Q:m:d:O:h" opt
do
   case "$opt" in
        b ) bam="$OPTARG" ;;
        f ) ref="$OPTARG" ;;
        g ) contig="$OPTARG" ;;
        s ) id="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        q ) minbq="$OPTARG" ;;
        Q ) minmq="$OPTARG" ;;
        m ) maxmismatch="$OPTARG" ;;
        d ) minsequencingdepth="$OPTARG" ;;
        O ) outdir="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;;
   esac
done

# Check if required parameters are empty
if [ -z "$bam" ] || [ -z "$ref" ] || [ -z "$id" ] || [ -z "$outdir" ] || [ -z "$contig" ]; then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the required parameters are empty:"
    [ -z "$bam" ] && echo "  BAMFILE is missing"
    [ -z "$ref" ] && echo "  REF is missing"
    [ -z "$id" ] && echo "  SAMPLEID is missing"
    [ -z "$contig" ] && echo "  CONTIG is missing"
    [ -z "$outdir" ] && echo "  OUTPUT directory is missing"
    helpFunction
fi

# check if output directory exists
if [ ! -d "$outdir" ]; then
    mkdir -p $outdir
fi

# echo parameters
echo "Parameters:"
echo "bam: $bam"
echo "ref: $ref"
echo "contig: $contig"
echo "id: $id"
echo "nprocs: $nprocs"
echo "minbq: $minbq"
echo "minmq: $minmq"
echo "maxmismatch: $maxmismatch"
echo "minsequencingdepth: $minsequencingdepth"
echo "outdir: $outdir"

# run python script
while read chrtmp starttmp endtmp
do
    tmpConfig=${chrtmp}_${starttmp}_${endtmp}
    echo "extract cell number and depth of ${tmpConfig}"
    python -u ${scMutrace_src_dir}/../src/getPosCellNumDepthBulk.py --bam ${bam} --ref ${ref} --contig ${tmpConfig} --minmq ${minmq} --minbq ${minbq} --mindp ${minsequencingdepth} --maxnm ${maxmismatch} --id ${id} --nprocs ${nprocs} --outdir ${outdir}
done < $contig