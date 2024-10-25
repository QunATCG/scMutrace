#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

# help info
helpFunction(){
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m                Step2 filter SNPs                  \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "Filter SNPs based on scMutrace Model"
    echo ""
    echo "Usage:"
    echo "$0 -V VCFDIR[str] -s SAMPLEID[str] -p nprocs[int] -k iteration[int] -e backp[float] -E errorp[float] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -V \t[vcf] vicf file directory to filter SNPs"
    echo -e "\t -s \t[id] sample ID"
    echo -e "\t -p \t[nprocs] number of processes"
    echo -e "\t -k \t[iteration] number of iterations"
    echo -e "\t -e \t[backp] background mutant cells threshold used to filter SNPs[0.03-0.06, default: 0.05]"
    echo -e "\t -E \t[errorp] background error threshold used to filter SNPs[default: 0.006]"
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
iteration=1000
backp=0.05
errorp=0.006


# get parameters
while getopts "V:s:p:k:e:E:O:h" opt
do
   case "$opt" in
        V ) vcfdir="$OPTARG" ;;
        s ) id="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        k ) iteration="$OPTARG" ;;
        e ) backp="$OPTARG" ;;
        E ) errorp="$OPTARG" ;;
        O ) outdir="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$vcfdir" ] || [ -z "$id" ] || [ -z "$outdir" ]
then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the parameters are empty"
    [ -z "$vcfdir" ] && echo "VCFDIR: ${vcfdir}"
    [ -z "$id" ] && echo "SAMPLEID: ${id}"
    [ -z "$outdir" ] && echo "outdir: ${outdir}"
    helpFunction
fi

# echo all parameters
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
echo -e "\033[32m               Step2 filter SNPs                   \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "VCFDIR: ${vcfdir}"
echo "SAMPLEID: ${id}"
echo "nprocs: ${nprocs}"
echo "iteration: ${iteration}"
echo "backp: ${backp}"
echo "errorp: ${errorp}"
echo "outdir: ${outdir}"
echo ""
echo ""

# main steps
python -u ${scMutrace_src_dir}/src/2_runStep0_FilterSNP.py --outdir ${outdir} --vcfdir ${vcfdir} --id ${id} --nprocs ${nprocs} --backp ${backp} --errorp ${errorp} --model3