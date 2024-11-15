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
    echo "$0 -b BAMFILE[str] -f REF[str] -v VCFPath[str] -s SAMPLEID[str] -c CELLBARCODEFILE[str] -C CELLBARCODETAG[str] -l BINSIZE[int] -L CLUSTERSIZE[int] -d MINDP[int] -D MINALTDP[int] -n MINCELLNUM[int] -N MINALTCELLNUM[int] -B LOWBQ[int] -q minBQ[int] -Q minMQ[str] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -b \t[bam path] bam path"
    echo -e "\t -f \t[ref path] reference path"
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
    echo -e "\t -C \t[cbtag] cell barcode tag (default CB)"
    echo -e "\t -q \t[min_bq] min base quality (default 0)"
    echo -e "\t -Q \t[min_mq] min mapping quality (default 255)"
    echo -e "\t -p \t[nprocs] number of processors to use (default 1)"
    echo -e "\t -h \t[help] display this help message"
    echo -e "\t Author: QunLi (liqun95@163.com)"
    exit 1
}

# Function to check and correct headers
function validate_and_correct_headers {
    local file="$1"
    local header=$(head -n 1 "$file")
    local first_column=$(echo "$header" | awk '{print $1}')
    local second_column=$(echo "$header" | awk '{print $2}')

    # Check and correct the first column header
    if [ "$first_column" != "Barcode" ]; then
        echo "WARNING: First column header is not 'Barcode'. Found '$first_column' instead. Adding 'Barcode' as the header."
        if [[ "$OSTYPE" == "darwin"* ]]; then
            sed -i '' '1s/^/Barcode\tCellType\n/' "$file"
        else
            sed -i '1s/^/Barcode\tCellType\n/' "$file"
        fi
    fi

    # Check and correct the second column header
    if [ "$second_column" != "CellType" ]; then
        echo "WARNING: Second column header is not 'CellType'. Found '$second_column' instead. Adding 'CellType' as the header and setting all values to 'CellType1'."
        awk 'BEGIN{OFS="\t"} NR==1{$2="CellType"} NR>1{$2="CellType1"}1' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
    fi
}

# Function to ensure two columns throughout the file
function validate_column_count {
    local file="$1"
    if [ $(awk '{print NF}' "$file" | sort -nu | tail -n 1) -ne 2 ]; then
        echo "ERROR: Some rows do not have two columns. Setting second column to 'CellType1' where missing."
        awk 'BEGIN{OFS="\t"} NR==1{$2="CellType"} NR>1{$2="CellType1"}1' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
    fi
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
min_bq=0
min_mq=255
cbtag="CB"
nprocs=1

# get parameters
while getopts "b:f:v:O:s:c:l:L:d:D:n:N:B:C:q:Q:p:h" opt
do
   case "$opt" in
        b ) bamfile="$OPTARG" ;;
        f ) ref="$OPTARG" ;;
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
        C ) cbtag="$OPTARG" ;;
        q ) min_bq="$OPTARG" ;;
        Q ) min_mq="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$vcfpath" ] || [ -z "$outdir" ] || [ -z "$cellbarcode" ] || [ -z "$id" ] || [ -z "$bamfile" ] || [ -z "$ref" ]
then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the parameters are empty"
    [ -z "$bamfile" ] && echo "BAMFILE: ${bamfile}"
    [ -z "$ref" ] && echo "REF: ${ref}"
    [ -z "$vcfpath" ] && echo "BAMFILE: ${vcfpath}"
    [ -z "$cellbarcode" ] && echo "REF: ${cellbarcode}"
    [ -z "$id" ] && echo "SAMPLEID: ${id}"
    [ -z "$outdir" ] && echo "OUTPUT: ${outdir}"
    helpFunction
fi

# Check if the barcode file exists
if [ ! -f "${cellbarcode}" ]; then
    echo "The barcode file does not exist. Please provide a valid barcode file."
    exit 1
fi

cp "$cellbarcode" "${cellbarcode}.bak"
validate_and_correct_headers "$cellbarcode"
validate_column_count "$cellbarcode"

echo "Barcode file header and format are correct."

# echo all parameters
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
echo -e "\033[32m            Step4 filter noisy sites               \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "BAMFILE: ${bamfile}"
echo "REF: ${ref}"
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
python -u ${scMutrace_src_dir}/src/3_runStep0_noisy.py --vcfpath ${vcfpath} --id ${id} --barcode ${cellbarcode} --binsize ${binsize} --clustersize ${clustersize} --dpsize ${minsequencingdepth} --altdpsize ${minaltdpsize} --cellsize ${mincellnumber} --altcellsize ${minaltcellnumber} --lowBaseCutoff ${lowBaseCutoff} --outdir ${outdir} --bam ${bamfile} --ref ${ref} --min_bq ${min_bq} --min_mq ${min_mq} --cbtag ${cbtag} --nprocs ${nprocs}