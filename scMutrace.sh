#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

helpFunction() {
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m==================================================\033[0m"
    echo -e "\033[32m             __  __       _                       \033[0m"
    echo -e "\033[32m   ___  ___|  \\/  |_   _| |_ _ __ __ _  ___ ___  \033[0m"
    echo -e "\033[32m  / __|/ __| |\\/| | | | | __| '__/ _\` |/ __/ _ \\ \033[0m"
    echo -e "\033[32m  \\__ \\ (__| |  | | |_| | |_| | | (_| | (_|  __/ \033[0m"
    echo -e "\033[32m  |___/\\___|_|  |_|\\__,_|\\__|_|  \\__,_|\\___\\___| \033[0m"
    echo -e "\033[32m==================================================\033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m                OneStep Process                    \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "scMutrace is a single-cell mutation tracing tool designed to identify and analyze mutations from"
    echo "single-cell sequencing data. It offers an efficient, one-step process for background filtering,"
    echo "SNP extraction, noise filtering, and population structure analysis." 
    echo ""
    echo "Usage:"
    echo "$0 -b BAMFILE[str] -f REF[str] -c CELLBARCODE[str] -s SAMPLEID[str] -g CONTIG[str] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo "Required Arguments:"
    echo -e "\t -b \t[bam:str] Path to the BAM file to exclude background cells."
    echo -e "\t -f \t[ref:str] Path to the reference file."
    echo -e "\t -c \t[barcode:str] Cell barcode file path (tab-separated with header: Barcode CellType). If CellType is not provided, all cells will be defined as 'CellType1'."
    echo -e "\t -C \t[barcodetag:str] Cell barcode tag in the BAM file (default: CB)."
    echo -e "\t -s \t[id:str] Sample ID."
    echo -e "\t -O \t[outdir:str] Output directory."
    echo -e "\t -g \t[contig:str] Contig information in tab-separated format: Chr\tstart\tend."
    echo "Advanced Arguments:"
    echo -e "\t -r \t[removeitems:str] Database files to be excluded (tab-separated format: path id)."
    echo -e "\t -i \t[includeitems:str] Database files to be included (tab-separated format: path id)."
    echo "Optional Arguments:"
    echo -e "\t -x \t[barcodecheck:str] if barcode check, split bam based on cell barcodes (false or true, default: true.)."
    echo -e "\t -u \t[umi:str] umi to be used to remove duplicated reads (default: CB, if XX, no umi)."
    echo -e "\t -t \t[ncut:int] Number of read bases to be trimed cut (default: 10)."
    echo -e "\t -p \t[nprocs:int] Number of processes to use (default: 1)."
    echo -e "\t -q \t[minbq:int] Minimum base quality (default: 20)."
    echo -e "\t -Q \t[minmq:int] Minimum mapping quality (default: 255, scATAC: 30). Check your mapping software's quality score range (e.g., STAR, BWA, HISAT2)."
    echo -e "\t -n \t[mincellnumber:int] Minimum cell number (default: 3)."
    echo -e "\t -N \t[minaltcellnumber:int] Minimum altered cell number (default: 2)."
    echo -e "\t -m \t[maxmismatch:int] Maximum mismatch number (default: 4)."
    echo -e "\t -d \t[minsequencingdepth:int] Minimum sequencing depth (default: 3)."
    echo -e "\t -D \t[minaltsequencingdepth:int] Minimum alternative base depth (default: 2)."
    echo -e "\t -k \t[iteration:int] Number of iterations (default: 1000)."
    echo -e "\t -e \t[backp:float] Background mutant cells threshold to filter SNPs (0.03-0.06, default: 0.05)."
    echo -e "\t -E \t[errorp:float] Background error threshold to filter SNPs (0.001-0.006, default: 0.006)."
    echo -e "\t -l \t[binsize:int] Base number between two SNPs (default: 20)."
    echo -e "\t -L \t[clustersize:int] Number of SNPs in one cluster (default: 3)."
    echo -e "\t -B \t[lowBaseCutoff:int] Cutoff for low base quality (default: 20)."
    echo -e "\t -R \t[ratio] cutoff to filter low freq mutations."
    echo -e "\t -h \t[help] Display this help message."
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

# Default parameters
ncut=10
nprocs=1
minbq=20
minmq=255
mincellnumber=3
minaltcellnumber=2
maxmismatch=4
minsequencingdepth=3
minaltsequencingdepth=2
iteration=1000
backp=0.05
errorp=0.006
binsize=20
clustersize=3
lowBaseCutoff=20
removeitem="NotAvailable"
includeitem="NotAvailable"
cellbarcodetag="CB"
barcodecheck="true"
ratioUser=0.005
umi="CB"

# Check if no arguments are passed
if [ $# -eq 0 ]; then
    echo "Failed to parse arguments."
    helpFunction
fi

# get parameters
while getopts "b:f:x:c:C:s:O:g:r:i:t:p:q:Q:R:n:N:m:d:D:k:e:E:l:L:B:u:h" opt
do
    case "$opt" in
        b ) bam="$OPTARG" ;;
        f ) ref="$OPTARG" ;;
        x ) barcodecheck="$OPTARG" ;;
        c ) barcode="$OPTARG" ;;
        C ) cellbarcodetag="$OPTARG" ;;
        s ) id="$OPTARG" ;;
        O ) outdir="$OPTARG" ;;
        g ) contig="$OPTARG" ;;
        r ) removeitem="$OPTARG" ;;
        i ) includeitem="$OPTARG" ;;
        t ) ncut="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        q ) minbq="$OPTARG" ;;
        Q ) minmq="$OPTARG" ;;
        n ) mincellnumber="$OPTARG" ;;
        N ) minaltcellnumber="$OPTARG" ;;
        m ) maxmismatch="$OPTARG" ;;
        d ) minsequencingdepth="$OPTARG" ;;
        D ) minaltsequencingdepth="$OPTARG" ;;
        k ) iteration="$OPTARG" ;;
        e ) backp="$OPTARG" ;;
        E ) errorp="$OPTARG" ;;
        l ) binsize="$OPTARG" ;;
        L ) clustersize="$OPTARG" ;;
        B ) lowBaseCutoff="$OPTARG" ;;
        u ) umi="$OPTARG" ;;
        R ) ratioUser="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;;
    esac
done

cat << 'EOF'
==================================================
             __  __       _
   ___  ___|  \/  |_   _| |_ _ __ __ _  ___ ___
  / __|/ __| |\/| | | | | __| '__/ _` |/ __/ _ \
  \__ \ (__| |  | | |_| | |_| | | (_| | (_|  __/
  |___/\___|_|  |_|\__,_|\__|_|  \__,_|\___\___|
==================================================
EOF

# check barcodecheck
if [ "${barcodecheck}" != "true" ] && [ "${barcodecheck}" != "false" ]; then
    echo "Invalid value for -x. Please provide 'true' or 'false'."
    exit 1
fi

# Check if required parameters are empty
if [ -z "$bam" ] || [ -z "$ref" ] || [ -z "$barcode" ] || [ -z "$id" ] || [ -z "$outdir" ] || [ -z "$contig" ]; then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the required parameters are empty:"
    [ -z "$bam" ] && echo "  BAMFILE is missing"
    [ -z "$ref" ] && echo "  REF is missing"
    [ -z "$barcode" ] && echo "  CELLBARCODE is missing"
    [ -z "$id" ] && echo "  SAMPLEID is missing"
    [ -z "$contig" ] && echo "  CONTIG is missing"
    [ -z "$outdir" ] && echo "  OUTPUT directory is missing"
    helpFunction
fi

# print all parameters
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
echo -e "\033[32m                OneStep Process                    \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "BAMFILE: ${bam}"
echo "REF: ${ref}"
echo "CELLBARCODECheck: ${barcodecheck}"
echo "CELLBARCODE: ${barcode}"
echo "CELLBARCODETAG: ${cellbarcodetag}"
echo "ncut: ${ncut}"
echo "SAMPLEID: ${id}"
echo "nprocs: ${nprocs}"
echo "OUTPUT: ${outdir}"
echo "CONTIG: ${contig}"
echo "removeitem: ${removeitem}"
echo "includeitem: ${includeitem}"
echo "minbq: ${minbq}"
echo "minmq: ${minmq}"
echo "mincellnumber: ${mincellnumber}"
echo "minaltcellnumber: ${minaltcellnumber}"
echo "maxmismatch: ${maxmismatch}"
echo "minsequencingdepth: ${minsequencingdepth}"
echo "minaltsequencingdepth: ${minaltsequencingdepth}"
echo "iteration: ${iteration}"
echo "backp: ${backp}"
echo "errorp: ${errorp}"
echo "binsize: ${binsize}"
echo "clustersize: ${clustersize}"
echo "lowBaseCutoff: ${lowBaseCutoff}"
echo "UMI: ${umi}"
echo "ratioUser: ${ratioUser}"
echo ""
echo ""

# Check if the output directory exists
if [ ! -d "${outdir}" ]; then
    mkdir -p "${outdir}"
fi

# Check if the barcode file exists
if [ ! -f "${barcode}" ]; then
    echo "The barcode file does not exist. Please provide a valid barcode file."
    exit 1
fi

# Check if the reference file exists
if [ ! -f "${ref}" ]; then
    echo "The reference file does not exist. Please provide a valid reference file."
    exit 1
fi

# Check if the BAM file exists
if [ ! -f "${bam}" ]; then
    echo "The BAM file does not exist. Please provide a valid BAM file."
    exit 1
fi

# Check if the contig file exists
if [ ! -f "${contig}" ]; then
    echo "The contig file does not exist. Please provide a valid contig file."
    exit 1
fi

cp "$barcode" "${barcode}.bak"
validate_and_correct_headers "$barcode"
validate_column_count "$barcode"

echo "Barcode file header and format are correct."

# main steps

# step0_excludeBGCell.sh
echo "making the output directory.."${id}
echo "running step0_excludeBGCell.sh"
${scMutrace_src_dir}/step0_excludeBGCell.sh -b ${bam} -f ${ref} \
    -u ${umi} -x ${barcodecheck} -c ${barcode} -C ${cellbarcodetag} \
    -g ${contig} -t ${ncut} -s ${id} -p ${nprocs} -O ${outdir}

# if step0_excludeBGCell.sh failed, exit
if [ $? -ne 0 ]; then
    echo "Failed to run step0_excludeBGCell.sh. Exiting..."
    exit 1
fi

# step1_extractSNP.sh
echo "running step1_extractSNP.sh"
${scMutrace_src_dir}/step1_extractSNP.sh -b ${outdir}/ExcludeBG_${id}.sort.rmdupicard.bam -f ${ref} -C ${cellbarcodetag} \
    -g ${contig} -s ${id} -p ${nprocs} -R ${ratioUser} -q ${minbq} -Q ${minmq} \
    -n ${mincellnumber} -N ${minaltcellnumber} \
    -m ${maxmismatch} -d ${minsequencingdepth} -D ${minaltsequencingdepth} \
    -O ${outdir}

# if step1_extractSNP.sh failed, exit
if [ $? -ne 0 ]; then
    echo "Failed to run step1_extractSNP.sh. Exiting..."
    exit 1
fi

# step2_filterSNPModel.sh
echo "running step2_filterSNPModel.sh"
${scMutrace_src_dir}/step2_filterSNPModel.sh -V ${outdir}/VCFPOS/ -s ${id} -p ${nprocs} -k ${iteration} -e ${backp} -E ${errorp} -O ${outdir}

# if step2_filterSNPModel.sh failed, exit
if [ $? -ne 0 ]; then
    echo "Failed to run step2_filterSNPModel.sh. Exiting..."
    exit 1
fi

# step3_filterPopulation.sh
echo "running step3_filterPopulation.sh"
${scMutrace_src_dir}/step3_filterPopulation.sh -O ${outdir} -V ${outdir}/filterVCF/ -s ${id} -r ${removeitem} -i ${includeitem} -p 1

# if step3_filterPopulation.sh failed, exit
if [ $? -ne 0 ]; then
    echo "Failed to run step3_filterPopulation.sh. Exiting..."
    exit 1
fi

# step4_filterNoisy.sh
echo "running step4_filterNoisy.sh"
${scMutrace_src_dir}/step4_filterNoisy.sh -b ${outdir}/ExcludeBG_${id}.sort.rmdupicard.bam \
    -f ${ref} -v ${outdir}/${id}_scmutrace.vcf -O ${outdir}/ \
    -s ${id} -c ${barcode} -l ${binsize} -L ${clustersize} \
    -d ${minsequencingdepth} -D ${minaltsequencingdepth} \
    -n ${mincellnumber} -N ${minaltcellnumber} -B ${lowBaseCutoff} \
    -C ${cellbarcodetag} -q 0 -Q ${minmq} -p ${nprocs}

# if step4_filterNoisy.sh failed, exit
if [ $? -ne 0 ]; then
    echo "Failed to run step4_filterNoisy.sh. Exiting..."
    exit 1
fi