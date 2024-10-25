#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

helpFunction() {
    echo ""
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
    echo -e "\t -s \t[id:str] Sample ID."
    echo -e "\t -O \t[outdir:str] Output directory."
    echo -e "\t -g \t[contig:str] Contig information in tab-separated format: Chr\tstart\tend."
    echo "Advanced Arguments:"
    echo -e "\t -r \t[removeitems:str] Database files to be excluded (tab-separated format: path id)."
    echo -e "\t -i \t[includeitems:str] Database files to be included (tab-separated format: path id)."
    echo "Optional Arguments:"
    echo -e "\t -t \t[ncut:int] Number of read bases to be trimed cut (default: 10)."
    echo -e "\t -p \t[nprocs:int] Number of processes to use (default: 1)."
    echo -e "\t -q \t[minbq:int] Minimum base quality (default: 20)."
    echo -e "\t -Q \t[minmq:int] Minimum mapping quality (default: 255). Check your mapping software's quality score range (e.g., STAR, BWA, HISAT2)."
    echo -e "\t -n \t[mincellnumber:int] Minimum cell number (default: 3)."
    echo -e "\t -N \t[minaltcellnumber:int] Minimum altered cell number (default: 2)."
    echo -e "\t -m \t[maxmismatch:int] Maximum mismatch number (default: 4)."
    echo -e "\t -d \t[minsequencingdepth:int] Minimum sequencing depth (default: 3)."
    echo -e "\t -D \t[minaltdepth:int] Minimum alternative base depth (default: 2)."
    echo -e "\t -k \t[iteration:int] Number of iterations (default: 1000)."
    echo -e "\t -e \t[backp:float] Background mutant cells threshold to filter SNPs (0.03-0.06, default: 0.05)."
    echo -e "\t -E \t[errorp:float] Background error threshold to filter SNPs (0.001-0.006, default: 0.006)."
    echo -e "\t -l \t[binsize:int] Base number between two SNPs (default: 20)."
    echo -e "\t -L \t[clustersize:int] Number of SNPs in one cluster (default: 3)."
    echo -e "\t -B \t[lowBaseCutoff:int] Cutoff for low base quality (default: 20)."
    echo -e "\t -h \t[help] Display this help message."
    exit 1
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
minaltdepth=2
iteration=1000
backp=0.05
errorp=0.006
binsize=20
clustersize=3
lowBaseCutoff=20
removeitem="NotAvailable"
includeitem="NotAvailable"

# Check if no arguments are passed
if [ $# -eq 0 ]; then
    echo "Failed to parse arguments."
    helpFunction
fi

# get parameters
while getopts "b:f:c:s:O:g:r:i:t:p:q:Q:n:N:m:d:D:k:e:E:l:L:B:h" opt
do
    case "$opt" in
        b ) bam="$OPTARG" ;;
        f ) ref="$OPTARG" ;;
        c ) barcode="$OPTARG" ;;
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
        D ) minaltdepth="$OPTARG" ;;
        k ) iteration="$OPTARG" ;;
        e ) backp="$OPTARG" ;;
        E ) errorp="$OPTARG" ;;
        l ) binsize="$OPTARG" ;;
        L ) clustersize="$OPTARG" ;;
        B ) lowBaseCutoff="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;;
    esac
done

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
echo "CELLBARCODE: ${barcode}"
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
echo "minaltdepth: ${minaltdepth}"
echo "iteration: ${iteration}"
echo "backp: ${backp}"
echo "errorp: ${errorp}"
echo "binsize: ${binsize}"
echo "clustersize: ${clustersize}"
echo "lowBaseCutoff: ${lowBaseCutoff}"
echo ""
echo ""

# main steps
# step0_excludeBGCell.sh
echo "making the output directory.."${id}
echo "running step0_excludeBGCell.sh"
bash ${scMutrace_src_dir}/step0_excludeBGCell.sh -b ${bam} -f ${ref} -c ${barcode} -t ${ncut} -s ${id} -p ${nprocs} -O ${outdir}

# step1_extractSNP.sh
echo "running step1_extractSNP.sh"
bash ${scMutrace_src_dir}/step1_extractSNP.sh -b ${outdir}/ExcludeBG_${id}.sort.rmdupicard.bam -f ${ref} -g ${contig} -s ${id} -p ${nprocs} -q ${minbq} -Q ${minmq} -n ${mincellnumber} -N ${minaltcellnumber} -m ${maxmismatch} -d ${minsequencingdepth} -O ${outdir}

# step2_filterSNPModel.sh
echo "running step2_filterSNPModel.sh"
bash ${scMutrace_src_dir}/step2_filterSNPModel.sh -V ${outdir}/VCFPOS/ -s ${id} -p ${nprocs} -k ${iteration} -e ${backp} -E ${errorp} -O ${outdir}

# step3_filterPopulation.sh
echo "running step3_filterPopulation.sh"
bash ${scMutrace_src_dir}/step3_filterPopulation.sh -O ${outdir} -V ${outdir}/filterVCF/ -s ${id} -r ${removeitem} -i ${includeitem} -p ${nprocs}

# step4_filterNoisy.sh
echo "running step4_filterNoisy.sh"
bash ${scMutrace_src_dir}/step4_filterNoisy.sh -v ${outdir}/${id}_scmutrace.vcf --id ${id} --barcode ${barcode} --binsize ${binsize} --clustersize ${clustersize} --dpsize ${minsequencingdepth} --altdpsize ${minaltdpsize} --cellsize ${mincellnumber} --altcellsize ${minaltcellnumber} --lowBaseCutoff ${lowBaseCutoff} --outdir ${outdir}