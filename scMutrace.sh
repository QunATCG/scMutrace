#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

# Help function
helpFunction(){
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m                OneStep Process                    \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "scMutrace OneStep Process"
    echo ""
    echo "Usage:"
    echo "$0 -b BAMFILE[str] -f REF[str] -c CELLBARCODE[str] -t ncut[int] -s SAMPLEID[str] -p nprocs[int] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    # step0_excludeBGCell.sh
    echo -e "\t -b \t[bam] bam file path to exclude background cells"
    echo -e "\t -f \t[ref] reference file path"
    echo -e "\t -c \t[barcode] cell barcode to be included"
    echo -e "\t -t \t[ncut] number of reads bases will be cut (default: 10)"
    echo -e "\t -s \t[id] sample id"
    echo -e "\t -p \t[nprocs] number of processes (default: 1)"
    echo -e "\t -O \t[outdir] output directory"
    # step1_extractSNP.sh
    echo -e "\t -g \t[contig] contig information, tab separated format for each line: Chr\tstart\tend"
    echo -e "\t -q \t[minbq] min base quality"
    echo -e "\t -Q \t[minmq] min mapping quality"
    echo -e "\t -n \t[mincellnumber] min cell number"
    echo -e "\t -N \t[minaltcellnumber] min alteration cell number"
    echo -e "\t -m \t[maxmismatch] number max mismatch number"
    echo -e "\t -d \t[minsequencingdepth] min sequecing depth"
    # step2_filterSNPModel.sh
    echo -e "\t -k \t[iteration] number of iterations"
    echo -e "\t -e \t[backp] background mutant cells threshold used to filter SNPs[0.03-0.06, default: 0.05]"
    echo -e "\t -E \t[errorp] background error threshold used to filter SNPs[default: 0.006]"
    # step3_filterPopulation.sh
    echo -e "\t -r \t[removeitems] files to be excluded, tab seperated format: path\tid"
    echo -e "\t -i \t[includeitems] files to be included, tab seperated format: path\tid"
    # step4_filterNoisy.sh
    echo -e "\t -l \t[binsize] base number between two snps (MultiAlleles)"
    echo -e "\t -L \t[clustersize] number of snps in one cluster"
    echo -e "\t -D \t[min altdpsize] min alternative bases depth"
    echo -e "\t -B \t[lowBaseCutoff] cutoff of low base quality (default 20)"
    echo -e "\t -h \t[help] display this help message"
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

# Check if no arguments are passed
if [ $# -eq 0 ]; then
    echo "Failed to parse arguments."
    helpFunction
fi

# get parameters
while getopts "b:f:c:t:s:p:O:g:q:Q:n:N:m:d:V:k:e:E:r:i:v:l:L:D:B:h" opt
do
    case "$opt" in
        b ) bam="$OPTARG" ;;
        f ) ref="$OPTARG" ;;
        c ) barcode="$OPTARG" ;;
        t ) ncut="$OPTARG" ;;
        s ) id="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        O ) outdir="$OPTARG" ;;

        g ) contig="$OPTARG" ;;
        q ) minbq="$OPTARG" ;;
        Q ) minmq="$OPTARG" ;;
        n ) mincellnumber="$OPTARG" ;;
        N ) minaltcellnumber="$OPTARG" ;;
        m ) maxmismatch="$OPTARG" ;;
        d ) minsequencingdepth="$OPTARG" ;;


        V ) casevcfdir="$OPTARG" ;;
        k ) iteration="$OPTARG" ;;
        e ) backp="$OPTARG" ;;
        E ) errorp="$OPTARG" ;;

        r ) removeitem="$OPTARG" ;;
        i ) includeitem="$OPTARG" ;;
        l ) binsize="$OPTARG" ;;
        L ) clustersize="$OPTARG" ;;
        D ) minaltdepth="$OPTARG" ;;
        B ) lowBaseCutoff="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
    esac
done

# Check if required parameters are empty
if [ -z "$bam" ] || [ -z "$ref" ] || [ -z "$barcode" ] || [ -z "$id" ] || [ -z "$outdir" ]; then
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m                  Log info                         \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo "Some or all of the required parameters are empty:"
    [ -z "$bam" ] && echo "  BAMFILE is missing"
    [ -z "$ref" ] && echo "  REF is missing"
    [ -z "$barcode" ] && echo "  CELLBARCODE is missing"
    [ -z "$id" ] && echo "  SAMPLEID is missing"
    [ -z "$outdir" ] && echo "  OUTPUT directory is missing"
    helpFunction  # Display the help info and exit
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
bash ${scMutrace_src_dir}/step4_filterNoisy.sh -v ${outdir}/${id}_scmutrace.vcf --id ${id} --barcode ${cellbarcode} --binsize ${binsize} --clustersize ${clustersize} --dpsize ${minsequencingdepth} --altdpsize ${minaltdpsize} --cellsize ${mincellnumber} --altcellsize ${minaltcellnumber} --lowBaseCutoff ${lowBaseCutoff} --outdir ${outdir}