#!/bin/bash
scMutrace_src_dir=$(cd `dirname $0`; pwd)

# Help function
helpFunction(){
    echo ""
    echo -e "\033[32m###################################################\033[0m"
    echo -e "\033[32m             scMutrace Version 1.0                 \033[0m"
    echo -e "\033[32m             Step0 exclude BG cell                 \033[0m"
    echo -e "\033[32m###################################################\033[0m"
    echo ""
    echo "Description:"
    echo "Exclude background cells and trim sequencing reads"
    echo ""
    echo "Usage:"
    echo "$0 -b BAMFILE[str] -f REF[str] -c CELLBARCODE[str] -C CELLBARCODETAG[str] -t ncut[int] -s SAMPLEID[str] -p nprocs[int] -O OUTPUT[str]"
    echo ""
    echo "Options:"
    echo -e "\t -b \t[bam] bam file path to exclude background cells"
    echo -e "\t -f \t[ref] reference file path"
    echo -e "\t -x \t[barcode check] if barcode check, split bam based on cell barcodes (default: true)"
    echo -e "\t -c \t[barcode] cell barcode to be included"
    echo -e "\t -C \t[cellbarcodetag] cell barcode tag in bam file (default: CB)"
    echo -e "\t -t \t[ncut] number of reads bases will be cut (default: 10)"
    echo -e "\t -s \t[id] sample id"
    echo -e "\t -p \t[nprocs] number of processes (default: 1)"
    echo -e "\t -O \t[outdir] output directory"
    echo -e "\t -h \t[help] display this help message"
    echo -e "\t Author: QunLi (liqun95@163.com) "
    exit 1
}

# Default parameters
ncut=10
nprocs=1
cellbarcodetag="CB"
barcodecheck="true"

# Check if no arguments are passed
if [ $# -eq 0 ]; then
    echo "Failed to parse arguments."
    helpFunction
fi

# get parameters
while getopts "b:f:x:c:C:t:s:p:O:h" opt
do
   case "$opt" in
        b ) bam="$OPTARG" ;;
        f ) ref="$OPTARG" ;;
        x ) barcodecheck="$OPTARG" ;;
        c ) barcode="$OPTARG" ;;
        C ) cellbarcodetag="$OPTARG" ;;
        t ) ncut="$OPTARG" ;;
        s ) id="$OPTARG" ;;
        p ) nprocs="$OPTARG" ;;
        O ) outdir="$OPTARG" ;;
        h ) helpFunction ;;
        ? ) helpFunction ;; # Print helpFunction in case of invalid option
   esac
done

if [ "${barcodecheck}" != "true" ] && [ "${barcodecheck}" != "false" ]; then
    echo "Invalid value for -x. Please provide 'true' or 'false'."
    exit 1
fi

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
echo -e "\033[32m             Step0 exclude BG cell                 \033[0m"
echo -e "\033[32m###################################################\033[0m"
echo -e "\033[32mYour input parameters:\033[0m"
echo "BAMFILE: ${bam}"
echo "REF: ${ref}"
echo "barcode check: ${barcodecheck}"
echo "CELLBARCODE: ${barcode}"
echo "CELLBARCODETAG: ${cellbarcodetag}"
echo "ncut: ${ncut}"
echo "SAMPLEID: ${id}"
echo "nprocs: ${nprocs}"
echo "OUTPUT: ${outdir}"
echo ""

# main steps
# if cellbarcode check is true, then split bam file
if [ ${barcodecheck} == "true" ]; then
    echo "making the output directory.."${id}
    mkdir -p ${outdir}
    echo "Exclude background cells and trim reads.."${id}
    python -u ${scMutrace_src_dir}/src/0_runStep0_SplitBamOption3.py --bam ${bam} --ref ${ref} --barcode ${barcode} --splitCheck --ncut ${ncut} \
        --id ${id} --barcodetag ${cellbarcodetag} --nprocs ${nprocs} --outdir ${outdir}
    echo "samtools merge bam files.."${id}
    samtools merge -@ ${nprocs} -o ${outdir}/ExcludeBG_${id}.bam ${outdir}/${id}_*tmp.SMT.bam
    echo "samtools sort..${id}"
    samtools sort -@ ${nprocs} ${outdir}/ExcludeBG_${id}.bam -o ${outdir}/ExcludeBG_${id}.sort.bam
    echo "samtools index..${id}"
    samtools index -@ ${nprocs} ${outdir}/ExcludeBG_${id}.sort.bam
    rm ${outdir}/ExcludeBG_${id}.bam
    rm ${outdir}/*tmp.SMT.bam*
    samtools view -@ ${nprocs} -h ${outdir}/ExcludeBG_${id}.sort.bam | perl -pe 's/CB:Z:([A-Z]*)-1/CB:Z:$1/' | samtools view -@ ${nprocs} -b -o ${outdir}/ExcludeBG_${id}.sort.reCB.bam
    samtools index -@ ${nprocs} ${outdir}/ExcludeBG_${id}.sort.reCB.bam
    echo "picard rmdup.."${id}
    java -jar ${scMutrace_src_dir}/software/picard.jar MarkDuplicates --BARCODE_TAG ${cellbarcodetag} --REMOVE_DUPLICATES --INPUT ${outdir}/ExcludeBG_${id}.sort.reCB.bam --METRICS_FILE ${outdir}/ExcludeBG_${id}.picard_dup_metrics.txt --OUTPUT ${outdir}/ExcludeBG_${id}.sort.rmdupicard.bam
    echo "samtools index rmdup.."${id}
    samtools index -@ ${nprocs} ${outdir}/ExcludeBG_${id}.sort.rmdupicard.bam
    rm ${outdir}/ExcludeBG_${id}.sort.reCB.bam ${outdir}/ExcludeBG_${id}.sort.reCB.bam.bai
    echo "Done!"${id}
elif [ ${barcodecheck} == "false" ]; then
    samtools view -@ ${nprocs} -h ${bam} | perl -pe 's/CB:Z:([A-Z]*)-1/CB:Z:$1/' | samtools view -@ ${nprocs} -b -o ${outdir}/ExcludeBG_${id}.sort.reCB.bam
    samtools index -@ ${nprocs} ${outdir}/ExcludeBG_${id}.sort.reCB.bam
    echo "picard rmdup.."${id}
    java -jar ${scMutrace_src_dir}/software/picard.jar MarkDuplicates --BARCODE_TAG ${cellbarcodetag} --REMOVE_DUPLICATES --INPUT ${outdir}/ExcludeBG_${id}.sort.reCB.bam --METRICS_FILE ${outdir}/ExcludeBG_${id}.picard_dup_metrics.txt --OUTPUT ${outdir}/ExcludeBG_${id}.sort.rmdupicard.bam
    echo "samtools index rmdup.."${id}
    samtools index -@ ${nprocs} ${outdir}/ExcludeBG_${id}.sort.rmdupicard.bam
    rm ${outdir}/ExcludeBG_${id}.sort.reCB.bam ${outdir}/ExcludeBG_${id}.sort.reCB.bam.bai
    echo "Done!"${id}
fi