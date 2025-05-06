#!/bin/bash
#SBATCH -A naiss2025-5-138
#SBATCH --job-name bowtie2
#SBATCH --partition shared
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

source /cfs/klemming/home/q/qunli/.bashrc
conda activate scMutRNA
ml bowtie2/2.5.4
ml samtools/1.20


THREADS=16
GENOMEINDEX=/cfs/klemming/home/q/qunli/qundir/reference/homo/hg19/igenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
faFile=/cfs/klemming/home/q/qunli/qundir/reference/homo/hg19/igenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa


workingDir=/cfs/klemming/home/q/qunli/qunr/GSE52529/
hg19_contig=${workingDir}/hg19.contig
removeItems_hg19=${workingDir}/exclude_hg19.txt
includeItems_hg19=${workingDir}/include_hg19.txt
processDir=${workingDir}/processed/

mkdir -p ${processDir}

for sampleID in
do
	echo ${sampleID}
	mkdir -p ${processDir}/${sampleID}/
	FQ1=${workingDir}/Rawfastq/${sampleID}/${sampleID}_1.fastq.gz
	FQ2=${workingDir}/Rawfastq/${sampleID}/${sampleID}_2.fastq.gz
	bowtie2 -p ${THREADS} -x ${GENOMEINDEX} -1 ${FQ1} -2 ${FQ2} -S ${processDir}/${sampleID}/${sampleID}.sam
	samtools view -@ ${THREADS} -Sb -o ${processDir}/${sampleID}/${sampleID}.bam ${processDir}/${sampleID}/${sampleID}.sam
	samtools sort -@ ${THREADS} ${processDir}/${sampleID}/${sampleID}.bam -o ${processDir}/${sampleID}/${sampleID}_sorted.bam
	samtools index -@ ${THREADS} ${processDir}/${sampleID}/${sampleID}_sorted.bam
	python /cfs/klemming/home/q/qunli/Private/scmutrace/src/0_runStep0_setCellBarcode.py --bam ${processDir}/${sampleID}/${sampleID}_sorted.bam --barcode ${sampleID} --outpath ${processDir}/${sampleID}/${sampleID}_sorted_CB.bam
	samtools sort -@ ${THREADS} ${processDir}/${sampleID}/${sampleID}_sorted_CB.bam -o ${processDir}/${sampleID}/${sampleID}_sorted_CB_sorted.bam
	samtools index -@ ${THREADS} ${processDir}/${sampleID}/${sampleID}_sorted_CB_sorted.bam
	rm ${processDir}/${sampleID}/${sampleID}.sam ${processDir}/${sampleID}/${sampleID}.bam ${processDir}/${sampleID}/${sampleID}_sorted.bam ${processDir}/${sampleID}/${sampleID}_sorted.bam.bai ${processDir}/${sampleID}/${sampleID}_sorted_CB.bam
	mv ${processDir}/${sampleID}/${sampleID}_sorted_CB_sorted.bam ${processDir}/${sampleID}/${sampleID}_sorted_CB.bam
	mv ${processDir}/${sampleID}/${sampleID}_sorted_CB_sorted.bam.bai ${processDir}/${sampleID}/${sampleID}_sorted_CB.bam.bai
	echo -e "Barcode\tCellType\n${sampleID}\tcelltype1" > ${processDir}/${sampleID}/cellbarcode.txt
	outDir=${processDir}/${sampleID}/scMutrace/
	mkdir -p ${outDir}
	/cfs/klemming/home/q/qunli/Private/scmutrace/scMutrace.sh -b ${processDir}/${sampleID}/${sampleID}_sorted_CB.bam -f ${faFile} -c ${processDir}/${sampleID}/cellbarcode.txt -s ${sampleID} -g ${hg19_contig} -r ${removeItems_hg19} -i ${includeItems_hg19} -u XX -q 20 -Q 23 -p 16 -e 0.045 -n 1 -N 1 -O ${outDir}
done