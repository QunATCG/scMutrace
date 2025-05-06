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

workingDir=/cfs/klemming/home/q/qunli/qunr/GSE52529/
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
	samtools index -@ ${THREADS} ${processDir}/${sampleID}/${sampleID}.bam 
	python /cfs/klemming/home/q/qunli/Private/scmutrace/src/0_runStep0_setCellBarcode.py --bam ${processDir}/${sampleID}/${sampleID}.bam --barcode ${sampleID} --outpath ${processDir}/${sampleID}/${sampleID}_CB.bam 
	samtools sort -@ ${THREADS} ${processDir}/${sampleID}/${sampleID}_CB.bam -o ${processDir}/${sampleID}/${sampleID}_CB_sorted.bam 
	samtools index -@ ${THREADS} ${processDir}/${sampleID}/${sampleID}_CB_sorted.bam
done