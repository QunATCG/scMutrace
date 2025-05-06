#!/bin/bash
#SBATCH -A naiss2024-22-89
#SBATCH --job-name sam2bam
#SBATCH --partition core
#SBATCH -t 2-10:30:00
#SBATCH -n 8
#SBATCH --output /home/qunli/log/log.%j.%x
#SBATCH --error /home/qunli/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

module load bioinfo-tools
module load HISAT2/2.2.1
module load bowtie2/2.5.1
module load star/2.7.11a
module load rsem/1.3.3
module load bwa/0.7.17
module load samtools/1.19

THREADS=16
GENOMEINDEX=/home/qunli/qunr/refGenome/hg38/bow2Index/genome

FQ1=/home/qunli/qunr/callVariantsTestData/scATAC/1_ENCFF126AEH.fastq.gz
FQ2=/home/qunli/qunr/callVariantsTestData/scATAC/1_ENCFF989RBF.fastq.gz

#bowtie2 -X 2000 --threads $THREADS -x ${GENOMEINDEX} \
#	-1 $FQ1 -2 $FQ2 --sam-append-comment -k 5 -S bowtie2Out/1_ENCFF126AEH_ENCFF989RBF.sam


samtools view -@ $THREADS -Sb -o bowtie2Out/1_ENCFF126AEH_ENCFF989RBF.bam bowtie2Out/1_ENCFF126AEH_ENCFF989RBF.sam
