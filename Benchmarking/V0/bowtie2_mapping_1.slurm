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

ml bowtie2/2.5.4
ml samtools/1.20


THREADS=16
GENOMEINDEX=/cfs/klemming/home/q/qunli/qundir/reference/homo/hg19/igenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

FQ1=/home/qunli/qunr/callVariantsTestData/scATAC/1_ENCFF126AEH.fastq.gz
FQ2=/home/qunli/qunr/callVariantsTestData/scATAC/1_ENCFF989RBF.fastq.gz


bowtie2 -p ${THREADS} -x ${GENOMEINDEX} -1 read1.fastq -2 read2.fastq -S output.sam


samtools view -@ $THREADS -Sb -o bowtie2Out/1_ENCFF126AEH_ENCFF989RBF.bam bowtie2Out/1_ENCFF126AEH_ENCFF989RBF.sam


bowtie2/2.5.4