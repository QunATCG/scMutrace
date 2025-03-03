#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name samtools
#SBATCH --partition shared
#SBATCH -t 4-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

source /cfs/klemming/home/q/qunli/.bashrc
conda activate scMutRNA

ml load systemdefault/1.0.0
ml load samtools/1.20
ml load bcftools/1.20

refFA=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/hg19/noChr/genome.fa
workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/
for i in P2_cSCC P2_normal P3_cSCC P3_normal P4_cSCC P4_normal P5_cSCC P5_normal P6_cSCC P6_normal P7_cSCC P7_normal P8_cSCC P8_normal P9_cSCC P9_normal P10_cSCC P10_normal
do
	echo $i
	sampleID=$i
	inBam=${workingDir}/${sampleID}/scMutrace/ExcludeBG_${sampleID}.sort.rmdupicard.bam
	outDir=${workingDir}/${sampleID}/SAMtools/
	mkdir -p ${outDir}
	echo ${sampleID}
	bcftools mpileup -Q 20 -A -x -Ou -f ${refFA} ${inBam} | bcftools call -mv > ${outDir}/${sampleID}_samtools.vcf
done
