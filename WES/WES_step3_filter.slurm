#!/bin/bash
#SBATCH -A naiss2024-22-344
#SBATCH --job-name WESgatk
#SBATCH --partition shared
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

module load systemdefault/1.0.0
module load gatk/4.5.0.0

workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/WES/
refFA=/cfs/klemming/home/q/qunli/qundir/reference/homo/hg19/hg19.fa
commBiallelic=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/gatkResource/small_exac_common_3_Addchr.vcf.gz

for sampleID in P2 P3 P4 P5 P6 P7 P8 P9 P10
do
	gatk --java-options "-Xmx16G" GetPileupSummaries -V ${commBiallelic} -L ${commBiallelic} -I ${workingDir}/bamFile/${sampleID}_normal_WES.bam -O ${workingDir}/table/${sampleID}.normal_pileups.table
	gatk --java-options "-Xmx16G" GetPileupSummaries -V ${commBiallelic} -L ${commBiallelic} -I ${workingDir}/bamFile/${sampleID}_cSCC_WES.bam -O ${workingDir}/table/${sampleID}.tumor_pileups.table
	gatk --java-options "-Xmx16G" CalculateContamination -I ${workingDir}/table/${sampleID}.tumor_pileups.table -matched ${workingDir}/table/${sampleID}.normal_pileups.table -O ${workingDir}/table/${sampleID}.contamination.table
	gatk --java-options "-Xmx16G" FilterMutectCalls -R ${refFA} -V ${workingDir}/somatic/${sampleID}_somatic.vcf.gz --contamination-table ${workingDir}/table/${sampleID}.contamination.table -O ${workingDir}/somatic/${sampleID}_somatic_filter.vcf.gz
done
