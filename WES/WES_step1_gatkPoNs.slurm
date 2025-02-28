#!/bin/bash
#SBATCH -A naiss2024-22-1411
#SBATCH --job-name gatkPoN
#SBATCH --partition shared
#SBATCH -t 0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

module load systemdefault/1.0.0
module load gatk/4.5.0.0

workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/WES/
refFA=/cfs/klemming/home/q/qunli/qundir/reference/homo/hg19/hg19.fa
intervalLit=/cfs/klemming/home/q/qunli/qunr/gatkData/hg19.genome.interval_list
afonlygnomad=/cfs/klemming/home/q/qunli/qunr/gatkData/af-only-gnomad.raw.sites.vcf
SCRATCHDIR=/tmp/
PONDB=${SCRATCHDIR}pon_db/
#for sampleID in P2 P3 P8
#do
#	gatk --java-options "-Xmx40G" Mutect2 -R ${refFA} -I ${workingDir}/bamFile/${sampleID}_normal_WES.bam -max-mnp-distance 0 -O ${workingDir}/PoN/${sampleID}_normal_WES.vcf.gz
#done
gatk GenomicsDBImport -R ${refFA} -L ${intervalLit} --genomicsdb-workspace-path ${PONDB} -V ${workingDir}/PoN/P2_normal_WES.vcf.gz -V ${workingDir}/PoN/P3_normal_WES.vcf.gz -V ${workingDir}/PoN/P4_normal_WES.vcf.gz -V ${workingDir}/PoN/P5_normal_WES.vcf.gz -V ${workingDir}/PoN/P6_normal_WES.vcf.gz -V ${workingDir}/PoN/P7_normal_WES.vcf.gz -V ${workingDir}/PoN/P8_normal_WES.vcf.gz -V ${workingDir}/PoN/P9_normal_WES.vcf.gz -V ${workingDir}/PoN/P10_normal_WES.vcf.gz

gatk CreateSomaticPanelOfNormals -R ${refFA} -V gendb://${PONDB} -O ${workingDir}/pon.vcf.gz
