#!/bin/bash
#SBATCH -A naiss2024-22-1411
#SBATCH --job-name WESgatk
#SBATCH --partition shared
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

module load systemdefault/1.0.0
module load gatk/4.5.0.0

workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/WES/
refFA=/cfs/klemming/home/q/qunli/qundir/reference/homo/hg19/hg19.fa
afonlygnomad=/cfs/klemming/home/q/qunli/qunr/gatkData/af-only-gnomad.raw.sites.Addchr.vcf.gz
pon=${workingDir}/PoN/pon.vcf.gz
declare -A sampleNameNT
sampleNameNT["P2"]=22
sampleNameNT["P3"]=10
sampleNameNT["P4"]=20
sampleNameNT["P5"]=27
sampleNameNT["P6"]=28
sampleNameNT["P7"]=21
sampleNameNT["P9"]=26
sampleNameNT['P8']=1
sampleNameNT["P10"]=23
# P2 P3 P4 P5 P6 P7 P8 P9 P10
for sampleID in P9 P10
do
	NAMEtmp=${sampleNameNT[${sampleID}]}
	gatk --java-options "-Xmx70G" Mutect2 -R ${refFA} -I ${workingDir}/bamFile/${sampleID}_cSCC_WES.bam -I ${workingDir}/bamFile/${sampleID}_normal_WES.bam -tumor T${NAMEtmp} -normal N${NAMEtmp} --germline-resource ${afonlygnomad} --panel-of-normals ${pon} -O ${workingDir}/${sampleID}_somatic.vcf.gz 
done	

