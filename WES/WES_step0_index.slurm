#!/bin/bash
#SBATCH -A naiss2024-22-1411
#SBATCH --job-name step0_exculdeBG
#SBATCH --partition shared
#SBATCH -t 0-03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

module load samtools/1.20

for sampleID in P2_cSCC P2_normal P3_cSCC P3_normal P4_cSCC P4_normal P5_cSCC P5_normal P6_cSCC P6_normal P7_cSCC P7_normal P8_cSCC P8_normal P9_cSCC P9_normal P10_cSCC P10_normal
do
	samtools index -@ 8 /cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/WES/${sampleID}_WES.bam
done	

