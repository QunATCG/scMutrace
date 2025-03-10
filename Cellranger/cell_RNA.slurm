#!/bin/bash
#SBATCH -A naiss2024-22-89
#SBATCH --job-name cellCountATAC
#SBATCH --partition shared
#SBATCH -t 0-06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se


ml bioinfo-tools
ml cellranger/7.1.0

#refFA=/cfs/klemming/home/q/qunli/qundir/reference/homo/hg19_maskMito/hg19/
refFA=/cfs/klemming/home/q/qunli/qundir/reference/homo/GRCh38_maskMito/GRCh38/
workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/mito/fastq/
#CD34_500_Day08 CD34_500_Day14 CD34_800_Day08 CD34_800_Day14 CD34_800_Day20
for i in CD34_800_Day20
do
	sampleID=${i}
	cellranger count --id=${sampleID} --fastqs=${workingDir}/${sampleID}/  --transcriptome=${refFA} --sample=${sampleID} --localcores=16 --localmem=57
done
