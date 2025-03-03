#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name freebayes
#SBATCH --partition shared
#SBATCH -t 1-08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

source /cfs/klemming/home/q/qunli/.bashrc
conda activate scMutRNA
ml PDC/23.12
ml R/4.4.0
ml systemdefault/1.0.0
ml freebayes/1.3.7

workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/

# fasta_generate_regions.py ref.fa.fai 100000 > targetRegion_freebayes.txt
for sampleID in P10_normal
do
	BAM=${workingDir}/${sampleID}/scMutrace/ExcludeBG_${sampleID}.sort.rmdupicard.bam
	OUT_DIR=${workingDir}/${sampleID}/freebayes/
	mkdir -p ${OUT_DIR}
	freebayes-parallel ${workingDir}/targetRegion_freebayes.txt 10 -f /cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/hg19/noChr/genome.fa ${BAM} > ${OUT_DIR}/${sampleID}_freebayes.vcf
done
