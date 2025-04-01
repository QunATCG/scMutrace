#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name Delly
#SBATCH --partition shared
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=200GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

ml bioinfo-tools
ml delly/1.1.6

workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/MousescRNAscATACWGS/WGS/
refFA=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/mm10/mm10_nochr.fa
mapFA=${workingDir}/CNV/map/Mus_musculus.GRCm38.dna.primary_assembly.fa.r101.s501.blacklist.gz


for sampleID in EPSRC1 EPSRC2 EPSRC4
do
	bamIn=${workingDir}/${sampleID}_WGS/${sampleID}_WGS.bam
	delly cnv -g ${refFA} -m ${mapFA} -c ${workingDir}/CNV/${sampleID}_out.cov.gz -o ${workingDir}/CNV/${sampleID}_cnv.bcf ${bamIn}
done