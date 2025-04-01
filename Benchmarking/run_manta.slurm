#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name Manta
#SBATCH --partition shared
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

ml bioinfo-tools
ml manta/1.6.0

workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/MousescRNAscATACWGS/WGS/
refFA=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/mm10/Mus_musculus.GRCm38.dna.toplevel_filter.fa
mapFA=${workingDir}/CNV/map/Mus_musculus.GRCm38.dna.primary_assembly.fa.r101.s501.blacklist.gz

EPSRC1_bam=${workingDir}/EPSRC1_WGS/EPSRC1_WGS.bam
runDir=${workingDir}/CNV/EPSRC1/Manta/

mkdir -p ${runDir}

configManta.py --bam ${EPSRC1_bam} --referenceFasta ${refFA} --runDir ${runDir}

python ${runDir}/runWorkflow.py -j 8