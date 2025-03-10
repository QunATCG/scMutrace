#!/bin/bash
#SBATCH -A naiss2023-22-1268
#SBATCH --job-name cellarc
#SBATCH --partition core
#SBATCH -t 0-24:00:00
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH --output /home/qunli/log/log.%j.%x
#SBATCH --error /home/qunli/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

module load bioinfo-tools
module load cellranger-ARC/2.0.2


cellranger-arc count --id=RT00374N \
	--reference=/home/qunli/qunr/refGenome/GRCh38/ \
	--libraries=/home/qunli/qundir/callvariants/GSE207334/RT00374N.library \
	--localcores=16 \
	--localmem=57