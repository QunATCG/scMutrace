#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name epiAneufinder
#SBATCH --partition shared
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

ml load PDC/23.12
ml load R/4.4.0

export R_LIBS_USER=/cfs/klemming/home/q/qunli/qundir/R/libraries/

Rscript /cfs/klemming/home/q/qunli/qunl/Qun/callVariants/MousescRNAscATACWGS/scATAC/run_epiAneufinder.R