#!/bin/bash
#SBATCH -A naiss2024-22-89
#SBATCH --job-name step3_filterPopulation
#SBATCH --partition main
#SBATCH -t 0-10:00:00
#SBATCH --nodes=1
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

/cfs/klemming/home/q/qunli/qunl/Qun/strelka/runWorkflow.py -m local -j 20
