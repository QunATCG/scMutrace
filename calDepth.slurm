#!/bin/bash
#SBATCH -A naiss2024-22-89
#SBATCH --job-name calDepth
#SBATCH --partition shared
#SBATCH -t 0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/WES/bamFile/

for tmpBam in P2_cSCC_WES.bam P3_cSCC_WES.bam P4_cSCC_WES.bam P5_cSCC_WES.bam P6_cSCC_WES.bam P7_cSCC_WES.bam P8_cSCC_WES.bam P9_cSCC_WES.bam P10_cSCC_WES.bam P2_normal_WES.bam P3_normal_WES.bam P4_normal_WES.bam P5_normal_WES.bam P6_normal_WES.bam P7_normal_WES.bam P8_normal_WES.bam P9_normal_WES.bam P10_normal_WES.bam
do
    inputBam=${workingDir}/${tmpBam}
    /cfs/klemming/home/q/qunli/Private/scmutrace/calDepth.sh ${inputBam}
done