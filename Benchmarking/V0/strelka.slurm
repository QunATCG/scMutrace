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


# somatic
/cfs/klemming/home/q/qunli/Software/benchmarkingTools/strelka-2.9.2.centos6_x86_64/bin//configureStrelkaSomaticWorkflow.py \
--normalBam HCC1187BL.bam \
--tumorBam HCC1187C.bam \
--referenceFasta hg19.fa \
--indelCandidates ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz \
--runDir ${STRELKA_ANALYSIS_PATH}

# germline
/cfs/klemming/home/q/qunli/Software/benchmarkingTools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
--bam NA12878.bam \
--referenceFasta hg19.fa \
--runDir ${STRELKA_ANALYSIS_PATH}

/cfs/klemming/home/q/qunli/qunl/Qun/strelka/runWorkflow.py -m local -j 20
