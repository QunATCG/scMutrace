#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name scAllele
#SBATCH --partition shared
#SBATCH -t 02-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

source /cfs/klemming/home/q/qunli/.bashrc
conda activate compareTools

faFile=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/hg19/noChr/genome.fa
outDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/

#for i in P10DxR1 P10DxR2 P10DxR3 P10DxR4 P10ReR1 P10ReR2
# P2_cSCC P2_normal P3_cSCC P3_normal P4_cSCC P4_normal P5_cSCC P5_normal P6_cSCC P6_normal P7_cSCC P7_normal P8_cSCC P8_normal P9_cSCC P9_normal P10_cSCC P10_normal
#for sampleID in P2_cSCC P2_normal P3_cSCC P3_normal P8_normal P10_cSCC
#do
#        echo ${sampleID}
# 	 mkdir -p ${outDir}/${sampleID}/scAllele/
#        bampath=${outDir}/${sampleID}/scMutrace//ExcludeBG_${sampleID}.sort.rmdupicard.bam
#        for j in {1..22}
#        do
#                scAllele  -b $bampath -g $faFile -o ${outDir}/${sampleID}/scAllele/${sampleID}_${j} --AC=3 --DP=3 -c $j
#        done
#done

for sampleID in P8_normal
do
        echo ${sampleID}
        mkdir -p ${outDir}/${sampleID}/scAllele/
        bampath=${outDir}/${sampleID}/scMutrace//ExcludeBG_${sampleID}.sort.rmdupicard.bam
        for j in {11..22}
        do
                scAllele  -b $bampath -g $faFile -o ${outDir}/${sampleID}/scAllele/${sampleID}_${j} --AC=3 --DP=3 -c $j
        done
done

for sampleID in P10_cSCC
do
        echo ${sampleID}
        mkdir -p ${outDir}/${sampleID}/scAllele/
        bampath=${outDir}/${sampleID}/scMutrace//ExcludeBG_${sampleID}.sort.rmdupicard.bam
        for j in {1..10}
        do
                scAllele  -b $bampath -g $faFile -o ${outDir}/${sampleID}/scAllele/${sampleID}_${j} --AC=3 --DP=3 -c $j
        done
done

