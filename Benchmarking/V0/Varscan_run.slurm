#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name Varscan
#SBATCH --partition shared
#SBATCH -t 2-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

source /cfs/klemming/home/q/qunli/.bashrc
conda activate scMutRNA

module load samtools/1.20

path_to_varscan=/cfs/klemming/home/q/qunli/Software/benchmarkingTools/VarScan.v2.4.6.jar
faFile=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/hg19/noChr/genome.fa
workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/
# P2_cSCC P2_normal P3_cSCC P3_normal P4_cSCC P4_normal P5_cSCC P5_normal P6_cSCC P6_normal P7_cSCC P7_normal P8_cSCC P8_normal P9_cSCC P9_normal P10_cSCC P10_normal
for i in P4_normal P5_cSCC P5_normal P6_cSCC P6_normal P7_cSCC P7_normal P8_cSCC P8_normal P9_cSCC P9_normal P10_cSCC P10_normal
do
	echo $i
	sampleID=$i
	inBam=${workingDir}/${sampleID}/scMutrace/ExcludeBG_${sampleID}.sort.rmdupicard.bam
	outDir=${workingDir}/${sampleID}/Varscan/
	mkdir -p ${outDir}
	samtools mpileup -f $faFile $inBam | java -jar $path_to_varscan mpileup2snp --min-coverage 3 --min-reads2 1 --output-vcf 1 --strand-filter 0 --p-value 0.1 > ${outDir}/${sampleID}_varscan.vcf
done
