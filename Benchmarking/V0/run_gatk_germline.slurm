#!/bin/bash
#SBATCH -A naiss2024-22-1411
#SBATCH --job-name gatkGermLine
#SBATCH --partition shared
#SBATCH -t 0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

module load systemdefault/1.0.0
module load gatk/4.5.0.0
module load fastp/0.24.0
module load bwa/0.7.18

THREADS=8
workingDir=/cfs/klemming/home/q/qunli/qundir/A549_5K/WGS/
trimDir=/cfs/klemming/home/q/qunli/qundir/A549_5K/WGS/trimmed/
tmpDir=/cfs/klemming/home/q/qunli/qundir/A549_5K/WGS/tmp/
VCFOUT=/cfs/klemming/home/q/qunli/qundir/A549_5K/WGS/VCF/
refFa=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/hg38/addChr/bwaIndex/hg38.fa

DBSNP=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/gatkResource/Homo_sapiens_assembly38.dbsnp138.vcf
DBSNP1000G=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/gatkResource/1000G_phase1.snps.high_confidence.hg38.vcf
DBINDEL1000G=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/gatkResource/Homo_sapiens_assembly38.known_indels.vcf

mkdir -p ${tmpDir}
mkdir -p ${trimDir}
mkdir -p ${VCFOUT}

# trim fastq
for sampleID in DRR016843 DRR016844 DRR016845 DRR016846 ENCSR521ELB
do
	echo -e "\033[32m======QC......${sampleID}======\033[0m"
	R1_FASTQ=${workingDir}/rawFQ/${sampleID}_1.fastq.gz
	R2_FASTQ=${workingDir}/rawFQ/${sampleID}_2.fastq.gz
	fastp -i ${R1_FASTQ} -o ${trimDir}/${sampleID}_clean_1.fastq.gz -I ${R2_FASTQ} -O ${trimDir}/${sampleID}_clean_2.fastq.gz -h ${tmpDir}/${sampleID}_fastp.html -j ${tmpDir}/${sampleID}_fastp.json
done

# Mapping and Sort
for sampleID in DRR016843 DRR016844 DRR016845 DRR016846 ENCSR521ELB
do
	bwa mem -t ${THREADS} -R "@RG\tID:${sampleID}\tPL:illumina\tLB:WGS\tSM:${sampleID}" ${refFa} ${trimDir}/${sampleID}_clean_1.fastq.gz ${trimDir}/${sampleID}_clean_2.fastq.gz > ${trimDir}/${sampleID}.sam
	gatk --java-options "-Xmx32G" SortSam -I ${trimDir}/${sampleID}.sam --TMP_DIR ${tmpDir} -O ${trimDir}/${sampleID}.sort.bam -SO coordinate --CREATE_INDEX true
done

# Markdup
for sampleID in DRR016843 DRR016844 DRR016845 DRR016846 ENCSR521ELB
do
	gatk --java-options "-Xmx32G" MarkDuplicates -I ${trimDir}/${sampleID}.sort.bam  --TMP_DIR ${tmpDir} -O ${trimDir}/${sampleID}.sorted.marked.bam -M ${trimDir}/${sampleID}.metrics
done

# Calling
for sampleID in DRR016843 DRR016844 DRR016845 DRR016846 ENCSR521ELB
do
	gatk --java-options "-Xmx32G" BaseRecalibrator -R ${refFa} --tmp-dir ${tmpDir} -I ${trimDir}/${sampleID}.sorted.marked.bam -O ${trimDir}/${sampleID}.recal_data.table --known-sites ${DBSNP} --known-sites ${DBSNP1000G} --known-sites ${DBINDEL1000G}
	gatk --java-options "-Xmx32G" ApplyBQSR -R ${refFa} --tmp-dir ${tmpDir} -I ${trimDir}/${sampleID}.sorted.marked.bam --bqsr-recal-file ${trimDir}/${sampleID}.recal_data.table -O ${trimDir}/${sampleID}.sorted.marked.BQSR.bam
	gatk --java-options "-Xmx32G" HaplotypeCaller -R ${refFa} -I ${trimDir}/${sampleID}.sorted.marked.BQSR.bam -ERC GVCF --tmp-dir ${tmpDir} -O ${VCFOUT}/${sampleID}.g.vcf
	gatk --java-options "-Xmx32G" GenotypeGVCFs -R ${refFa} -V ${VCFOUT}/${sampleID}.g.vcf -O ${VCFOUT}/${sampleID}.vcf.gz
done

# filter
for sampleID in DRR016843 DRR016844 DRR016845 DRR016846 ENCSR521ELB
do
	gatk --java-options "-Xmx32G" VariantFiltration -R ${refFa} -V ${VCFOUT}/${sampleID}.vcf.gz -O ${VCFOUT}/${sampleID}_filter.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || QUAL < 30" --filter-name "QunBenchmarking"
done