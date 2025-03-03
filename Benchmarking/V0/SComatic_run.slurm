#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name scomatic
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

faFile=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/hg19/noChr/genome.fa
workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/
path_to_scomatic=/cfs/klemming/home/q/qunli/Software/benchmarkingTools/SComatic
editing=$path_to_scomatic/RNAediting/AllEditingSites.hg19.nochr.txt
PON=$path_to_scomatic/PoNs/PoN.scRNAseq.hg19.nochr.lift.tsv
nprocs=1

# P2_cSCC P3_cSCC P4_cSCC P5_cSCC P6_cSCC P7_cSCC P8_cSCC P9_cSCC P10_cSCC
for sampleID in P9_cSCC P10_cSCC
do
	inBam=${workingDir}/${sampleID}/scMutrace/ExcludeBG_${sampleID}.sort.rmdupicard.bam
	barcodeData=${workingDir}/metadata/scomatic/${sampleID}_excludeBG.barcode
	# step1
	echo ">>scomatic step1"
	SCOMATIC=${workingDir}/${sampleID}/SComatic/
	output_dir=$SCOMATIC/results
	mkdir -p $output_dir
	output_dir1=$output_dir/Step1_BamCellTypes
	mkdir -p $output_dir1
	python $path_to_scomatic/scripts/SplitBam/SplitBamCellTypes.py --bam $inBam \
		--meta ${barcodeData} \
		--id ${sampleID} \
		--n_trim 5 \
		--max_nM 5 \
		--max_NH 1 \
		--outdir $output_dir1

	# step2
	echo ">scomatic step2..."
	output_dir2=$output_dir/Step2_BaseCellCounts
	mkdir -p $output_dir2
	for bam in $(ls -d $output_dir1/*bam);do
		# Cell type
		echo $bam
		cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
		temp=$output_dir2/temp_${cell_type}
		mkdir -p $temp
		python $path_to_scomatic/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
			--ref $faFile \
			--chrom all \
			--out_folder $output_dir2 \
			--min_bq 20 \
			--tmp_dir $temp \
			--nprocs ${nprocs}
		rm -rf $temp
	done
	
	# step3
	echo ">scomatic step3..."
	output_dir3=$output_dir/Step3_BaseCellCountsMerged
	mkdir -p $output_dir3
	python $path_to_scomatic/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
		--outfile ${output_dir3}/${sampleID}.BaseCellCounts.AllCellTypes.tsv

	# step4
	echo ">scomatic step4..."
	output_dir4=$output_dir/Step4_VariantCalling
	mkdir -p $output_dir4
	python $path_to_scomatic/scripts/BaseCellCalling/BaseCellCalling.step1.py \
		--infile ${output_dir3}/${sampleID}.BaseCellCounts.AllCellTypes.tsv \
		--outfile ${output_dir4}/${sampleID} \
		--ref $faFile
	python $path_to_scomatic/scripts/BaseCellCalling/BaseCellCalling.step2.py \
		--infile ${output_dir4}/${sampleID}.calling.step1.tsv \
		--outfile ${output_dir4}/${sampleID} \
		--editing $editing \
		--pon $PON
done
