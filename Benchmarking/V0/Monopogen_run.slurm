#!/bin/bash
#SBATCH -A naiss2024-5-32
#SBATCH --job-name Monopogen
#SBATCH --partition shared
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=200GB
#SBATCH --output /cfs/klemming/home/q/qunli/qunl/log/log.%j.%x
#SBATCH --error /cfs/klemming/home/q/qunli/qunl/log/error.%j.%x
#SBATCH --mail-type ALL
#SBATCH --mail-user qun.li@ki.se

source /cfs/klemming/home/q/qunli/.bashrc
conda activate scMutRNA
ml load PDC/23.12
ml load R/4.4.0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cfs/klemming/home/q/qunli/Software/benchmarkingTools/Monopogen/apps/
export R_LIBS_USER=/cfs/klemming/home/q/qunli/qunl/R/library/

MonopogenPath=/cfs/klemming/home/q/qunli/Software/benchmarkingTools/Monopogen/src/
MonopogenAPP=/cfs/klemming/home/q/qunli/Software/benchmarkingTools/Monopogen/apps/
LDIMPUTATION=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/1KG3Hg19/
workingDir=/cfs/klemming/home/q/qunli/qunl/Qun/callVariants/cSCC/20241022/
refFA=/cfs/klemming/home/q/qunli/qunl/Qun/publicData/reference/hg19/addChr/genome.fa
regionLST=${workingDir}/region.lst
nprocs=5

for sampleID in P2_cSCC P3_cSCC P4_cSCC P5_cSCC P6_cSCC P7_cSCC P8_cSCC P9_cSCC P10_cSCC
do
	outDir=${workingDir}/${sampleID}/Monopogen/
    barcode=${workingDir}/metadata/monopogen/${sampleID}_excludeBG.barcode.csv
    bamPath=${workingDir}/${sampleID}/scMutrace/ExcludeBG_${sampleID}.sort.rmdupicard.bam
	echo "somatic featureInfo"
    python  ${MonopogenPath}/Monopogen.py somatic -a ${MonopogenAPP} \
         -r ${workingDir}/region.lst  -t ${nprocs} -i ${outDir} -l ${barcode} \
          -s featureInfo  -g ${refFA}

    echo "somatic cellScan"
    python  ${MonopogenPath}/Monopogen.py somatic -a ${MonopogenAPP} \
          -r ${workingDir}/region.lst  -t ${nprocs} -i ${outDir}  -l ${barcode} \
          -s cellScan -g ${refFA}

    echo "LDrefinement"
    python  ${MonopogenPath}/Monopogen.py somatic -a ${MonopogenAPP} \
            -r ${workingDir}/region.lst -t ${nprocs} -i ${outDir}  -l ${barcode} \
             -s LDrefinement -g ${refFA}
done

#for sampleID in P2_cSCC P3_cSCC P4_cSCC P5_cSCC P6_cSCC P7_cSCC P8_cSCC P9_cSCC P10_cSCC
#do
#	outDir=${workingDir}/${sampleID}/Monopogen/
#	barcode=${workingDir}/metadata/monopogen/${sampleID}_excludeBG.barcode.csv
#	bamPath=${workingDir}/${sampleID}/scMutrace/ExcludeBG_${sampleID}.sort.rmdupicard.bam
#	#echo "${sampleID},${bamPath}" >> ${workingDir}/metadata/monopogen/bam_${sampleID}.lst
#	mkdir -p ${outDir}
#
#	echo "preProcess"
#	python ${MonopogenPath}/Monopogen.py preProcess -b ${workingDir}/metadata/monopogen/bam_${sampleID}.lst \
#		-o ${outDir} -a ${MonopogenAPP} -t 1
#
#	echo "main steps"
#	python ${MonopogenPath}/Monopogen.py germline  -a ${MonopogenAPP} \
#		-g ${refFA} -m 3 -s all -o ${outDir} \
#		-p ${LDIMPUTATION} -r ${regionLST} -t ${nprocs}
#
#	echo "somatic featureInfo"
#	python  ${MonopogenPath}/Monopogen.py somatic -a ${MonopogenAPP} \
#		-r ${workingDir}/region.lst  -t ${nprocs} -i ${outDir} -l ${barcode} \
#		-s featureInfo  -g ${refFA}
#
#	echo "somatic cellScan"
#	python  ${MonopogenPath}/Monopogen.py somatic -a ${MonopogenAPP} \
#		-r ${workingDir}/region.lst  -t ${nprocs} -i ${outDir}  -l ${barcode} \
#		-s cellScan -g ${refFA}
#
#	echo "LDrefinement"
#	python  ${MonopogenPath}/Monopogen.py somatic -a ${MonopogenAPP} \
#		-r ${workingDir}/region.lst -t ${nprocs} -i ${outDir}  -l ${barcode} \
#		-s LDrefinement -g ${refFA}
#done
