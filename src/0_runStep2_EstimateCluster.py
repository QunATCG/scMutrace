import os
import timeit
import argparse
import subprocess
from py_utils import *
import scipy.stats as stats
import multiprocessing as mp

def initialize_parser():
    parser = argparse.ArgumentParser(description='filter SNP based on public database')
    parser.add_argument('--outdir', type=str, default='', help='out file path', required = True)
    parser.add_argument('--bam', type=str, default='', help='bam file path', required = True)
    parser.add_argument('--ref', type=str, default='', help='reference file path', required = True)
    parser.add_argument('--vcfpath', type=str, default='', help='Vcf path', required = True)
    parser.add_argument('--id', type=str, default='', help='sample id', required = True)
    parser.add_argument('--binsize', type=int, default=10, help='base number between two snps (MultiAlleles)', required = False)
    parser.add_argument('--dpsize', type=int, default=10, help='min depth', required = False)
    parser.add_argument('--altdpsize', type=int, default=5, help='min alternative depth', required = False)
    parser.add_argument('--cellsize', type=int, default=10, help='min cell number', required = False)
    parser.add_argument('--altcellsize', type=int, default=5, help='min mut cell number', required = False)
    parser.add_argument('--min_bq', type=int, default=0, help='cutoff base quality', required = False)
    parser.add_argument('--min_mq', type=int, default=255, help='cutoff mapping quality', required = False)
    parser.add_argument('--cbtag', type=str, default='CB', help='cell barcode tag', required = False)
    parser.add_argument('--nprocs', type=int, default=1, help='number of threads', required = False)
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

def assign_alt_gt(refBase, sitemutinfo):
    indelBase = {"D", "I"}
    snpBase = {"A", "T", "C", "G"}
    basicBase = {"A", "T", "C", "G", "D", "I"}
    countDic = {base: 0 for base in basicBase}
    check_base = basicBase - {refBase}
    sitemutinfoList = list(map(int, sitemutinfo.split(":")))
    for base, count in zip(["A", "T", "C", "G", "D", "I"], sitemutinfoList):
        countDic[base] += count
    gtList = [base for base in check_base if countDic[base] != 0]
    indelintersection = indelBase & set(gtList)
    snpintersection = snpBase & set(gtList)
    altGT = ",".join(f"{base}" for base in snpintersection | indelintersection)
    return altGT

def check_noisy_bg(targetVcfPath, altDp, altCell):
    data = open(targetVcfPath, "r")
    mutCount = 0
    for lineeach in data:
        if lineeach.startswith("#"):
            continue
        line = lineeach.strip().split("\t")
        if int(line[5]) >= altDp and int(line[8]) >= altCell:
            mutCount = mutCount + 1
    return mutCount


def count_cluster_snp(vcfPathIn, vcfPathOut, bamPath, refPath, outdirTmp, min_bq, min_mq, min_cov, min_alt, min_cell, min_altcell, max_nm, cbtag, binsize, id):
    if os.path.exists(vcfPathOut):
        try:
            subprocess.run(['rm', vcfPathOut], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error removing file {vcfPathOut}: {e}")
    with open(vcfPathIn, "r") as vcfIn, open(vcfPathOut, "w") as vcfOut:
        for lineeach in vcfIn:
            linesets = lineeach.strip().split("\t")
            chr = linesets[0].strip()
            pos = linesets[1].strip()
            target_extractsnp_fromBAM(bamPath, refPath, outdirTmp, [chr, pos, pos], min_bq, min_mq, True, True, True, False, min_cov, min_alt, min_cell, max_nm, cbtag, binsize, id)
            targetVcfPath = os.path.join(outdirTmp, f"VCF_{chr}_{pos}_{pos}_{id}_random.vcf")
            mutCount= check_noisy_bg(targetVcfPath, min_alt, min_altcell)
            sqlinline = f"{chr}\t{pos}\t{mutCount}\n"
            vcfOut.writelines(sqlinline)
            vcfOut.flush()

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outpathtmp = args.outdir
    id = args.id
    outpath = outpathtmp + "/" + id
    vcfpath = args.vcfpath
    binsize = args.binsize
    dpsize = args.dpsize
    altsize = args.altdpsize
    cellsize = args.cellsize
    altcellsize = args.altcellsize
    bamPath = args.bam
    refPath = args.ref
    outdirTmp = outpathtmp + "/tmp/"
    min_bq = args.min_bq
    min_mq = args.min_mq
    nprocs = args.nprocs
    cbtag = args.cbtag
    COREUSE = nprocs
    finalOutPath = outpath + ".scmutrace.cluster.count"
    if not os.path.exists(outdirTmp):
        os.makedirs(outdirTmp)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    lines_per_file = int(get_Vcf_rownum(vcfpath) / COREUSE)
    split_file_by_lines(vcfpath, "clean_" + id + "_", lines_per_file, outpathtmp)
    tmpSplitFiles = [f for f in os.listdir(outpathtmp) if f.startswith("splitVcf_clean_" + id + "_")]
    if COREUSE > 1:
        print('Parent process %s.' % os.getpid())
        p = mp.get_context("fork").Pool(COREUSE)
        print("Your cores is " + str(COREUSE))
        for tmpSplitFile in tmpSplitFiles:
            p.apply_async(count_cluster_snp, args=(outpathtmp + "/" + tmpSplitFile, outpathtmp + "/" + tmpSplitFile + ".vcf", bamPath, refPath, outdirTmp, min_bq, min_mq, dpsize, altsize, cellsize, altcellsize, 5, cbtag, binsize, id,), error_callback=err_call_back)
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')
    else:
        for tmpSplitFile in tmpSplitFiles:
            count_cluster_snp(outpathtmp + "/" + tmpSplitFile, outpathtmp + "/" + tmpSplitFile + ".vcf", bamPath, refPath, outdirTmp, min_bq, min_mq, dpsize, altsize, cellsize, altcellsize, 5, cbtag, binsize, id)

    with open(finalOutPath, "w") as finalOut:
        for tmpSplitFile in tmpSplitFiles:
            tmpVcfPath = outpathtmp + "/" + tmpSplitFile + ".vcf"
            with open(tmpVcfPath, "r") as tmpVcf:
                for line in tmpVcf:
                    finalOut.writelines(line)
    for tmpSplitFile in tmpSplitFiles:
        tmpVcfPath = outpathtmp + "/" + tmpSplitFile + ".vcf"
        try:
            subprocess.run(['rm', outpathtmp + "/" + tmpSplitFile], check=True)
            subprocess.run(['rm', tmpVcfPath], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error removing temporary files: {e}")
    subprocess.run(f'rm -rf {outdirTmp}/*', shell=True, check=True)

if __name__ == "__main__":
    start_time = timeit.default_timer()
    main()
    end_time = timeit.default_timer()
    print(f"Total running time: {end_time - start_time} seconds")