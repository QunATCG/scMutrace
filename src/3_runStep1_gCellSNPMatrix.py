import os
import timeit
import argparse
import pandas as pd
import numpy as np
from py_utils import *
import multiprocessing as mp

# Function to generate gCellSNP matrix
def initialize_parser():
    parser = argparse.ArgumentParser(description='generate gCellSNP matrix')
    parser.add_argument('--outdir', type=str, default='', help='out file directory', required = True)
    parser.add_argument('--jobconfig', type=str, default='', help='config file path, tab separate file, with the format [ID\tVCFPATH\tBAMPATH]', required = True)
    parser.add_argument('--reference', type=str, default='', help='reference path', required = True)
    parser.add_argument('--min_bq', type=int, default=0, help='cutoff base quality', required = False)
    parser.add_argument('--min_mq', type=int, default=255, help='cutoff mapping quality', required = False)
    parser.add_argument('--cbtag', type=str, default='CB', help='cell barcode tag', required = False)
    parser.add_argument('--nprocs', type=int, default=1, help='number of cores', required = False)
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')


# merge multiple mut info and remove duplicates
def merge_mut_info(jobconfigPath, mutinfoPath):
    tmpCheck = []
    with open(jobconfigPath, "r") as jobconfig, open(mutinfoPath, "w") as mutinfo:
        for line in jobconfig:
            lineSets = line.strip().split("\t")
            with open(lineSets[1], "r") as vcfTmp:
                for vcfLine in vcfTmp:
                    vcfLineSets = vcfLine.strip().split("\t")
                    if vcfLineSets[0] == "chr":
                        continue
                    if vcfLineSets[0] + "_" + vcfLineSets[1] not in tmpCheck:
                        tmpCheck.append(vcfLineSets[0] + "_" + vcfLineSets[1])
                        sql = str(vcfLineSets[0]) + "\t" + str(vcfLineSets[1]) + "\t" + str(vcfLineSets[2])
                        mutinfo.write(sql + "\n")

# check if a mutation is in a cell
def check_ifMut_oneCell(cellInfo, refBase):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    countDic = {base: 0 for base in basicBase}
    check_base = list(set(basicBase) - {refBase})
    cellInfoList = list(map(int, cellInfo.split("|")))
    for base, count in zip(basicBase, cellInfoList):
        countDic[base] += count
    numCheck = sum(countDic[base] for base in check_base)
    if numCheck >= 1:
        return "Mut"
    else:
        return "WT"

# get cell info of all mutations
def get_all_mut_info(mutinfoPath, jobconfigPath, outpathDir ,referencePath, min_bq, min_mq, cbtag):
    foutFile = open(mutinfoPath + "gCellSNPMatrix.step0.txt", "w")
    outdirTmp = outpathDir + "/tmp/"
    make_dir_check(outdirTmp)
    with open(mutinfoPath, "r") as mutinfo:
        for line in mutinfo:
            lineSets = line.strip().split("\t")
            chr = str(lineSets[0])
            pos = int(lineSets[1])
            refBase = str(lineSets[2])
            with open(jobconfigPath, "r") as jobconfig:
                for jobLine in jobconfig:
                    jobLineSets = jobLine.strip().split("\t")
                    idtmp = jobLineSets[0]
                    bamTmp = jobLineSets[2]
                    target_extractsnp_fromBAM(bamTmp, referencePath, outdirTmp, [chr, pos, pos], min_bq, min_mq, True, True, True, True, 1, 1, 1, 5, cbtag, 2, idtmp)
                    targetVcfPath = os.path.join(outdirTmp, f"VCF_{chr}_{pos}_{pos}_{idtmp}_random.vcf")
                    with open(targetVcfPath, "r") as targetVcf:
                        for vcfLine in targetVcf:
                            vcfLineSets = vcfLine.strip().split("\t")
                            posTmp = int(vcfLineSets[1])
                            if posTmp == pos:
                                cellsList = vcfLineSets[9].split(";")
                                cellsMutInfo = [check_ifMut_oneCell(i, refBase) for i in vcfLineSets[10].split(";")]
                                for cell, cellMutInfo in zip(cellsList, cellsMutInfo):
                                    foutFile.writelines(chr + "_" + str(pos) + "\t" + idtmp + "_" + cell + "\t" + str(cellMutInfo) + "\n")
                    # remove target vcf file
                    os.remove(targetVcfPath)
    foutFile.close()

# exclude control mutations
def exclude_control_mut(gCSMPath, jobconfigPath, outpathDir):
    foutFile = open(outpathDir + "/gCellSNPMatrix.step2.txt", "w")
    ControlIDs = []
    with open(jobconfigPath, "r") as jobconfig:
        for jobLine in jobconfig:
            jobLineSets = jobLine.strip().split("\t")
            if jobLineSets[3].strip() == "Control":
                ControlIDs.append(jobLineSets[0])
    SNPsControl = []
    with open(gCSMPath, "r") as gCSM:
        for line in gCSM:
            lineSets = line.strip().split("\t")
            if lineSets[1].split("_")[0] in ControlIDs and lineSets[2] == "Mut":
                SNPsControl.append(lineSets[0])
    with open(gCSMPath, "r") as gCSM:
        for line in gCSM:
            lineSets = line.strip().split("\t")
            if lineSets[0] not in SNPsControl:
                foutFile.write(line)
    foutFile.close()

# generate snp * sample mutation info
def generate_snp_sample_mutInfo(gCSMPath, outpath):
    tmpCheck = []
    with open(gCSMPath, "r") as gCSM, open(outpath, "w") as out:
        for line in gCSM:
            lineSets = line.strip().split("\t")
            pos = lineSets[0]
            idtmp = lineSets[1].split("_")[0]
            statustmp = lineSets[2]
            if pos + "_" + idtmp + "_" + statustmp not in tmpCheck:
                tmpCheck.append(pos + "_" + idtmp + "_" + statustmp)
                out.write(pos + "\t" + idtmp + "\t" + statustmp + "\n")

# generate gCellSNP matrix
def generate_gCellSNP_matrix(outpath, vcfpath):
    df = pd.read_csv(vcfpath, sep='\t', header=None, names=['Location', 'Cell', 'Status'])
    matrix = df.pivot_table(index='Location', columns='Cell', values='Status', aggfunc=lambda x: 'Mut' if 'Mut' in x.values else 'WT').fillna('NoCall')
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    matrix.to_csv(outpath, sep=',', index=True)

# generate gSNPSample matrix
def generate_gSNPSample_matrix(outpath, vcfpath):
    df = pd.read_csv(vcfpath, sep='\t', header=None, names=['Location', 'ID', 'Status'])
    matrix = df.pivot_table(index='Location', columns='ID', values='Status', aggfunc=lambda x: 'Mut' if 'Mut' in x.values else 'WT').fillna('NoCall')
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    matrix.to_csv(outpath, sep=',', index=True)

# Main function
def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outdir = args.outdir
    jobconfig = args.jobconfig
    min_bq = args.min_bq
    min_mq = args.min_mq
    reference = args.reference
    cbtag = args.cbtag
    COREUSE = args.nprocs
    make_dir_check(outdir)
    mutinfoPath = outdir + "/mutinfo.txt"
    print(">>> Merge all mutations info")
    merge_mut_info(jobconfig, mutinfoPath)
    totalMut = get_Vcf_rownum(mutinfoPath)
    split_file_by_lines(mutinfoPath, "split_mutInfo_", int(totalMut/COREUSE), outdir)
    split_fileList = [outdir + "/" + f for f in os.listdir(outdir) if "split_mutInfo_" in f]
    print(">>> Start to get all mutations info")
    if COREUSE > 1:
        print('Parent process %s.' % os.getpid())
        p = mp.get_context("fork").Pool(COREUSE)
        print("Your cores is " + str(COREUSE))
        for split_file in split_fileList:
            p.apply_async(get_all_mut_info, args=(split_file, jobconfig, outdir ,reference, min_bq, min_mq, cbtag,), error_callback=err_call_back)
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        print('All subprocesses done.')
    else:
        for split_file in split_fileList:
            get_all_mut_info(split_file, jobconfig, outdir ,reference, min_bq, min_mq, cbtag)
    # merge all mut info
    print(">>> Merge all mut info")
    for split_file in split_fileList:
        os.remove(split_file)
    merge_file = outdir + "/gCellSNPMatrix.step0.txt"
    split_gCellSM_fileList = [outdir + "/" + f for f in os.listdir(outdir) if "split_mutInfo_" in f]
    for split_file in split_gCellSM_fileList:
        with open(split_file, "r") as f:
            for line in f:
                with open(merge_file, "a") as merge_f:
                    merge_f.write(line)
        os.remove(split_file)
    # exclude control mutations
    print(">>> Exclude control mutations")
    exclude_control_mut(merge_file, jobconfig, outdir)
    os.remove(merge_file)
    print(">>> Generate gCellSNP matrix")
    finalResult = outdir + "/gCellSNPMatrix.step3.csv"
    generate_gCellSNP_matrix(finalResult, outdir + "/gCellSNPMatrix.step2.txt")
    print(">>> Generate gSNPSample matrix")
    outpath = outdir + "/gSNPSampleMatrix.step2.txt"
    generate_snp_sample_mutInfo(outdir + "/gCellSNPMatrix.step2.txt", outpath)
    finalResult = outdir + "/gSNPSampleMatrix.step3.csv"
    generate_gSNPSample_matrix(finalResult, outpath)


if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    end = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(end - start)) + ' seconds')