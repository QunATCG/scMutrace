import os
import re
import fnmatch
import timeit
import argparse
import pandas as pd
import numpy as np
from py_utils import *
import multiprocessing as mp
from collections import defaultdict

# Function to generate gSNPSample matrix
def initialize_parser():
    parser = argparse.ArgumentParser(description='generate gSNPSample matrix')
    parser.add_argument('--outdir', type=str, default='', help='out file directory', required = True)
    parser.add_argument('--jobconfig', type=str, default='', help='config file path, tab separate file, with the format [ID\tVCFPATH\tBAMPATH\tTAG]', required = True)
    parser.add_argument('--reference', type=str, default='', help='reference path', required = True)
    parser.add_argument('--min_bq', type=int, default=20, help='cutoff base quality', required = False)
    parser.add_argument('--min_mq', type=int, default=255, help='cutoff mapping quality', required = False)
    parser.add_argument('--cbtag', type=str, default='CB', help='cell barcode tag', required = False)
    parser.add_argument('--nprocs', type=int, default=1, help='number of cores', required = False)
    return (parser)

# define basic base
BASIC_BASE = ["A", "T", "C", "G", "D", "I"]

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

def merge_mut_info(jobconfigPath, mutinfoPath):
    # Use a set for O(1) average time complexity for lookups
    seen_mutations = set()
    pattern = r'([ATCGDI])\[(\d+)\]'
    try:
        with open(jobconfigPath, "r") as jobconfig, open(mutinfoPath, "w") as mutinfo:
            for line in jobconfig:
                lineSets = line.strip().split("\t")
                vcf_path = lineSets[1]
                try:
                    with open(vcf_path, "r") as vcfTmp:
                        for vcfLine in vcfTmp:
                            vcfLineSets = vcfLine.strip().split("\t")
                            # Skip header
                            if vcfLineSets[0] == "chr":
                                continue
                            mutation_key = f"{vcfLineSets[0]}_{vcfLineSets[1]}"
                            if mutation_key not in seen_mutations:
                                seen_mutations.add(mutation_key)
                                sql = f"{vcfLineSets[0]}\t{vcfLineSets[1]}\t{vcfLineSets[2]}"
                                mutBases = [f"{sql}\t{match[0]}" for match in re.findall(pattern, vcfLineSets[3])]
                                for mutBase in mutBases:
                                    mutinfo.write(mutBase + "\n")
                except IOError as e:
                    print(f"Unable to open VCF file: {vcf_path}, error: {e}")
    except IOError as e:
        print(f"Error opening job config file: {jobconfigPath}, error: {e}")

# mini check noisy bg
def minicheck_noisy_bg(targetVcfPath):
    data = open(targetVcfPath, "r")
    mutCount = 0
    for lineeach in data:
        line = lineeach.strip().split("\t")
        if int(line[5]) > 0:
            mutCount = mutCount + 1
    data.close()
    return mutCount

# check if a mutation is in a cell
def check_ifMut_oneCell(cellInfo, refBase):
    #basicBase = ["A", "T", "C", "G", "D", "I"]
    # Validate refBase
    if refBase not in BASIC_BASE:
        raise ValueError(f"Invalid reference base: {refBase}. Should be one of {BASIC_BASE}")
    # Initialize counts for each base
    countDic = {base: 0 for base in BASIC_BASE}
    # Create a list of checked bases
    check_base = list(set(BASIC_BASE) - {refBase})
    # Parse cell info into counts
    cellInfoList = list(map(int, cellInfo.split("|")))
    # Accumulate counts into the dictionary
    for base, count in zip(BASIC_BASE, cellInfoList):
        countDic[base] += count
    # Count the number of relevant mutations
    numCheck = sum(countDic[base] for base in check_base)
    # Determine mutation types
    mutype_list = [base for base in check_base if countDic[base] >= 1]
    mutype = "_".join(mutype_list)
    # Determine tag based on counts
    tag = "Mut" if numCheck >= 1 else "WT"
    return tag, mutype

# get cell info of all mutations
def get_all_mut_info(mutinfoPath, jobconfigPath, outpathDir, referencePath, min_bq, min_mq, cbtag):
    try:
        foutFile = open(mutinfoPath + "gCellSNPMatrix.step0.txt", "w")
        #foutFile_miniNoisy = open(mutinfoPath + "gCellSNPMatrix.step0_noisy.txt", "w")
    except IOError as e:
        print(f"Error opening output file: {e}")
        return
    outdirTmp = os.path.join(outpathDir, "tmp")
    make_dir_check(outdirTmp)
    try:
        with open(mutinfoPath, "r") as mutinfo:
            for line in mutinfo:
                lineSets = line.strip().split("\t")
                if len(lineSets) < 3:
                    print(f"Skipping invalid line in mutation info: {line}")
                    continue
                chr = str(lineSets[0])
                pos = int(lineSets[1])
                refBase = str(lineSets[2])
                with open(jobconfigPath, "r") as jobconfig:
                    for jobLine in jobconfig:
                        jobLineSets = jobLine.strip().split("\t")
                        # Validate jobLineSets
                        if len(jobLineSets) < 3:
                            print(f"Skipping invalid job line: {jobLine}")
                            continue
                        idtmp = jobLineSets[0]
                        bamTmp = jobLineSets[2]
                        target_extractsnp_fromBAM(bamTmp, referencePath, outdirTmp, [chr, pos, pos], min_bq, min_mq, True, True, True, True, 1, 0, 1, 5, cbtag, 5, idtmp)
                        #target_extractsnp_fromBAM(bamTmp, referencePath, outdirTmp, [chr, pos, pos], 0, min_mq, True, True, True, True, 1, 0, 1, 5, cbtag, 5, idtmp + "_noisy")
                        targetVcfPath = os.path.join(outdirTmp, f"VCF_{chr}_{pos}_{pos}_{idtmp}_random.vcf")
                        #targetVcfPath_noisy = os.path.join(outdirTmp, f"VCF_{chr}_{pos}_{pos}_{idtmp}_noisy_random.vcf")
                        #mutCount = minicheck_noisy_bg(targetVcfPath_noisy)
                        #f mutCount >= 2:
                        #    bgTag = "Noisy"
                        #else:
                        #    bgTag = "Clean"
                        #foutFile_miniNoisy.writelines(f"{chr}_{pos}_{refBase}\t{idtmp}\t{bgTag}\n")
                        #foutFile_miniNoisy.flush()
                        try:
                            with open(targetVcfPath, "r") as targetVcf:
                                for vcfLine in targetVcf:
                                    vcfLineSets = vcfLine.strip().split("\t")
                                    posTmp = int(vcfLineSets[1])
                                    # considering sequencing errors, 0.001 is the cutoff
                                    if int(vcfLineSets[5])/int(vcfLineSets[3]) < 0.001:
                                        continue
                                    if posTmp == pos:
                                        cellsList = vcfLineSets[9].split(";")
                                        cellsMutInfo = [check_ifMut_oneCell(i, refBase) for i in vcfLineSets[10].split(";")]
                                        for cell, cellMutInfo in zip(cellsList, cellsMutInfo):
                                            mutypes = cellMutInfo[1].split("_")
                                            for m in mutypes:
                                                if str(cellMutInfo[0]) == "Mut":
                                                    foutFile.writelines(f"{chr}_{pos}_{refBase}_{m}\t{idtmp}_{cell}\t{cellMutInfo[0]}\n")
                                                    foutFile.flush()
                                                if str(cellMutInfo[0]) == "WT":
                                                    foutFile.writelines(f"{chr}_{pos}_{refBase}_{refBase}\t{idtmp}_{cell}\t{cellMutInfo[0]}\n")
                                                    foutFile.flush()
                        except IOError as e:
                            print(f"Error reading VCF file {targetVcfPath}: {e}")
                        finally:
                            # Cleanup the VCF file
                            if os.path.exists(targetVcfPath):
                                os.remove(targetVcfPath)
                            #if os.path.exists(targetVcfPath_noisy):
                            #    os.remove(targetVcfPath_noisy)
    except IOError as e:
        print(f"Error reading mutation info file: {e}")
    finally:
        foutFile.close()

# focus on the mutation info
def filter_mut_info(inputpath, outputpath, mutinfoPath):
    mutinfo = set()
    with open(mutinfoPath, "r") as mutinfoFile:
        for line in mutinfoFile:
            chrom, pos, ref_base, mut_base = line.strip().split("\t")
            tag_mut = f"{chrom}_{pos}_{ref_base}_{mut_base}"
            tag_ref = f"{chrom}_{pos}_{ref_base}_{ref_base}"
            mutinfo.add(tag_mut)
            mutinfo.add(tag_ref)
    with open(inputpath, "r") as inputFile, open(outputpath, "w") as outputFile:
        for line in inputFile:
            fields = line.strip().split("\t")
            chrom, pos, ref_base, mut_base = fields[0].split("_")
            tag = f"{chrom}_{pos}_{ref_base}_{mut_base}"
            if tag in mutinfo:
                outputFile.write(line)

def modify_mutinfo(inputpath, outputpath):
    mutDic = defaultdict(list)
    print(">>> Processing input file")
    try:
        with open(inputpath, "r") as inputFile:
            for line in inputFile:
                fields = line.strip().split("\t")
                if len(fields) < 2:
                    continue
                try:
                    chrom, pos, ref_base, mut_base = fields[0].split("_")
                except ValueError:
                    print(f"Skipping line due to format issue: {line.strip()}")
                    continue
                loc = f"{chrom}_{pos}"
                if mut_base != ref_base:
                    if mut_base not in mutDic[loc]:
                        mutDic[loc].append(mut_base)
    except FileNotFoundError:
        print(f"Input file not found: {inputpath}")
        return
    print(">>> Generating output")
    try:
        with open(outputpath, "w") as outputFile:
            with open(inputpath, "r") as inputFile:
                for line in inputFile:
                    fields = line.strip().split("\t")
                    if len(fields) < 2:
                        continue
                    try:
                        chrom, pos, ref_base, mut_base = fields[0].split("_")
                    except ValueError:
                        print(f"Skipping line due to format issue: {line.strip()}")
                        continue
                    loc = f"{chrom}_{pos}"
                    for mut in mutDic[loc]:
                        output_tag = f"{loc}_{ref_base}_{mut}"
                        if mut_base == mut:
                            status = "Mut"
                        else:
                            status = "WT"
                        outputFile.write(f"{output_tag}\t{fields[1]}\t{status}\n")
    except Exception as e:
        print(f"An error occurred while writing to the output file: {e}")

# exclude control mutations
def exclude_control_mut(gCSMPath, jobconfigPath, outpathDir):
    foutFile = open(outpathDir + "/gCellSNPMatrix.step2.txt", "w")
    ControlIDs = []
    with open(jobconfigPath, "r") as jobconfig:
        for jobLine in jobconfig:
            jobLineSets = jobLine.strip().split("\t")
            if jobLineSets[3].strip() == "Control":
                ControlIDs.append(jobLineSets[0])
    SNPsControl = set()
    with open(gCSMPath, "r") as gCSM:
        for line in gCSM:
            lineSets = line.strip().split("\t")
            if lineSets[1].split("_")[0] in ControlIDs and lineSets[2] == "Mut":
                SNPsControl.add(lineSets[0])
    with open(gCSMPath, "r") as gCSM:
        for line in gCSM:
            lineSets = line.strip().split("\t")
            if lineSets[0] not in SNPsControl:
                foutFile.writelines(line)
                foutFile.flush()
    foutFile.close()

# generate snp * sample mutation info
def generate_snp_sample_mutInfo(gCSMPath, outpath):
    seen_entries = set()
    with open(gCSMPath, "r") as gCSM, open(outpath, "w") as out:
        for line in gCSM:
            lineSets = line.strip().split("\t")
            pos = lineSets[0]
            idtmp = "_".join(lineSets[1].split("_")[0:-1])
            statustmp = lineSets[2]
            entry = f"{pos}_{idtmp}_{statustmp}"
            if entry not in seen_entries:
                seen_entries.add(entry)
                out.write(f"{pos}\t{idtmp}\t{statustmp}\n")

# generate gCellSNP matrix
def generate_gCellSNP_matrix(outpath, vcfpath):
    df = pd.read_csv(vcfpath, sep='\t', header=None, names=['Location', 'Cell', 'Status'])
    matrix = df.pivot_table(index='Location', columns='Cell', values='Status', aggfunc=lambda x: 'Mut' if 'Mut' in x.values else 'WT').fillna('NoCall')
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    matrix.to_csv(outpath, sep=',', index=True)

def generate_gCellSNP_matrix_chunked_sparse(outpath, vcfpath, chunksize=100000):
    """
    分块读取VCF-like文件，逐步构建稀疏突变矩阵，避免内存溢出。

    参数：
    -------
    outpath : str
        输出csv路径
    vcfpath : str
        输入vcf-like文件路径（tab分隔，三列：Location, Cell, Status）
    chunksize : int
        每次读取行数，默认 1,000,000
    """
    # 用嵌套dict存稀疏矩阵，外层是Location，内层是Cell: Status
    matrix = defaultdict(dict)
    all_cells = set()
    # 分块读取
    for chunk in pd.read_csv(vcfpath, sep='\t', header=None, names=['Location', 'Cell', 'Status'], chunksize=chunksize):
        for row in chunk.itertuples(index=False):
            loc, cell, status = row
            all_cells.add(cell)
            # 如果当前已有'Mut'，跳过，优先保留'Mut'
            if matrix[loc].get(cell) == 'Mut':
                continue
            matrix[loc][cell] = status
    # 转DataFrame（稀疏存储）
    all_cells = sorted(all_cells)
    df = pd.DataFrame(index=matrix.keys(), columns=all_cells, dtype="string")
    for loc, cell_statuses in matrix.items():
        for cell, status in cell_statuses.items():
            df.at[loc, cell] = status
    # 填充缺失值为'NoCall'
    df = df.fillna('NoCall')
    # 转稀疏类型，节省内存
    for col in df.columns:
        df[col] = pd.Series(pd.arrays.SparseArray(df[col], fill_value='NoCall'))
    # 输出
    df.to_csv(outpath, sep=',', index=True)
    print(f"Done! Matrix shape: {df.shape}")

# generate gSNPSample matrix
def generate_gSNPSample_matrix(outpath, vcfpath):
    df = pd.read_csv(vcfpath, sep='\t', header=None, names=['Location', 'ID', 'Status'])
    matrix = df.pivot_table(index='Location', columns='ID', values='Status', aggfunc=lambda x: 'Mut' if 'Mut' in x.values else 'WT').fillna('NoCall')
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    matrix.to_csv(outpath, sep=',', index=True)

def generate_gSNPSample_matrix_chunked_sparse(outpath, vcfpath, chunksize=100000):
    """
    分块读取VCF-like文件，逐步构建稀疏突变矩阵（Location × Sample）。

    参数：
    -------
    outpath : str
        输出csv路径
    vcfpath : str
        输入vcf-like文件路径（tab分隔，三列：Location, ID, Status）
    chunksize : int
        每次读取行数，默认 1,000,000
    """
    # 用嵌套dict存稀疏矩阵，外层是Location，内层是ID: Status
    matrix = defaultdict(dict)
    all_ids = set()
    # 分块读取
    for chunk in pd.read_csv(vcfpath, sep='\t', header=None, names=['Location', 'ID', 'Status'], chunksize=chunksize):
        for row in chunk.itertuples(index=False):
            loc, sample_id, status = row
            all_ids.add(sample_id)
            # 如果当前已有'Mut'，跳过，优先保留'Mut'
            if matrix[loc].get(sample_id) == 'Mut':
                continue
            matrix[loc][sample_id] = status
    # 转DataFrame（稀疏存储）
    all_ids = sorted(all_ids)
    df = pd.DataFrame(index=matrix.keys(), columns=all_ids, dtype="string")
    for loc, id_statuses in matrix.items():
        for sample_id, status in id_statuses.items():
            df.at[loc, sample_id] = status
    # 填充缺失值为'NoCall'
    df = df.fillna('NoCall')
    # 转稀疏类型，节省内存
    for col in df.columns:
        df[col] = pd.Series(pd.arrays.SparseArray(df[col], fill_value='NoCall'))
    # 输出
    df.to_csv(outpath, sep=',', index=True)
    print(f"Done! Matrix shape: {df.shape}")

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
    split_gCellSM_fileList = [outdir + "/" + f for f in os.listdir(outdir) if fnmatch.fnmatch(f, "splitVcf_split_mutInfo_*gCellSNPMatrix.step0.txt")]
    os.system(f"cat {' '.join(split_gCellSM_fileList)} >> {merge_file}")
    for split_file in split_gCellSM_fileList:
        os.remove(split_file)
    # filter mut info
    print(">>> Filter mut info")
    filter_mut_info(merge_file, outdir + "/gCellSNPMatrix.step0_0.txt", mutinfoPath)
    # modify mut info
    print(">>> Modify mut info")
    modifyedFile = outdir + "/gCellSNPMatrix.step1.txt"
    modify_mutinfo(outdir + "/gCellSNPMatrix.step0_0.txt", modifyedFile)
    # exclude control mutations
    print(">>> Exclude control mutations")
    exclude_control_mut(modifyedFile, jobconfig, outdir)
    print(">>> Generate gSNPSample matrix")
    outpath1 = outdir + "/gSNPSampleMatrix.step2.1.txt"
    outpath2 = outdir + "/gSNPSampleMatrix.step2.2.txt"
    generate_snp_sample_mutInfo(outdir + "/gCellSNPMatrix.step2.txt", outpath1)
    generate_snp_sample_mutInfo(outdir + "/gCellSNPMatrix.step1.txt", outpath2)
    finalResult3 = outdir + "/gSNPSampleMatrix.step3.1_excludeControlSites.csv"
    finalResult4 = outdir + "/gSNPSampleMatrix.step3.2_includeControlSites.csv"
    generate_gSNPSample_matrix(finalResult3, outpath1)
    generate_gSNPSample_matrix(finalResult4, outpath2)


if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    end = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(end - start)) + ' seconds')