import os
import re
import timeit
import argparse
import multiprocessing as mp
from py_utils import *
import subprocess

def initialize_parser():
    parser = argparse.ArgumentParser(description='filter SNP based on 3Model')
    parser.add_argument('--outdir', type=str, help='outdir', required = True)
    parser.add_argument('--vcfdir', type=str, help='VcfPath', required = True)
    parser.add_argument('--id', help='Prefix used to name output file.', required = True)
    parser.add_argument('--a1', type=float, default=0.0056773121, help='alphaAltCountEst [Default: 0.0056773121]', required = False)
    parser.add_argument('--b1', type=float, default=1.3777933664, help='betaRefCountEst [Default: 1.3777933664]', required = False)
    parser.add_argument('--a2', type=float, default=0.0098542673, help='alphaMutCellCountEst [Default: 0.0098542673]', required = False)
    parser.add_argument('--b2', type=float, default=1.6730916145, help='betaRefCellCountEst [Default: 1.6730916145]', required = False)
    parser.add_argument('--maxcheck', action='store_true', help='A boolean flag for best performance checking, it will use full cpus')
    parser.add_argument('--model1', action='store_true', help='A boolean flag for model1, it will use samtools to identify snps')
    parser.add_argument('--model2', action='store_true', help='A boolean flag for model2, it will use beta distribution statistics to identify snps')
    parser.add_argument('--model3', action='store_true', help='A boolean flag for model3, it will use genotype mode to identify snps')
    parser.add_argument('--nprocs', default = 1, help='Number of processes [Default: 1]',required=False,type = int)
    parser.add_argument('--iternum', default = 1000, help='Number of iteration [Default: 1000]',required=False,type = int)
    parser.add_argument('--errorp', type=float, default=0.006, help='estimate sequecing error [Default: 0.006]', required = False)
    parser.add_argument('--backp', type=float, default=0.05, help='background mutation rate [Default: 0.03-0.05]', required = False)
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

def merge_filteredVcf_file(sampleID, filteredVcfPaths, dirPath):
    """
    Merges filtered VCF files for a given sample ID and chromosome, and removes the original files.

    Args:
        sampleID (str): The sample ID used to identify the VCF files.
        filteredVcfPaths (list of str): List of paths to the filtered VCF files.
        dirPath (str): Directory path where the merged VCF files will be saved.

    Raises:
        subprocess.CalledProcessError: If an error occurs while executing the shell commands.
    """
    chrs = list(set([i.split(sampleID + "_")[1].split(".")[0] for i in filteredVcfPaths]))
    for chr in chrs:
        tmpipaths = [j for j in filteredVcfPaths if re.search(sampleID + "_" + str(chr) + ".vcf", j) != None]
        tmpsql = " ".join(tmpipaths)
        cmdmerge = f"cat {tmpsql} > {dirPath}{sampleID}_{chr}_filter.vcf"
        cmdrm = f"rm {tmpsql}"
        try:
            subprocess.run(cmdmerge, shell=True, check=True)
            subprocess.run(cmdrm, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while executing command: {e.cmd}")
            print(f"Return code: {e.returncode}")
            print(f"Output: {e.output}")


def filter_snp_threeModel(VCFPath, outPath, a1, b1, a2, b2, model3ProbabilityDic, quainfos, iternum, errorP, backP, m1check, m2check, m3check):
    """
    Filters SNPs using a three-model approach and writes the results to an output file.

    Args:
        VCFPath (str): Path to the input VCF file.
        outPath (str): Path to the output file where filtered results will be written.
        a1 (float): Parameter a1 for the genotype assignment model.
        b1 (float): Parameter b1 for the genotype assignment model.
        a2 (float): Parameter a2 for the genotype assignment model.
        b2 (float): Parameter b2 for the genotype assignment model.
        model3ProbabilityDic (dict): Dictionary containing probabilities for the third model.
        quainfos (list): List of quality information.
        iternum (int): Number of iterations for the genotype assignment.
        errorP (float): Error probability.
        backP (float): Background probability.
        m1check (bool): Check flag for model 1.
        m2check (bool): Check flag for model 2.
        m3check (bool): Check flag for model 3.

    Raises:
        IndexError: If there is an indexing error while processing the VCF file.
        ValueError: If there is a value error while processing the VCF file.

    Returns:
        None
    """
    try:
        with open(VCFPath, "r") as vcfile, open(outPath, "w+") as outfile:
            fileRowNum = get_Vcf_rownum(VCFPath)
            counttmp = 0
            for line in vcfile:
                counttmp += 1
                if counttmp % 5000 == 0:
                    print("Processed " + str(counttmp) + "/" + str(fileRowNum) + " positions")
                lineSets = line.split("\t")
                refBase = lineSets[2]
                baseQualityInfoOnesite = lineSets[11]
                dptotal = int(lineSets[3])
                if dptotal <= 10:
                    testResult = assign_site_genotype(refBase, baseQualityInfoOnesite, a1, b1, a2, b2, False, model3ProbabilityDic, quainfos, 1000, iternum, errorP, backP, m1check, m2check, m3check)
                else:
                    testResult = assign_site_genotype(refBase, baseQualityInfoOnesite, a1, b1, a2, b2, True, model3ProbabilityDic, quainfos, 1000, iternum, errorP, backP, m1check, m2check, m3check)
                filterinfo = ";".join([str(testResult["m1"]), str(testResult["m2"]), str(testResult["m3"])])
                statpvalue = ";".join([str(testResult["p1"]), str(testResult["p2"]), str(testResult["p3"])])
                iteresampleNum = ";".join([str(testResult["ic1"]), str(testResult["ic2"]), str(testResult["ic3"]), str(iternum)])
                filtercheck = "NonePASS"
                supportEvidence = "NoneSupport"
                if testResult["m1"] + testResult["m2"] + testResult["m3"] >= 1:
                    filtercheck = "PASS_3Model"
                    supportEvidence = "Strong"
                if filtercheck == "NonePASS" and (testResult["c1"] + testResult["c2"] + testResult["c3"]) >= 1:
                    supportEvidence = "Weak"
                    filtercheck = "PASS_3Model"
                sql = line.strip() + "\t" + filterinfo + "\t" + filtercheck + "\t" + statpvalue + "\t" + supportEvidence + "\t" + iteresampleNum
                outfile.writelines(sql + "\n")
                outfile.flush()
            print("Done for " + VCFPath)
            cmd = ["rm", VCFPath]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while executing command: {e.cmd}")
                print(f"Return code: {e.returncode}")
                print(f"Output: {e.output}")
    except (IndexError, ValueError) as e:
        print(f"An error occurred: {e}")

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outPath = args.outdir
    VCFDir = args.vcfdir
    IDSample = args.id
    errorP = args.errorp
    backP = args.backp
    a1 = args.a1
    b1 = args.b1
    a2 = args.a2
    b2 = args.b2
    m1check = args.model1
    m2check = args.model2
    m3check = args.model3
    iternum = args.iternum
    maxcheck = args.maxcheck
    COREUSE=args.nprocs
    # try to use max cores
    if maxcheck:
        if COREUSE != 1 and mp.cpu_count() > COREUSE:
            COREUSE = mp.cpu_count() - 1
            print("Try to use max cores number...")
        if COREUSE > mp.cpu_count():
            COREUSE = mp.cpu_count() - 1
            print("Cores is over computer max cpus...")
    else:
        if COREUSE > mp.cpu_count():
            COREUSE = mp.cpu_count() - 1
            print("Cores is over computer max cpus...")
    filterPath = outPath + "filterVCF/"
    model3ProbabilityDic = create_model3_probabilityDic(errorP)
    quainfos = list(range(10, 100, 1))
    make_dir_check(filterPath)
    # split file
    print("Splitting files...")
    split_multi_files(COREUSE, VCFDir)
    vcfilestotal_pos_vcf = os.listdir(VCFDir)
    vcfilestotal = [i for i in vcfilestotal_pos_vcf if re.match("splitVcf_", i) != None]
    vcfilesInpaths = [VCFDir + i for i in vcfilestotal]
    print(vcfilesInpaths)
    vcfilesOutpaths = [filterPath + IDSample + "_" + i.split("_")[-1] + "_filter.vcf" for i in vcfilestotal]
    print(vcfilesOutpaths)
    NumberVCF = len(vcfilesInpaths)
    if COREUSE > 1:
        print('Parent process %s.' % os.getpid())
        p = mp.get_context("fork").Pool(COREUSE)
        print("Your cores is " + str(COREUSE))
        results = []
        for i in range(0, NumberVCF):
            result = p.apply_async(filter_snp_threeModel, args = (vcfilesInpaths[i], vcfilesOutpaths[i], a1, b1, a2, b2, model3ProbabilityDic, quainfos, iternum, errorP, backP, m1check, m2check, m3check,), error_callback=err_call_back)
            results.append(result)
        print('Waiting for all subprocesses done...')
        p.close()
        p.join()
        for result in results:
            print(result.get())
        print('All subprocesses done.')
    else:
        for i in range(0, len(vcfilesInpaths)):
            filter_snp_threeModel(vcfilesInpaths[i], vcfilesOutpaths[i], a1, b1, a2, b2, model3ProbabilityDic, quainfos, iternum, errorP, backP, m1check, m2check, m3check)
    print("Merging vcf files...")
    merge_filteredVcf_file(IDSample, vcfilesOutpaths, filterPath)

if __name__=='__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')
