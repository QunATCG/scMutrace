import os
import re
import sys
import timeit
import argparse
from py_utils import *
import subprocess
import multiprocessing as mp

def initialize_parser():
    parser = argparse.ArgumentParser(description='filter SNP based on public database')
    parser.add_argument('--outdir', type=str, default="", help='outdir', required = True)
    parser.add_argument('--casevcfdir', type=str, default="", help='Vcfdir', required = True)
    parser.add_argument('--caseid', type=str, default="", help='Prefix used to name output file', required = True)
    parser.add_argument('--exclude', type=str, default="", help='files to be excluded, tab seperated format: path\tid', required = False)
    parser.add_argument('--include', type=str, default="", help='files to be included, tab seperated format: path\tid', required = False)
    parser.add_argument('--maxcheck', action='store_true', help='A boolean flag for best performance checking, it will use full cpus')
    parser.add_argument('--nprocs', default = 1, help='Number of processes [Default: 1]', required=False, type = int)
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

def check_files_dir(vcfdir):
    """
    Checks the specified directory for files with names containing "filter_".
    If such files are found, they are deleted.

    Args:
        vcfdir (str): The directory to check for files.

    Raises:
        subprocess.CalledProcessError: If an error occurs while attempting to delete the files.

    Example:
        check_files_dir("/path/to/vcfdir")
    """
    vcfiles = os.listdir(vcfdir)
    filterCheck = [i for i in vcfiles if re.search("filter_", i) != None]
    if len(filterCheck) >= 1:
        print("clean vcfiles...")
        filterCheckPaths = [os.path.join(vcfdir, i) for i in filterCheck]
        try:
            subprocess.run(["rm"] + filterCheckPaths, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while cleaning vcfiles: {e}")

def merge_out_vcfiles(vcfdir, idTag):
    """
    Merges VCF files in a specified directory that match a given ID tag and removes the original files.

    Args:
        vcfdir (str): The directory containing the VCF files to be merged.
        idTag (str): The identifier tag used to filter the VCF files.

    Raises:
        subprocess.CalledProcessError: If an error occurs during the execution of the shell commands.

    Example:
        merge_out_vcfiles('/path/to/vcfdir', 'sampleID')
    """
    vcfiles = os.listdir(vcfdir)
    filterCheckPaths = [os.path.join(vcfdir, i) for i in vcfiles if re.search("filter_" + str(idTag), i) is not None]
    filesCheck = " ".join(filterCheckPaths)
    if len(filterCheckPaths) < 1:
        print("No files to merge.")
        cmdempty = f"touch {os.path.join(vcfdir, str(idTag) + '_scmutrace.vcf')}"
        subprocess.run(cmdempty, shell=True, check=True)
        return
    sqlCheckmerge = f"cat {filesCheck} > {os.path.join(vcfdir, str(idTag) + '_scmutrace.vcf')}"
    print(sqlCheckmerge)
    sqlCheckrm = f"rm {filesCheck}"
    try:
        subprocess.run(sqlCheckmerge, shell=True, check=True)
        subprocess.run(sqlCheckrm, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")

def bedtools_intersect_v(path1, path2, pathdatabase):
    """
    Perform an intersection between two BED files using bedtools and return common and unique positions.

    This function uses the bedtools intersect command to find overlapping regions between the input BED file
    (path1) and the database BED file (pathdatabase). The results are written to an output file (path2).
    It then parses the output to identify common and unique positions.

    Args:
        path1 (str): Path to the input BED file.
        path2 (str): Path to the output file where the intersection results will be written.
        pathdatabase (str): Path to the database BED file.

    Returns:
        tuple: A tuple containing two lists:
            - locListCommon (list): List of positions that are common between the input and database BED files.
            - locListUniq (list): List of positions that are unique to the input BED file.

    Raises:
        FileNotFoundError: If any of the specified files (path1 or pathdatabase) do not exist.
        subprocess.CalledProcessError: If the bedtools intersect command fails.
    """
    if not os.path.exists(path1):
        raise FileNotFoundError(f"File not found: {path1}")
    if not os.path.exists(pathdatabase):
        raise FileNotFoundError(f"File not found: {pathdatabase}")
    # Check if bedtools is installed
    if subprocess.run(["which", "bedtools"], capture_output=True, text=True).returncode != 0:
        raise EnvironmentError("bedtools is not installed or not found in the system PATH.")
    cmd = ["bedtools", "intersect", "-a", path1, "-b", pathdatabase, "-wa", "-wb"]
    try:
        with open(path2, "w") as output_file:
            subprocess.run(cmd, stdout=output_file, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running bedtools intersect: {e}")
        return [], []
    locListCommon = []
    try:
        with open(path2, "r") as bedFile:
            for line in bedFile:
                lineSets = line.split("\t")
                sourcePos = lineSets[1]
                targetPosStart = lineSets[4]
                targetPosEnd = lineSets[5]
                if int(sourcePos) < int(targetPosStart) or int(sourcePos) > int(targetPosEnd):
                    continue
                if sourcePos not in locListCommon:
                    locListCommon.append(lineSets[1])
    except FileNotFoundError as e:
        print(f"Error occurred while reading {path2}: {e}")
        return [], []
    locListUniq = []
    try:
        with open(path1, "r") as bedFile2:
            for line in bedFile2:
                lineSets = line.split("\t")
                if lineSets[1] not in locListCommon and lineSets[1] not in locListUniq:
                    locListUniq.append(lineSets[1])
    except FileNotFoundError as e:
        print(f"Error occurred while reading {path1}: {e}")
        return locListCommon, []
    return locListCommon, locListUniq

def sort_vcf_loc(PathVCFIn, PathVCFOut):
    """
    Sorts a VCF file by chromosome and position.

    Args:
        PathVCFIn (str): Path to the input VCF file.
        PathVCFOut (str): Path to the output sorted VCF file.

    Returns:
        None
    """
    cmd = f"sort -k1,1 -k2,2n {PathVCFIn} > {PathVCFOut}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while sorting VCF file: {e}")

def filter_SNP_f0(VCFile, pathDir, pathOut):
    """
    Filters SNPs based on the scmutrace model by removing lines containing 'NonePASS'.

    Args:
        VCFile (str): The name of the VCF file to be filtered.
        pathDir (str): The directory where the VCF file is located.
        pathOut (str): The path to the output file where the filtered results will be written.

    Returns:
        None
    """
    input_path = os.path.join(pathDir, VCFile)
    cmd = f"grep -v NonePASS {input_path} > {pathOut}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while filtering SNPs: {e}")

def filter_SNP_f1(PathVCFIn, PathVCFOut, PathVCFInBed, excludeFilePath=None, includeFilePath=None):
    # filter snps based on public databases
    ex_dic = {}
    in_dic = {}
    # exclude
    if excludeFilePath:
        with open(excludeFilePath, "r") as excludeFile_total:
            for line in excludeFile_total:
                exPath = line.split("\t")[0].strip()
                exID = line.split("\t")[1].strip()
                print(exID)
                tmpOut_ex = PathVCFInBed + "." + str(exID)
                _, ex_locList = bedtools_intersect_v(PathVCFInBed, tmpOut_ex, exPath)
                ex_dic[exID] = ex_locList
                subprocess.run(["rm", tmpOut_ex], check=True)
    # include
    if includeFilePath:
        with open(includeFilePath, "r") as includeFile_total:
            for line in includeFile_total:
                inPath = line.split("\t")[0].strip()
                inID = line.split("\t")[1].strip()
                print(inID)
                tmpOut_in = PathVCFInBed + "." + str(inID)
                in_locList, _ = bedtools_intersect_v(PathVCFInBed, tmpOut_in, inPath)
                in_dic[inID] = in_locList
                subprocess.run(["rm", tmpOut_in], check=True)
    with open(PathVCFIn, "r") as VCFile, open(PathVCFOut, "w") as VCFout:
        for line in VCFile:
            tag = []
            loc = line.split("\t")[1]
            # exclude
            if not ex_dic:
                tag.append("EmptyEx")
            else:
                for tmp_id_ex, tmp_list_ex in ex_dic.items():
                    if loc in tmp_list_ex:
                        tag.append("PASS_" + tmp_id_ex)
                    else:
                        tag.append("NonePASS_" + tmp_id_ex)
            # include
            if not in_dic:
                tag.append("EmptyIn")
            else:
                for tmp_id_in, tmp_list_in in in_dic.items():
                    if loc in tmp_list_in:
                        tag.append("In_" + tmp_id_in)
                    else:
                        tag.append("Not_in_" + tmp_id_in)
            tagsql = ";".join(tag)
            sql = line.strip() + "\t" + tagsql + "\n"
            VCFout.writelines(sql)

def filter_multi_process(outDir, pathDir, VCFileName, excludeFilePath, includeFilePath):
    """
    Filters SNPs from a VCF file using a multi-step process and moves the final output to a specified directory.

    Args:
        outDir (str): The directory where the final filtered VCF file will be moved.
        pathDir (str): The directory where intermediate files will be stored.
        VCFileName (str): The name of the VCF file to be processed.
        excludeFilePath (str): The path to the file containing SNPs to be excluded.
        includeFilePath (str): The path to the file containing SNPs to be included.

    Returns:
        None
    """
    print("processing..", VCFileName)
    tmpsampleID = VCFileName.split("_filter.vcf")[0]
    Path_f0_out = os.path.join(pathDir, f"filter_{tmpsampleID}.f0.vcf")
    Path_f0_bed = os.path.join(pathDir, f"filter_{tmpsampleID}.f0.bed")
    Path_f1_out = os.path.join(pathDir, f"filter_{tmpsampleID}.f1.vcf")
    filter_SNP_f0(VCFileName, pathDir, Path_f0_out)
    with open(Path_f0_out, "r") as vcff0, open(Path_f0_bed, "w") as bedf0:
        for line in vcff0:
            linesets = line.split("\t")
            sql = f"{linesets[0].strip()}\t{linesets[1].strip()}\t{linesets[1].strip()}"
            bedf0.write(sql + "\n")
    filter_SNP_f1(Path_f0_out, Path_f1_out, Path_f0_bed, excludeFilePath, includeFilePath)
    try:
        subprocess.run(["mv", Path_f1_out, outDir], check=True)
        subprocess.run(["rm", Path_f0_out, Path_f0_bed], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while moving or removing files: {e}")

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outDir = args.outdir
    make_dir_check(outDir)
    VCFDir = args.casevcfdir
    caseid = args.caseid
    check_files_dir(VCFDir)
    COREUSE = args.nprocs
    excludeFilePath = args.exclude
    includeFilePath = args.include
    if excludeFilePath == 1 or excludeFilePath == "" or excludeFilePath == "NotAvailable":
        excludeFilePath = None
    if includeFilePath == 1 or includeFilePath == "" or includeFilePath == "NotAvailable":
        includeFilePath = None
    maxcheck = args.maxcheck
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
    vcfilestotal_raw = os.listdir(VCFDir)
    vcfilestotal = [i for i in vcfilestotal_raw if re.search("vcf", i) != None]
    if COREUSE > 1:
        print('Parent process %s.' % os.getpid())
        p = mp.get_context("fork").Pool(COREUSE)
        print("Your cores is " + str(COREUSE))
        for vcfile in vcfilestotal:
            p.apply_async(filter_multi_process, args = (outDir, VCFDir, vcfile, excludeFilePath, includeFilePath,), error_callback = err_call_back)
    else:
        for vcfile in vcfilestotal:
            filter_multi_process(outDir, VCFDir, vcfile, excludeFilePath, includeFilePath)
    print("Merging files...")
    merge_out_vcfiles(outDir, caseid)

if __name__=='__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')