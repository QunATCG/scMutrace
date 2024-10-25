import os
import re
import timeit
import argparse
import subprocess

# get argparse
def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to obtain a list of base and cell counts in scRNA bam file')
    parser.add_argument('--vcfilelist', type=str, default=1, help='path file of VCFPOS folders in different samples [one path per line in path file]', required = True)
    parser.add_argument('--mincut', type=int, default=2, help='min number to get PoNs locs', required = False)
    parser.add_argument('--outdir', type=str, default=1, help='output dir', required = True)
    return (parser)

# get loc tmp file
def get_loc_tmp(pathIn, pathOut):
    with open(pathIn, "r") as vcfileIn, open(pathOut, "w") as vcfileOut:
        for line in vcfileIn:
            lineSets = line.split("\t")
            sql = lineSets[0] + "_" + lineSets[1] + "_" + lineSets[2]
            vcfileOut.writelines(sql + "\n")

# merge file
def merge_vcf_fromlist(vcfilelistPath, locMergepath):
    with open(vcfilelistPath, "r") as vcfilelist:
        for line in vcfilelist:
            line = line.strip()
            print(line)
            filesTotal1 = os.listdir(line)
            loctmpTotal = [i for i in filesTotal1 if re.search("loc.tmp", i) is not None]
            loctmpTotalPaths = [os.path.join(line, i) for i in loctmpTotal]
            print(loctmpTotal)
            if loctmpTotal:
                print("remove tmp files")
                loctmpTotalPathsjoin = " ".join(loctmpTotalPaths)
                cmdclean = ["rm"] + loctmpTotalPathsjoin
                print(" ".join(cmdclean))
                try:
                    subprocess.run(cmdclean, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error occurred while removing tmp files: {e}")
            filesTotal = os.listdir(line)
            VCFtmpTotal = [i for i in filesTotal if re.search("VCF_", i) is not None]
            VCFtmpTotalPaths = [os.path.join(line, i) for i in VCFtmpTotal]
            for vcfpath in VCFtmpTotalPaths:
                print(vcfpath)
                pathOut = f"{vcfpath}.loc.tmp"
                get_loc_tmp(vcfpath, pathOut)  
    loctmpPaths = []
    with open(vcfilelistPath, "r") as vcfilelist:
        for line in vcfilelist:
            line = line.strip()
            print(line)
            filesTotal = os.listdir(line)
            loctmpTotal = [i for i in filesTotal if re.search("loc.tmp", i) is not None]
            loctmpTotalPaths = [os.path.join(line, i) for i in loctmpTotal]
            loctmpPaths.extend(loctmpTotalPaths)
    jointPath = " ".join(loctmpPaths)
    cmdmerge = f"cat {jointPath} > {locMergepath}"
    try:
        subprocess.run(cmdmerge, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while merging files: {e}")
    cmdrmtmp = ["rm"] + loctmpPaths
    try:
        subprocess.run(cmdrmtmp, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while removing tmp files: {e}")

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    vcfilelistpath = args.vcfilelist
    cutoff = args.mincut
    outdir = args.outdir
    locMergepath = os.path.join(outdir, "PoN.merge")
    merge_vcf_fromlist(vcfilelistpath, locMergepath)
    
    # count number duplicated
    countPath = os.path.join(outdir, "PoN.merge.count")
    countPathfinal = os.path.join(outdir, "PoN.txt")
    cmdcount = f"sort {locMergepath} | uniq -c | awk '{{if($1>={cutoff}) print $0}}' | awk -v OFS='\t' '{{print $2}}' > {countPath}"
    try:
        subprocess.run(cmdcount, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while counting duplicates: {e}")
    
    with open(countPathfinal, "w") as outfile, open(countPath, "r") as ff:
        for line in ff:
            line = line.strip()
            if line:
                lineSets = line.split("_")
                sql = f"{lineSets[0]}\t{lineSets[1]}\t{lineSets[1]}"
                outfile.writelines(sql + "\n")
    
    cmdrm = f"rm {locMergepath} {countPath}"
    try:
        subprocess.run(cmdrm, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while removing temporary files: {e}")

if __name__=='__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')
