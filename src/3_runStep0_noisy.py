import os
import timeit
import argparse
import subprocess
import scipy.stats as stats

def initialize_parser():
    parser = argparse.ArgumentParser(description='filter SNP based on public database')
    parser.add_argument('--outdir', type=str, default='', help='out file path', required = True)
    parser.add_argument('--vcfpath', type=str, default='', help='Vcf path', required = True)
    parser.add_argument('--id', type=str, default='', help='sample id', required = True)
    parser.add_argument('--barcode', type=str, default='', help='barcode file [tab separated file with header: Barcode CellType]', required = True)
    parser.add_argument('--binsize', type=int, default=5, help='base number between two snps (MultiAlleles)', required = False)
    parser.add_argument('--clustersize', type=int, default=4, help='number of MultiAlleles snps', required = False)
    parser.add_argument('--dpsize', type=int, default=10, help='min depth', required = False)
    parser.add_argument('--altdpsize', type=int, default=5, help='min alternative depth', required = False)
    parser.add_argument('--cellsize', type=int, default=10, help='min cell number', required = False)
    parser.add_argument('--altcellsize', type=int, default=5, help='min mut cell number', required = False)
    parser.add_argument('--lowBaseCutoff', type=int, default=20, help='cutoff low base quality', required = False)
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

# filter pos list based on bin size and cluster size
def filter_pos_binCluster(posList, binsize, clustersize):
    posListSort = sorted(set(int(i) for i in posList))
    excludeList = []
    posListmp = []
    for pos in posListSort:
        if not posListmp:
            posListmp.append(pos)
        elif pos - posListmp[-1] <= binsize:
            posListmp.append(pos)
        else:
            if len(posListmp) >= clustersize:
                excludeList.extend(posListmp)
            posListmp = [pos]
    if len(posListmp) >= clustersize:
        excludeList.extend(posListmp)
    return excludeList

# assign indel and snp
def assign_mutype_snpindex(sitemutinfo, refBase):
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
    gtListAllele = [base for base in check_base if countDic[base] > 1]
    gtListAlleleEqul1 = [base for base in check_base if countDic[base] == 1]
    allelenum = len(gtListAllele)
    allelereadnum = sum(countDic[base] for base in gtListAllele)
    altGT = ",".join(f"{base}[{countDic[base]}]" for base in snpintersection | indelintersection)
    mutypeindel = "INDEL" if indelintersection else "None"
    mutypesnp = "SNP" if snpintersection else "None"
    allelecheck = "MultiAlleles" if (len(gtListAlleleEqul1) / (len(gtListAlleleEqul1) + allelereadnum) >= 0.2 and len(gtListAlleleEqul1) + allelereadnum >= 2) or (len(gtListAlleleEqul1) >= 2 or allelenum >= 2) else "None"
    readsupportcheck = "PASS_reads" if allelenum != 0 else "NonePASS_reads"
    return mutypeindel, mutypesnp, allelecheck, readsupportcheck, altGT


# filter vcfile without barcode
def filter_VCF_withoutCB(vcfilePathIn, vcfilePathOut, binsize, clustersize, dpsize, altsize):
    vcfilePathOutmp0 = vcfilePathOut + ".f0"
    vcfilePathOutmp1 = vcfilePathOut + ".f1"
    # Remove existing output files if they exist
    for path in [vcfilePathOutmp0, vcfilePathOutmp1]:
        if os.path.exists(path):
            try:
                subprocess.run(['rm', path], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error removing file {path}: {e}")
    cmd_sort = f"cat {vcfilePathIn} | sort -k1,1 -k2,2n > {vcfilePathOutmp0}"
    try:
        subprocess.run(cmd_sort, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command '{cmd_sort}': {e}")
    # Filter noisy cluster mutations
    posList = {}
    locExcludeTotal = {}
    with open(vcfilePathOutmp0, "r") as vcfileIn:
        for line in vcfileIn:
            chr, pos = line.split("\t")[0], line.split("\t")[1]
            posList.setdefault(chr, []).append(pos)
            locExcludeTotal.setdefault(chr, [])
            if len(posList[chr]) > 1 and int(pos) - int(posList[chr][0]) <= binsize:
                locExcludeTotal[chr].extend([pos, posList[chr][0]])
            posList[chr] = [pos]
    locExclude = {chr: filter_pos_binCluster(posList, binsize, clustersize) for chr, posList in locExcludeTotal.items()}
    # Filter based on depth and alternative size
    with open(vcfilePathOutmp0, "r") as vcfileIn, open(vcfilePathOutmp1, "w") as vcfileOut:
        for line in vcfileIn:
            dp, alt = int(line.split("\t")[3]), int(line.split("\t")[5])
            if dp < dpsize or alt < altsize:
                continue
            chr, pos, refBase, sitmutinfo = line.split("\t")[0], line.split("\t")[1], line.split("\t")[2], line.split("\t")[4]
            m1, m2, m3, m4, m5 = assign_mutype_snpindex(sitmutinfo, refBase)
            cluster_status = "in_cluster" if int(pos) in locExclude.get(chr, []) else "not_in_cluster"
            vcfileOut.write(f"{line.strip()}\t{cluster_status}\t{m1}\t{m2}\t{m3}\t{m4}\t{m5}\n")
    # Remove temporary files
    for path in [vcfilePathOutmp0]:
        if os.path.exists(path):
            try:
                subprocess.run(['rm', path], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error removing file {path}: {e}")

# generate cell type dic
def generate_cellType_dic(barcodePath, excludeList):
    cellTypeDic = {}
    with open(barcodePath, "r") as barcodeIn:
        for line in barcodeIn:
            numCheck = line.split("\t")
            if len(numCheck) >= 2:
                barcodetmp = numCheck[0].strip()
                cellTypetmp = numCheck[1].strip()
            else:
                barcodetmp = numCheck[0].strip()
                cellTypetmp = "CellType1"
            if cellTypetmp in excludeList:
                continue
            if cellTypetmp not in cellTypeDic:
                cellTypeDic[cellTypetmp] = [barcodetmp]
            else:
                cellTypeDic[cellTypetmp].append(barcodetmp)
    return cellTypeDic

# check one cell if mut
def check_ifMut_oneCell(cellInfo, refBase):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    countDic = {base: 0 for base in basicBase}
    check_base = list(set(basicBase) - {refBase})
    cellInfoList = list(map(int, cellInfo.split("|")))
    for base, count in zip(basicBase, cellInfoList):
        countDic[base] += count
    numCheck = sum(countDic[base] for base in check_base)
    return numCheck >= 1

# get one cell genotype  
def check_MutGT_oneCell(cellInfo, refBase):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    countDic = {base: 0 for base in basicBase}
    check_base = list(set(basicBase) - {refBase})
    cellInfoList = list(map(int, cellInfo.split("|")))
    for base, count in zip(basicBase, cellInfoList):
        countDic[base] += count
    genotype = "|".join(f"{base}-{countDic[base]}" for base in check_base if countDic[base] != 0)
    return genotype

# generate cell type count dic
def generate_cellType_countDic(cellTypeDic, totalCellBarcodeList, mutcellBarcodeList, mutcellMutInfoList, refBase):
    CellTypeCountDic = {}
    altCellType = []
    MultiCellType = "None"
    for key in cellTypeDic.keys():
        interactionCellBarcodeMut = list(set(cellTypeDic[key]).intersection(set(mutcellBarcodeList)))
        interactionCellBarcodeTotal = list(set(cellTypeDic[key]).intersection(set(totalCellBarcodeList)))
        interactionCellMutInfoMut = ":".join([f"{i}-{check_MutGT_oneCell(mutcellMutInfoList[mutcellBarcodeList.index(i)], refBase)}" for i in interactionCellBarcodeMut])
        numMut = len(interactionCellBarcodeMut)
        if numMut > 0:
            altCellType.append(key)
        numTotal = len(interactionCellBarcodeTotal)
        CellTypeCountDic[key] = f"{numMut}|{numTotal};{interactionCellMutInfoMut}"
    if len(altCellType) >= 2:
        MultiCellType = "MultiCellType"
    return CellTypeCountDic, altCellType, MultiCellType

def integrate_baseQuality_onesite(oneSiteBaseQua):
    # Example: 32||||20|;27:32||||20|;32:37|||||
    basicBase = ["A", "T", "C", "G", "D", "I"]
    baseQuaDic = {base: [] for base in basicBase}
    cellsBaseQuaList = oneSiteBaseQua.split(";")
    for cell in cellsBaseQuaList:
        cellBaseList = cell.split("|")
        for index_base, cellBase in enumerate(cellBaseList):
            if cellBase:
                base = basicBase[index_base]
                baseQuaDic[base].append(cellBase)
    integrate_baseQuality = "|".join(":".join(baseQuaDic[base]) for base in basicBase)
    return integrate_baseQuality

# cal lowquality base
def split_baseQua_lowqualitynum(refBase, oneCellBaseQua, lowBaseCutoff):
    # Example: 32:32||||20|
    basicBase = ["A", "T", "C", "G", "D", "I"]
    baseQuainfoList = oneCellBaseQua.split("|")
    lowrefBaseNum = len([int(i) for i in baseQuainfoList[basicBase.index(refBase)].split(":") if i and int(i) <= lowBaseCutoff])
    altBase = [base for base in basicBase if base != refBase and baseQuainfoList[basicBase.index(base)]]
    lowaltBaseNum = sum(len([int(j) for j in baseQuainfoList[basicBase.index(base)].split(":") if j and int(j) <= lowBaseCutoff]) for base in altBase)
    return lowrefBaseNum, lowaltBaseNum

# get cell mutation info 
def get_cellmutInfo_forTest(vcfPathIn, vcfPathOut, cellTypeDic, cellsize, altcellsize, lowBaseCutoff):
    if os.path.exists(vcfPathOut):
        try:
            subprocess.run(['rm', vcfPathOut], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error removing file {vcfPathOut}: {e}")
    vcfIn = open(vcfPathIn, "r")
    vcfOut = open(vcfPathOut, "w")
    sqlHeaderList = ["chr", "pos", "refbase", "altbase", "mode", "altCellType", "totalCellNumber", "totalReadsNumber", "isMultiAlleles", "isMultiCellType", "isCluster", "SNP", "INDEL", "otherFilter"] + [i for i in cellTypeDic.keys()]
    sqlHeader = "\t".join(sqlHeaderList)
    vcfOut.writelines(sqlHeader + "\n")
    for line in vcfIn:
        linsets = line.split("\t")
        totalCellBarcodeList = []
        mutCellBarcodeList = []
        mutCellMutInfoList = []
        chr = linsets[0].strip()
        pos = linsets[1].strip()
        refBase = linsets[2].strip()
        totalReadsnum = linsets[3].strip()
        altReadsnum = linsets[5].strip()
        cellbarcodeList = linsets[9].split(";")
        cellmutInfoList = linsets[10].split(";")
        mode = linsets[19]
        totalCellnum = len(cellbarcodeList)
        MultiAllelesCheck = linsets[25].strip()
        SNPCheck = linsets[24].strip()
        INDELCheck = linsets[23].strip()
        altGT = linsets[27].strip()
        clusterCheck = linsets[22].strip()
        baseQualityInfoOnesite = linsets[11].strip()
        integrate_baseQua = integrate_baseQuality_onesite(baseQualityInfoOnesite)
        lowwtBasenum, lowmutBasenum = split_baseQua_lowqualitynum(refBase, integrate_baseQua, lowBaseCutoff)
        # fisher test for lowQualityBases
        p_fisher_lowBase = stats.fisher_exact([[int(lowmutBasenum), int(lowwtBasenum)], [int(altReadsnum) - int(lowmutBasenum), int(totalReadsnum) - int(altReadsnum) - int(lowwtBasenum)]])[1]
        if p_fisher_lowBase <= 0.05:
            tag_fisher_lowBase = "NonePASS_fisherLB"
        else:
            tag_fisher_lowBase = "PASS_fisherLB"
        frq_lowMutbase_totalMutbase = int(lowmutBasenum)/int(altReadsnum)
        if int(totalReadsnum) - int(altReadsnum) > 0:
            frq_lowWTbase_totalMutbase = int(lowwtBasenum)/(int(totalReadsnum) - int(altReadsnum))
        else:
            frq_lowWTbase_totalMutbase = int(lowwtBasenum)/(int(totalReadsnum) - int(altReadsnum) + 1)
        if frq_lowMutbase_totalMutbase >= 0.2 or frq_lowWTbase_totalMutbase >= 0.2:
            tag_lowBaseQuality = "NonePASS_NLB"
        else:
            tag_lowBaseQuality = "PASS_NLB"
        otherfilterInfo = ";".join([linsets[21].strip(), tag_fisher_lowBase, tag_lowBaseQuality])
        if totalCellnum < cellsize:
            continue
        mutCellnum = 0
        for i in range(0, totalCellnum):
            totalCellBarcodeList.append(cellbarcodeList[i])
            checkMut = check_ifMut_oneCell(cellmutInfoList[i], refBase)
            if checkMut:
                mutCellnum += 1
                mutCellBarcodeList.append(cellbarcodeList[i])
                mutCellMutInfoList.append(cellmutInfoList[i])
        if mutCellnum < altcellsize:
            continue
        CellTypeCountDic, altCellType, MultiCellType = generate_cellType_countDic(cellTypeDic, totalCellBarcodeList, mutCellBarcodeList, mutCellMutInfoList, refBase)
        sqlinlineList = [chr, pos, refBase, altGT, mode, ",".join(altCellType), str(mutCellnum) + "|" + str(totalCellnum), str(altReadsnum) + "|" + str(totalReadsnum), MultiAllelesCheck, MultiCellType, clusterCheck, SNPCheck, INDELCheck, otherfilterInfo] + [i for i in CellTypeCountDic.values()]
        sqlinline = "\t".join(sqlinlineList)
        vcfOut.writelines(sqlinline + "\n")


def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outpathtmp = args.outdir
    id = args.id
    outpath = outpathtmp + "/" + id
    vcfpath = args.vcfpath
    barcodePath = args.barcode
    binsize = args.binsize
    clustersize = args.clustersize
    dpsize = args.dpsize
    altsize = args.altdpsize
    cellsize = args.cellsize
    altcellsize = args.altcellsize
    lowBaseCutoff = args.lowBaseCutoff
    CellTypeDictmp = generate_cellType_dic(barcodePath, ["CellType"])
    finalOutPath = outpath + ".scmutrace.clean.vcf"
    filter_VCF_withoutCB(vcfpath, outpath, binsize, clustersize, dpsize, altsize)
    get_cellmutInfo_forTest(outpath + ".f1", finalOutPath, CellTypeDictmp, cellsize, altcellsize, lowBaseCutoff)
    cmdrm = "rm " + outpath + ".f1"
    try:
        subprocess.run(cmdrm, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command '{cmdrm}': {e}")


if __name__=='__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')