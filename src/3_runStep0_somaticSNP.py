import os
import timeit
import argparse
import scipy.stats as stats
#import multiprocessing as mp

def initialize_parser():
    parser = argparse.ArgumentParser(description='filter SNP based on public database')
    parser.add_argument('--outpath', type=str, default=1, help='out file path', required = True)
    parser.add_argument('--vcfpath', type=str, default=1, help='Vcf path', required = True)
    parser.add_argument('--id', type=str, default=1, help='sample id', required = True)
    parser.add_argument('--barcode', type=str, default=1, help='barcode file [tab separated file with header: Barcode CellType]', required = True)
    parser.add_argument('--mode', type=str, help='keep Weak or Strong evidence candidate snps [strong or weak]', required = False)
    parser.add_argument('--binsize', type=int, default=5, help='base number between two snps (MultiAlleles)', required = False)
    parser.add_argument('--clustersize', type=int, default=4, help='number of MultiAlleles snps', required = False)
    parser.add_argument('--dpsize', type=int, default=10, help='min depth', required = False)
    parser.add_argument('--altdpsize', type=int, default=5, help='min alternative depth', required = False)
    parser.add_argument('--cellsize', type=int, default=10, help='min cell number', required = False)
    parser.add_argument('--altcellsize', type=int, default=5, help='min mut cell number', required = False)
    parser.add_argument('--bamcoverage', type=str, default=1, help='coverage file of bam', required = False)
    parser.add_argument('--distance', type=int, default=10, help='distance between mut and start/end pos', required = False)
    parser.add_argument('--sexcheck', action='store_true', help='A boolean flag for exclude sex chromosome')
    return (parser)

def err_call_back(err):
    print(f'apply_async error:{str(err)}')

# filter pos list based on bin size and cluster size
def filter_pos_binCluster(posList, binsize, clustersize):
    posListSort = [int(i) for i in list(set(posList))]
    posListSort.sort()
    excludeList = []
    posListmp = []
    for pos in posListSort:
        pos = int(pos)
        if len(posListmp) == 0:
            posListmp.append(pos)
        elif pos - posListmp[-1] <= binsize:
            posListmp.append(pos)
        elif pos - posListmp[-1] > binsize and len(posListmp) >= clustersize:
            excludeList += posListmp
            posListmp.clear()
            posListmp.append(pos)
        elif pos - posListmp[-1] > binsize and len(posListmp) < clustersize:
            posListmp.clear()
            posListmp.append(pos)
        else:
            pass
    if len(posListmp) >= clustersize:
        excludeList += posListmp
    return excludeList

# assign indel and snp
def assign_mutype_snpindex(sitemutinfo, refBase):
    indelBase = ["D", "I"]
    snpBase = ["A", "T", "C", "G"]
    basicBase = ["A", "T", "C", "G", "D", "I"]
    countDic = {"A": 0, "T": 0, "C": 0, "G": 0, "D": 0, "I": 0}
    check_base = list(set(basicBase).difference(set([refBase])))
    sitemutinfoList = sitemutinfo.split(":")
    countDic['A'] += int(sitemutinfoList[0])
    countDic['T'] += int(sitemutinfoList[1])
    countDic['C'] += int(sitemutinfoList[2])
    countDic['G'] += int(sitemutinfoList[3])
    countDic['D'] += int(sitemutinfoList[4])
    countDic['I'] += int(sitemutinfoList[5])
    gtList = [i for i in check_base if countDic[i] != 0]
    indelintersection = list(set(gtList).intersection(set(indelBase)))
    snpintersection = list(set(gtList).intersection(set(snpBase)))
    gtListAllele = [i for i in check_base if countDic[i] > 1]
    gtListAlleleEqul1 = [i for i in check_base if countDic[i] == 1]
    allelenum = len(gtListAllele)
    altGT = ",".join(list(i + "[" + str(countDic[i]) + "]" for i in list(set(snpintersection + indelintersection))))
    mutypeindel = "None"
    mutypesnp = "None"
    allelecheck = "None"
    readsupportcheck = "PASS_reads"
    if len(indelintersection) >= 1:
        mutypeindel = "INDEL"
    if len(snpintersection) >= 1:
        mutypesnp = "SNP"
    if len(gtListAlleleEqul1) >= 2 or allelenum >= 2:
        allelecheck = "MultiAlleles"
    #if len(gtListAlleleEqul1) + allelenum >= 2:
    #    allelecheck = "MultiAlleles"
    if allelenum == 0:
        readsupportcheck = "NonePASS_reads"
    return mutypeindel, mutypesnp, allelecheck, readsupportcheck, altGT

# filter vcfile without barcode
def filter_VCF_withoutCB(vcfilePathIn, vcfilePathOut, mode, binsize, clustersize, dpsize, altsize, sexcheck):
    vcfilePathOutmp0 = vcfilePathOut + ".0.fNonePass_sexcheck"
    vcfilePathOutmp1 = vcfilePathOut + ".1.fNormal"
    vcfilePathOutmp2 = vcfilePathOut + ".2.fdpsize"
    if os.path.exists(vcfilePathOut):
        cmd = "rm " + vcfilePathOut
        os.system(cmd)
    if os.path.exists(vcfilePathOutmp0):
        cmd = "rm " + vcfilePathOutmp0
        os.system(cmd)
    if os.path.exists(vcfilePathOutmp1):
        cmd = "rm " + vcfilePathOutmp1
        os.system(cmd)
    # sex check
    if sexcheck:
        cmd0 = "cat " + vcfilePathIn + " | grep -v X | grep -v Y | sort -k1 -k2 -n > " + vcfilePathOutmp0
        os.system(cmd0)
    else:
        cmd0 = "cat " + vcfilePathIn + " | sort -k1 -k2 -n > " + vcfilePathOutmp0
        os.system(cmd0)
    # filter strong or weak
    if mode == "strong":
        print("strong mode...")
        cmd = "cat " + vcfilePathOutmp0 + "| grep Strong | sort -k1 -k2 -n > " + vcfilePathOutmp1
        os.system(cmd)
    if mode == "weak":
        print("weak mode...")
        cmd = "cat " + vcfilePathOutmp0 + "| sort -k1 -k2 -n > " + vcfilePathOutmp1
        os.system(cmd)
    # filter noisy cluster mutations
    vcfileIn = open(vcfilePathOutmp1, "r")
    posList = {}
    locExcludeTotal = {}
    for line in vcfileIn:
        chr = line.split("\t")[0]
        if chr not in posList.keys():
            posList[chr] = []
        if chr not in locExcludeTotal.keys():
            locExcludeTotal[chr] = []
        pos = line.split("\t")[1]
        # get locs will be excluded
        if len(posList[chr]) == 0:
            posList[chr] += [pos]
        else:
            if int(pos) - int(posList[chr][0]) <= binsize:
                locExcludeTotal[chr] += [pos]
                locExcludeTotal[chr] += [posList[chr][0]]
            else:
                pass
            posList[chr].pop()
            posList[chr].append(pos)
    vcfileIn.close()
    locExclude = {}
    #print(locExcludeTotal)
    for exclude_chr, exclude_posList in locExcludeTotal.items():
        if exclude_chr not in locExclude.keys():
            locExclude[exclude_chr] = []
        locExclude[exclude_chr] = filter_pos_binCluster(exclude_posList, binsize, clustersize)
    # filter 
    vcfileIn = open(vcfilePathOutmp1, "r")
    vcfileOut = open(vcfilePathOutmp2, "w")
    for line in vcfileIn:
        dp = line.split("\t")[3]
        alt = line.split("\t")[5]
        if int(dp) < int(dpsize):
            continue
        if int(alt) < int(altsize):
            continue
        chr = line.split("\t")[0]
        pos = line.split("\t")[1]
        refBase = line.split("\t")[2]
        sitmutinfo = line.split("\t")[4]
        m1, m2, m3, m4, m5,= assign_mutype_snpindex(sitmutinfo, refBase)
        if int(pos) in locExclude[chr]:
            vcfileOut.writelines(line.strip() + "\t" + "in_cluster" + "\t" + m1 + "\t" + m2 + "\t" + m3 + "\t" + m4 + "\t" + m5 +"\n")
        else:
            vcfileOut.writelines(line.strip() + "\t" + "not_in_cluster" + "\t" + m1 + "\t" + m2 + "\t" + m3 + "\t" + m4 + "\t" + m5 + "\n")
    vcfileIn.close()
    vcfileOut.close()

# filter bam coverage
def filter_bam_coverage(vcfPathIn, vcfPathOut, bamCoveragePath, distance):
    vcfilePathOutmp1 = vcfPathOut + ".3.fcoverage1"
    vcfilePathOutmp2 = vcfPathOut + ".3.fcoverage2"
    vcfilePathOutmp3 = vcfPathOut + ".3.fcoverage"
    if os.path.exists(vcfilePathOutmp1):
        cmd = "rm " + vcfilePathOutmp1
        os.system(cmd)
    if os.path.exists(vcfilePathOutmp2):
        cmd = "rm " + vcfilePathOutmp2
        os.system(cmd)
    if os.path.exists(vcfilePathOutmp3):
        cmd = "rm " + vcfilePathOutmp3
        os.system(cmd)
    outfiletmp = open(vcfilePathOutmp1, "w")
    # process dpsize vcf
    with open(vcfPathIn, "r") as fileIn:
        for line in fileIn:
            linsets = line.split("\t")
            sql = linsets[0] + "\t" + linsets[1] + "\t" + linsets[1]
            outfiletmp.writelines(sql + "\n")
    fileIn.close()
    outfiletmp.close()
    cmd = "bedtools intersect -a " + bamCoveragePath + " -b " + vcfilePathOutmp1 + " -wa -wb > " + vcfilePathOutmp2
    os.system(cmd)
    excludePoslist = []
    #outfiletmp3 = open(vcfilePathOutmp3, "w")
    with open(vcfilePathOutmp2, "r") as fileIn:
        for line in fileIn:
            lineSets = line.split("\t")
            chrtmp = lineSets[0]
            startpos = int(lineSets[1])
            endpos = int(lineSets[2])
            mutpos = int(lineSets[4])
            dsm = abs(startpos - mutpos)
            dem = abs(endpos - mutpos)
            if dsm < distance or dem < distance:
                excludePoslist.append(str(chrtmp) + "_" + str(mutpos))
    fileIn.close()
    fileOut = open(vcfilePathOutmp3, "w")
    with open(vcfPathIn, "r") as fileIn:
        for line in fileIn:
            lineSets = line.split('\t')
            chrtchecktmp = lineSets[0]
            mutposchecktmp = lineSets[1]
            checksql = str(chrtchecktmp) + "_" + str(mutposchecktmp)
            if checksql in excludePoslist:
                dCheck = "near_pair"
            else:
                dCheck = "not_near_pair"
            sqlout = line.strip() + "\t" + dCheck + "\n"
            fileOut.writelines(sqlout)
    fileIn.close()
    fileOut.close()
    cmdrm = "rm " + vcfilePathOutmp1 + " " + vcfilePathOutmp2
    os.system(cmdrm)

# generate cell type dic
def generate_cellType_dic(barcodePath, excludeList):
    cellTypeDic = {}
    barcodeIn = open(barcodePath, "r")
    for line in barcodeIn:
        barcodetmp = line.split("\t")[0].strip()
        cellTypetmp = line.split("\t")[1].strip()
        if cellTypetmp in excludeList:
            continue
        if cellTypetmp not in cellTypeDic.keys():
            cellTypeDic[cellTypetmp] = [barcodetmp]
        else:
            cellTypeDic[cellTypetmp] += [barcodetmp]
    return cellTypeDic

# check one cell if mut
def check_ifMut_oneCell(cellInfo, refBase):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    countDic = {"A": 0, "T": 0, "C": 0, "G": 0, "D": 0, "I": 0}
    check_base = list(set(basicBase).difference(set([refBase])))
    cellInfoList = cellInfo.split("|")
    countDic['A'] += int(cellInfoList[0])
    countDic['T'] += int(cellInfoList[1])
    countDic['C'] += int(cellInfoList[2])
    countDic['G'] += int(cellInfoList[3])
    countDic['D'] += int(cellInfoList[4])
    countDic['I'] += int(cellInfoList[5])
    numCheck = sum([int(countDic[i]) for i in countDic.keys() if i in check_base])
    if numCheck >= 1:
        return True
    else:
        return False

# get one cell genotype  
def check_MutGT_oneCell(cellInfo, refBase):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    countDic = {"A": 0, "T": 0, "C": 0, "G": 0, "D": 0, "I": 0}
    check_base = list(set(basicBase).difference(set([refBase])))
    cellInfoList = cellInfo.split("|")
    countDic['A'] += int(cellInfoList[0])
    countDic['T'] += int(cellInfoList[1])
    countDic['C'] += int(cellInfoList[2])
    countDic['G'] += int(cellInfoList[3])
    countDic['D'] += int(cellInfoList[4])
    countDic['I'] += int(cellInfoList[5])
    genotype = "|".join([i + "-" + str(countDic[i]) for i in countDic.keys() if i in check_base and countDic[i] != 0])
    return genotype

# generate cell type count dic
def generate_cellType_countDic(cellTypeDic, totalCellBarcodeList, mutcellBarcodeList, mutcellMutInfoList, refBase):
    CellTypeCountDic = {}
    altCellType = []
    MultiCellType = "None"
    for key in cellTypeDic.keys():
        interactionCellBarcodeMut = list(set(cellTypeDic[key]).intersection(set(mutcellBarcodeList)))
        interactionCellBarcodeTotal = list(set(cellTypeDic[key]).intersection(set(totalCellBarcodeList)))
        interactionCellMutInfoMut = ":".join(list([i + "-" + check_MutGT_oneCell(mutcellMutInfoList[mutcellBarcodeList.index(i)], refBase) for i in interactionCellBarcodeMut]))
        numMut = len(interactionCellBarcodeMut)
        if numMut > 0:
            altCellType.append(key)
        numTotal = len(interactionCellBarcodeTotal)
        CellTypeCountDic[key] = str(numMut) + "|" + str(numTotal) + ";" + interactionCellMutInfoMut
    if len(altCellType) >= 2:
        MultiCellType = "MultiCellType"
    return CellTypeCountDic, altCellType, MultiCellType

def integrate_baseQuality_onesite(oneSiteBaseQua):
    # Example: 32||||20|;27:32||||20|;32:37|||||
    basicBase = ["A", "T", "C", "G", "D", "I"]
    baseQuaDic = {"A": [], "T": [], "C": [], "G": [], "D": [], "I": []}
    cellsBaseQuaList = oneSiteBaseQua.split(";")
    for cell in cellsBaseQuaList:
        cellBaseList = cell.split("|")
        for index_base, cellBase in enumerate(cellBaseList):
            if cellBase != "":
                indexBase = basicBase[index_base]
                baseQuaDic[indexBase].append(cellBase)
    integrate_baseQuality = "|".join([":".join(baseQuaDic["A"]), ":".join(baseQuaDic["T"]), ":".join(baseQuaDic["C"]), ":".join(baseQuaDic["G"]), ":".join(baseQuaDic["D"]), ":".join(baseQuaDic["I"])])
    return integrate_baseQuality

# cal lowquality base
def split_baseQua_lowqualitynum(refBase, oneCellBaseQua, lowBaseCutoff):
    # Example: 32:32||||20|
    basicBase = ["A", "T", "C", "G", "D", "I"]
    baseQuainfoList = oneCellBaseQua.split("|")
    base = []
    base = [basicBase[i] for i in range(0, 6) if baseQuainfoList[i] != ""]
    lowrefBaseNum =len([int(i) for i in baseQuainfoList[basicBase.index(refBase)].split(":") if i != "" and int(i) <= lowBaseCutoff])
    altBase = list(set(base).difference(set([refBase])))
    lowaltBaseNum = sum([int(len([int(j) for j in baseQuainfoList[basicBase.index(i)].split(":") if j != "" and int(j) <= lowBaseCutoff])) for i in altBase])
    return lowrefBaseNum, lowaltBaseNum

# get cell mutation info 
def get_cellmutInfo_forTest(vcfPathIn, vcfPathOut, cellTypeDic, cellsize, altcellsize):
    vcfIn = open(vcfPathIn, "r")
    if os.path.exists(vcfPathOut):
        cmdrm = "rm " + vcfPathOut
        os.system(cmdrm)
    vcfOut = open(vcfPathOut, "w")
    sqlHeaderList = ["chr", "start", "end", "refbase", "altbase", "mode", "altCellType", "totalCellNumber", "totalReadsNumber", "isMultiAlleles", "isMultiCellType", "isCluster", "isCOSMIC","SNP", "INDEL", "isNear", "otherFilter"] + [i for i in cellTypeDic.keys()]
    sqlHeader = "\t".join(sqlHeaderList)
    vcfOut.writelines(sqlHeader + "\n")
    for line in vcfIn:
        linsets = line.split("\t")
        if linsets[-8] == "in_normal" or linsets[-9] == "in_paired_normal":
            continue
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
        MultiAllelesCheck = linsets[-4].strip()
        SNPCheck = linsets[-5].strip()
        INDELCheck = linsets[-6].strip()
        altGT = linsets[-2].strip()
        clusterCheck = linsets[-7].strip()
        COSMICcheck = linsets[-10].strip()
        nearCheck = linsets[-1].strip()
        lowBaseCutoff = 20
        baseQualityInfoOnesite = linsets[11].strip()
        integrate_baseQua = integrate_baseQuality_onesite(baseQualityInfoOnesite)
        lowwtBasenum, lowmutBasenum = split_baseQua_lowqualitynum(refBase, integrate_baseQua, lowBaseCutoff)
        # fisher test for lowQualityBases
        p_fisher_lowBase = stats.fisher_exact([[int(lowmutBasenum), int(lowwtBasenum)], [int(altReadsnum) - int(lowmutBasenum), int(totalReadsnum) - int(altReadsnum) - int(lowwtBasenum)]]).pvalue
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
        otherfilterInfo = ";".join([linsets[-3].strip(), linsets[-11].strip(), linsets[-12].strip(), linsets[-13].strip(), linsets[-14].strip(), linsets[-15].strip(), linsets[-16].strip(), linsets[-17].strip(), linsets[-18].strip(), tag_fisher_lowBase, tag_lowBaseQuality])
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
        sqlinlineList = [chr, pos, pos, refBase, altGT, mode, ",".join(altCellType), str(mutCellnum) + "|" + str(totalCellnum), str(altReadsnum) + "|" + str(totalReadsnum), MultiAllelesCheck, MultiCellType, clusterCheck, COSMICcheck, SNPCheck, INDELCheck, nearCheck, otherfilterInfo] + [i for i in CellTypeCountDic.values()]
        sqlinline = "\t".join(sqlinlineList)
        vcfOut.writelines(sqlinline + "\n")

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outpathtmp = args.outpath
    id = args.id
    outpath = outpathtmp + "/" + id
    vcfpath = args.vcfpath
    barcodePath = args.barcode
    mode = args.mode
    binsize = args.binsize
    clustersize = args.clustersize
    dpsize = args.dpsize
    altsize = args.altdpsize
    cellsize = args.cellsize
    altcellsize = args.altcellsize
    sexcheck = args.sexcheck
    bamCoveragePath = args.bamcoverage
    distance = args.distance
    CellTypeDictmp = generate_cellType_dic(barcodePath, ["CellType"])
    finalOutPath = outpath + ".somatic.vcf"
    filter_VCF_withoutCB(vcfpath, outpath, mode, binsize, clustersize, dpsize, altsize, sexcheck)
    filter_bam_coverage(outpath + ".2.fdpsize", outpath, bamCoveragePath, distance)
    get_cellmutInfo_forTest(outpath + ".3.fcoverage", finalOutPath, CellTypeDictmp, cellsize, altcellsize)
    cmdrm = "rm " + outpath + ".0.fNonePass_sexcheck " + outpath + ".1.fNormal " + outpath + ".2.fdpsize " + outpath + ".3.fcoverage"
    os.system(cmdrm)

if __name__=='__main__':
    start = timeit.default_timer()
    #19_1_59128983
    k = main()
    print(k)
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')
