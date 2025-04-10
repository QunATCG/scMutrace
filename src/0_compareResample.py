import random
import math
import timeit
import argparse

def initialize_parser():
    parser = argparse.ArgumentParser(description='filter SNP based on 3Model')
    parser.add_argument('--outpath', type=str, default=1, help='outpath', required = True)
    parser.add_argument('--vcfpath', type=str, default=1, help='vcfpath', required = True)
    return (parser)

def cal_pError_fromBaseQuality(baseQuality):
    pError = 10 ** (-int(baseQuality)/10)
    return(pError)

def define_accuracy_calling(baseQuality):
    # Default 10x single cell mapping quality is 255. Also, should keep accuracy != 0
    # in our analysis, we filtered reads with low quality of mapping
    callingErrorP = cal_pError_fromBaseQuality(baseQuality)
    accuracySequencingP = 1 - callingErrorP
    if accuracySequencingP == 1:
        accuracySequencingP = 0.999999
    return accuracySequencingP

def shuffle_reads_features(refBase,baseQualityInfoOnesite):
    # Example input: 32||||20|;27:32||||20|;32:37|||||;22:41:32:32:37:27|||||;41|||||;|||27|20|
    basicBase = ["A", "T", "C", "G", "D", "I"]
    cellNum = []
    tmpDic = {"A": [], "T": [], "C":[], "G": [], "D": [], "I":[]}
    tmpDicP = {"A": [], "T": [], "C":[], "G": [], "D": [], "I":[]}
    num = len(baseQualityInfoOnesite.split(";"))
    for i in range(0, num):
        baseQuaOneCell = baseQualityInfoOnesite.split(";")[i].split("|")
        readcount = 0
        for index_base, j in enumerate(baseQuaOneCell):
            if j != "":
                readcount += len(j.split(":"))
                tmpDic[basicBase[index_base]] += [m for m in j.split(":")]
                #tmpDicP[basicBase[index_base]] += [define_accuracy_calling(int(m)) for m in j.split(":")]
                tmpDicP[basicBase[index_base]] += [int(m) for m in j.split(":")]
        cellNum.append(readcount)
    baseWeight = [len(value) for value in tmpDic.values()]
    baseNumTotal = sum(baseWeight)
    # check baseWeight, when there is no refbase, should consider the refbase exist
    if baseWeight[basicBase.index(refBase)] == 0:
        refBaseWeight = cal_pError_fromBaseQuality(min(tmpDic["A"] + tmpDic["T"] + tmpDic["C"] + tmpDic["G"] + tmpDic["D"] + tmpDic["I"]))
        baseWeight[basicBase.index(refBase)] = refBaseWeight * baseNumTotal
        tmpDicP[refBase] = [1]
    readcountsNum = list(set(cellNum))
    readcountsWeight = [cellNum.count(i) for i in readcountsNum]
    readcountsNumNormal = readcountsNum
    if max(readcountsNum) >= 10:
        readcountsNumNormal = [math.ceil((10 * i) / max(readcountsNum)) for i in readcountsNum]
    return readcountsNumNormal, baseWeight, tmpDic, tmpDicP, readcountsWeight

def shuffle_reads_baseonFeatures(readcountsNum, baseWeight, tmpDic, tmpDicP, readcountsWeight):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    readsNumSelect = random.choices(readcountsNum, weights = readcountsWeight, k = 1)
    baseSelect = random.choices(basicBase, weights = baseWeight, k = readsNumSelect[0])
    baseQuaA = []
    baseQuaT = []
    baseQuaC = []
    baseQuaG = []
    baseQuaD = []
    baseQuaI = []
    if len(tmpDic["A"]) != 0:
        baseQuaA = random.choices(tmpDic["A"], weights = tmpDicP["A"], k = baseSelect.count("A"))
    if len(tmpDic["T"]) != 0:
        baseQuaT = random.choices(tmpDic["T"], weights = tmpDicP["T"], k = baseSelect.count("T"))
    if len(tmpDic["C"]) != 0:
        baseQuaC = random.choices(tmpDic["C"], weights = tmpDicP["C"], k = baseSelect.count("C"))
    if len(tmpDic["G"]) != 0:
        baseQuaG = random.choices(tmpDic["G"], weights = tmpDicP["G"], k = baseSelect.count("G"))
    if len(tmpDic["D"]) != 0:
        baseQuaD = random.choices(tmpDic["D"], weights = tmpDicP["D"], k = baseSelect.count("D"))
    if len(tmpDic["I"]) != 0:
        baseQuaI = random.choices(tmpDic["I"], weights = tmpDicP["I"], k = baseSelect.count("I"))
    integrate_baseQuality = "|".join([":".join(baseQuaA), ":".join(baseQuaT), ":".join(baseQuaC), ":".join(baseQuaG), ":".join(baseQuaD), ":".join(baseQuaI)])
    return integrate_baseQuality

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

def shuffle_reads_simple(baseQualityInfoOnesite, times):
    cells = baseQualityInfoOnesite.split(";")
    cell = random.choices(cells, k = times)
    return cell

def countOnecellAverageBQ(cellinfo):
    tmp0 = []
    tmp1 = []
    for i in cellinfo.split("|"):
        tmp0 += i.split(":")
    for i in tmp0:
        if i != "":
            tmp1.append(int(i))
    return sum(tmp1)/len(tmp1)

def shuffle_reads_simple_weight(baseQualityInfoOnesite, times):
    cells = baseQualityInfoOnesite.split(";")
    weights = [countOnecellAverageBQ(cell) for cell in cells]
    cell = random.choices(cells, weights= weights, k = times)
    return cell

def countfraq(baseQualityInfoOnesite, thresholdBQ):
    # 41|||||41:41;||41|||;|||27||;||41|||;|||27||
    tmp0 = []
    tmp1 = []
    tmp2 = []
    for cellinfo in baseQualityInfoOnesite.split(";"):
        tmp0 += cellinfo.split("|")
    for i in tmp0:
        tmp1 += i.split(":")
    for i in tmp1:
        if i != "":
            tmp2.append(i)
    baselessthanthreshold = [i for i in tmp2 if int(i) <= int(thresholdBQ)]
    frq = 100 * len(baselessthanthreshold)/len(tmp2)    
    return frq

def get_base_qual(baseQualityInfoOnesite):
    # 41|||||41:41;||41|||;|||27||;||41|||;|||27||
    tmp0 = []
    tmp1 = []
    tmp2 = []
    for cellinfo in baseQualityInfoOnesite.split(";"):
        tmp0 += cellinfo.split("|")
    for i in tmp0:
        tmp1 += i.split(":")
    for i in tmp1:
        if i != "":
            tmp2.append(i)
    sql = ",".join(tmp2)
    return sql

def countCellfeature(cellinfo, refBase):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    altBase = list(set(basicBase).difference(set([refBase])))
    tmpDic = {"A": [], "T": [], "C":[], "G": [], "D": [], "I":[]}
    tmpDic["A"] += cellinfo.split("|")[0].split(":")
    tmpDic["T"] += cellinfo.split("|")[1].split(":")
    tmpDic["C"] += cellinfo.split("|")[2].split(":")
    tmpDic["G"] += cellinfo.split("|")[3].split(":")
    tmpDic["D"] += cellinfo.split("|")[4].split(":")
    tmpDic["I"] += cellinfo.split("|")[5].split(":")
    if tmpDic[refBase][0] != "":
        refCount = len(tmpDic[refBase])
    else:
        refCount = 0
    altCount = 0
    for i in altBase:
        if tmpDic[i][0] != "":
            tmpaltCount = len(tmpDic[i])
        else:
            tmpaltCount = 0
        altCount += tmpaltCount
    return refCount, altCount, refCount + altCount

def countOnesitefeature(oneSiteinfo, refBase):
    cells = oneSiteinfo.split(";")
    cells_reads = []
    refcount = 0
    altcount = 0
    for cell in cells:
        a, b, c = countCellfeature(cell, refBase)
        refcount += a
        altcount += b
        cells_reads.append(str(c))
    return str(refcount), str(altcount), str(len([i for i in cells_reads if int(i) <= 1]))

def check_ifMut_oneCell(cellInfo, refBase):
    basicBase = ["A", "T", "C", "G", "D", "I"]
    countDic = {"A": 0, "T": 0, "C": 0, "G": 0, "D": 0, "I": 0}
    check_base = list(set(basicBase).difference(set([refBase])))
    cellInfoList = cellInfo.split("|")
    if cellInfoList[0] != "":
        countDic['A'] += 1
    if cellInfoList[1] != "":
        countDic['T'] += 1
    if cellInfoList[2] != "":
        countDic['C'] += 1
    if cellInfoList[3] != "":
        countDic['G'] += 1
    if cellInfoList[4] != "":
        countDic['D'] += 1
    if cellInfoList[5] != "":
        countDic['I'] += 1
    numCheck = sum([int(countDic[i]) for i in countDic.keys() if i in check_base])
    if numCheck >= 1:
        return True
    else:
        return False
    
def count_mut_cellNum(baseQualityInfoOnesite, refBase):
    cells = baseQualityInfoOnesite.split(";")
    totalCellsNum = len(cells)
    altCellNum = 0
    for cell in cells:
        tag = check_ifMut_oneCell(cell, refBase)
        if tag:
            altCellNum += 1
    return str(altCellNum), str(totalCellsNum)

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    vcfpath = args.vcfpath
    outpath = args.outpath
    outfile = open(outpath, "w")
    with open(vcfpath, "r") as ff:
        for line in ff:
            if line != "":
                baseQ = 30
                lineSets = line.split("\t")
                readsNum = lineSets[3]
                if int(readsNum) >= 10:
                    refBase = lineSets[2]
                    # Sampling
                    # original reads 
                    baseQualityInfoOnesite = lineSets[11]

                    # simple resample
                    tmp_resampleSimple_100 = shuffle_reads_simple(baseQualityInfoOnesite, 100)
                    tmp_resampleSimple_1000 = shuffle_reads_simple(baseQualityInfoOnesite, 1000)
                    
                    # simple weight resample
                    tmp_resampleSimpleWei_100 = shuffle_reads_simple_weight(baseQualityInfoOnesite, 100)
                    tmp_resampleSimpleWei_1000 = shuffle_reads_simple_weight(baseQualityInfoOnesite, 1000)

                    # feature resample
                    tmp_resampleFeature_100 = []
                    tmp_resampleFeature_1000 = []
                    readcountsNum, baseWeight, tmpDic, tmpDicP, readcountsWeight = shuffle_reads_features(refBase,baseQualityInfoOnesite)
                    for i in range(0, 100):
                        sample_random_baseQuaInfo_tmp = shuffle_reads_baseonFeatures(readcountsNum, baseWeight, tmpDic, tmpDicP, readcountsWeight)
                        tmp_resampleFeature_100.append(sample_random_baseQuaInfo_tmp)
                    for i in range(0, 1000):
                        sample_random_baseQuaInfo_tmp = shuffle_reads_baseonFeatures(readcountsNum, baseWeight, tmpDic, tmpDicP, readcountsWeight)
                        tmp_resampleFeature_1000.append(sample_random_baseQuaInfo_tmp)

                    # Onesiteinfo
                    # original
                    onesiteinfo_original = baseQualityInfoOnesite
                    # simple resample
                    onesiteinfo_simple_100 = ";".join(tmp_resampleSimple_100)
                    onesiteinfo_simple_1000 = ";".join(tmp_resampleSimple_1000)
                    # simple resample
                    onesiteinfo_simpleWei_100 = ";".join(tmp_resampleSimpleWei_100)
                    onesiteinfo_simpleWei_1000 = ";".join(tmp_resampleSimpleWei_1000)
                    # feature resample
                    onesiteinfo_feature_100 = ";".join(tmp_resampleFeature_100)
                    onesiteinfo_feature_1000 = ";".join(tmp_resampleFeature_1000)

                    # Base quality
                    base_quality_original = str(countfraq(onesiteinfo_original, baseQ))
                    base_quality_simple100 = str(countfraq(onesiteinfo_simple_100, baseQ))
                    base_quality_simple1000 = str(countfraq(onesiteinfo_simple_1000, baseQ))
                    base_quality_simple100Wei = str(countfraq(onesiteinfo_simpleWei_100, baseQ))
                    base_quality_simple1000Wei = str(countfraq(onesiteinfo_simpleWei_1000, baseQ))
                    base_quality_feature100 = str(countfraq(onesiteinfo_feature_100, baseQ))
                    base_quality_feature1000 = str(countfraq(onesiteinfo_feature_1000, baseQ))

                    # Reads feature
                    reads_ref_original, reads_alt_original, reads_num_original = countOnesitefeature(onesiteinfo_original, refBase)
                    reads_ref_simple100, reads_alt_simple100, reads_num_simple100 = countOnesitefeature(onesiteinfo_simple_100, refBase)
                    reads_ref_simple1000, reads_alt_simple1000, reads_num_simple1000 = countOnesitefeature(onesiteinfo_simple_1000, refBase)
                    reads_ref_simple100Wei, reads_alt_simple100Wei, reads_num_simple100Wei = countOnesitefeature(onesiteinfo_simpleWei_100, refBase)
                    reads_ref_simple1000Wei, reads_alt_simple1000Wei, reads_num_simple1000Wei = countOnesitefeature(onesiteinfo_simpleWei_1000, refBase)
                    reads_ref_feature100, reads_alt_feature100, reads_num_feature100 = countOnesitefeature(onesiteinfo_feature_100, refBase)
                    reads_ref_feature1000, reads_alt_feature1000, reads_num_feature1000 = countOnesitefeature(onesiteinfo_feature_1000, refBase)

                    # Cell number
                    # original
                    altCellNum_original, totalCellsNum_original = count_mut_cellNum(baseQualityInfoOnesite, refBase)
                    # simple resample
                    altCellNum_simple100, totalCellsNum_simple100 = count_mut_cellNum(onesiteinfo_simple_100, refBase)
                    altCellNum_simple1000, totalCellsNum_simple1000 = count_mut_cellNum(onesiteinfo_simple_1000, refBase)
                    # simple resample Wei
                    altCellNum_simple100Wei, totalCellsNum_simple100Wei = count_mut_cellNum(onesiteinfo_simpleWei_100, refBase)
                    altCellNum_simple1000Wei, totalCellsNum_simple1000Wei = count_mut_cellNum(onesiteinfo_simpleWei_1000, refBase)
                    # feature resample
                    altCellNum_feature100, totalCellsNum_feature100 = count_mut_cellNum(onesiteinfo_feature_100, refBase)
                    altCellNum_feature1000, totalCellsNum_feature1000 = count_mut_cellNum(onesiteinfo_feature_1000, refBase)

                    sql = "\t".join([base_quality_original, base_quality_simple100, base_quality_simple1000, base_quality_simple100Wei, base_quality_simple1000Wei, base_quality_feature100, base_quality_feature1000, reads_ref_original, reads_alt_original, reads_num_original, reads_ref_simple100, reads_alt_simple100, reads_num_simple100, reads_ref_simple1000, reads_alt_simple1000, reads_num_simple1000, reads_ref_simple100Wei, reads_alt_simple100Wei, reads_num_simple100Wei, reads_ref_simple1000Wei, reads_alt_simple1000Wei, reads_num_simple1000Wei, reads_ref_feature100, reads_alt_feature100, reads_num_feature100, reads_ref_feature1000, reads_alt_feature1000, reads_num_feature1000, str(len(baseQualityInfoOnesite.split(";"))),
                                     altCellNum_original, totalCellsNum_original, altCellNum_simple100, totalCellsNum_simple100, altCellNum_simple1000, totalCellsNum_simple1000, altCellNum_simple100Wei, totalCellsNum_simple100Wei, altCellNum_simple1000Wei, totalCellsNum_simple1000Wei, altCellNum_feature100, totalCellsNum_feature100, altCellNum_feature1000, totalCellsNum_feature1000])
                    #sql = "\t".join([base_quality_original, base_quality_simple100, base_quality_simple1000, base_quality_simple100Wei, base_quality_simple1000Wei, base_quality_feature100, base_quality_feature1000, reads_ref_original, reads_alt_original, reads_num_original, reads_ref_simple100, reads_alt_simple100, reads_num_simple100, reads_ref_simple1000, reads_alt_simple1000, reads_num_simple1000, reads_ref_simple100Wei, reads_alt_simple100Wei, reads_num_simple100Wei, reads_ref_simple1000Wei, reads_alt_simple1000Wei, reads_num_simple1000Wei, reads_ref_feature100, reads_alt_feature100, reads_num_feature100, reads_ref_feature1000, reads_alt_feature1000, reads_num_feature1000, str(len(baseQualityInfoOnesite.split(";")))])
                    outfile.writelines(sql + "\n")

if __name__=='__main__':
    start = timeit.default_timer()
    #19_1_59128983
    main()
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')