import os
import re
import timeit
import argparse
import numpy as np
from scipy.stats import beta
from scipy.optimize import fmin
from py_utils import *

def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to obtain a list of base and cell counts in scRNA bam file')
    parser.add_argument('--outdir', type=str, default=1, help='outdir', required = True)
    parser.add_argument('--tmpvcfdir', type=str, default=1, help='Path to random select vcfcount path', required = True)
    parser.add_argument('--id', help='Prefix used to name output file', required = True)
    return (parser)

def merge_VCF_files(outPath, VCFtmpPath, IDSample):
    #OutDir = "/Users/liqun/BrainRNA/SNP/testData/scRNA/tmp/subam/"
    #VCFtmpDir = "/Users/liqun/BrainRNA/SNP/testData/scRNA/tmp/subam/tmpVCF/"
    #IDSample = "P6_cSCC"
    filesTotal = os.listdir(VCFtmpPath)
    VCFtmpTotal = [i for i in filesTotal if re.match("random.*", i) != None]
    VCFtmpTotalPaths = [VCFtmpPath + "/" + i for i in VCFtmpTotal]
    PathsSql = " ".join(VCFtmpTotalPaths)
    sql = "cat " + PathsSql + " > " + outPath + "VCF_" + IDSample + ".vcf"
    os.system(sql)

def betaNLL(param,*args):
    '''Negative log likelihood function for beta
    <param>: list for parameters to be fitted.
    <args>: 1-element array containing the sample data.

    Return <nll>: negative log-likelihood to be minimized.
    '''

    a,b=param
    data=args[0]
    pdf=beta.pdf(data,a,b,loc=0,scale=1)
    lg=np.log(pdf)
    #-----Replace -inf with 0s------
    lg=np.where(lg==-np.inf,0,lg)
    nll=-1*np.sum(lg)
    return nll

def fit_using_moment(data):
    mean=np.mean(data)
    var=np.var(data,ddof=1)
    alpha=mean**2*(1-mean)/var-mean
    beta=alpha*(1-mean)/mean
    return alpha, beta


def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outPath = args.outdir
    VCFtmpPath = args.tmpvcfdir
    IDSample = args.id
    merge_VCF_files(outPath, VCFtmpPath, IDSample)
    tmpoutrandom = outPath + "randomSelect_" + IDSample + ".tmp"
    random_extractsnp_fromVCF(outPath + "VCF_" + IDSample + ".vcf", tmpoutrandom, 1000000)
    # return success and trailArray
    # successesSingleCellReads, trialsSingleCellReads, successesCellNumber, trialsCellNumber, successesSiteTotalReads, trialsSiteTotalReads
    a, b, c, d, e, f = generate_success_trailArray(tmpoutrandom)
    data_altReadsFreOneCell = a/b
    data_altCellNumFreOneSite = c/d
    data_altReadsFreOneSite = e/f
    #----------------Fit using moments----------------
    alpha_data1, beta_data1 = fit_using_moment(data_altReadsFreOneCell)
    alpha_data2, beta_data2 = fit_using_moment(data_altCellNumFreOneSite)
    alpha_data3, beta_data3 = fit_using_moment(data_altReadsFreOneSite)
    '''
    #------------------Fit using mle------------------
    result1_data1=fmin(betaNLL,[1,1],args=(data1,))
    result1_data2=fmin(betaNLL,[1,1],args=(data2,))
    alpha2_data1,beta2_data1=result1_data1
    alpha2_data2,beta2_data2=result1_data2
    #----------------Fit using beta.fit----------------
    alpha3_data1,beta3_data1,xx_data1,yy_data1=beta.fit(data1)
    alpha3_data2,beta3_data2,xx_data2,yy_data2=beta.fit(data2)
    '''
    print("#############################################")
    print ('\n# alpha,beta from moments [altReadsFreOneCell]:', alpha_data1, beta_data1)
    print ('# alpha, beta from moments [altCellNumFreOneSite]:',alpha_data2, beta_data2)
    print ('# alpha,beta from moments [altReadsFreOneSite]:',alpha_data3, beta_data3)
'''
def update_from_data_by_moment(tries, success):
    #estimate alpha, beta using moment estimation
    mean, var = compute_moment(tries, success)
    #print 'mean and variance: ', mean, var
    #self.alpha = mean*(mean*(1-mean)/(var+0.000001)-1)
    alpha = (mean+0.000001) * ((mean+0.000001) * (1.000001 - mean) / (var+0.000001) - 1)
    #self.beta = (1-mean)*(mean*(1-mean)/(var+0.000001)-1)
    beta = (1.000001 - mean) * ((mean+0.000001) * (1.000001 - mean) / (var+0.000001) - 1)
    return alpha, beta

def compute_moment(tries, success):
    #moment estimation
    ctr_list = []
    var = 0.0
    for i in range(len(tries)):
        ctr_list.append(float(success[i])/tries[i])
    mean = sum(ctr_list)/len(ctr_list)
    for ctr in ctr_list:
        var += pow(ctr-mean, 2)
    return mean, var/(len(ctr_list)-1)

def main():
    parser = initialize_parser()
    args = parser.parse_args()
    outPath = args.outpath
    VCFtmpPath = args.tmpvcfpath
    IDSample = args.id
    merge_VCF_files(outPath, VCFtmpPath, IDSample)
    tmpoutrandom = outPath + "randomSelect_" + IDSample + ".tmp"
    random_extractsnp_fromVCF(outPath + "VCF_" + IDSample + ".vcf", tmpoutrandom, 100000)
    a, b, c, d = generate_success_trailArray(tmpoutrandom)
    print(a, b, c, d)
    ratio_ab = a/b
    ratio_cd = c/d
    print(beta.fit(ratio_ab))
    print(beta.fit(ratio_cd))
    l1, m1 = estimate_shape_parameters_mle(a, b)
    l2, m2 = update_from_data_by_moment(b, a)
    o1, p1 = estimate_shape_parameters_mle(c, d)
    o2, p2 = update_from_data_by_moment(d, c)
    return l1, m1, o1, p1, l2, m2, o2, p2
'''


if __name__=='__main__':
    start = timeit.default_timer()
    #19_1_59128983
    k = main()
    print(k)
    stop = timeit.default_timer()
    print('All processes done.')
    print ("Computation time: " + str(round(stop - start)) + ' seconds')