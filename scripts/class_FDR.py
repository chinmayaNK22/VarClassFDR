#### This section of this script is direct adaptation of BayesClassSpecificFDR.py from "https://github.com/yafeng/proteogenomics_python". However, it has been modified for Xcorr from Proteome Discoverer

import sys
import os
import getopt
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt

def plot_temp_scatter(decoy_dic, search_score):
    x = []
    y = []
    for score in decoy_dic:
        x.append(score)
        frac=float(decoy_dic[score][0]/decoy_dic[score][1])
        y.append(frac)
    fig=plt.figure()
    ax = fig.add_subplot(111) 
    ax.scatter(x,y,label = "Real data",s=1)
    ax.legend()
    plt.xlabel(search_score)
    plt.ylabel('Proportion (Variant decoy/All decoys)')
    plt.show()

def save_scatter_plot(x, y, coefs, search_score, low_threshold, upper_threshold, inpath):
    fit = poly.polyval(x, coefs)

    fig=plt.figure()
    ax = fig.add_subplot(111) 
    ax.scatter(x,y,label = "Real data",s=1)
    ax.plot(x,fit,label = "Polynomial with order=1", color='C1')
    ax.legend()
    plt.xlabel(search_score)
    plt.ylabel('Proportion (Variant decoy/All decoys)')
    fig_name = os.path.join(inpath , "Least_square_regression_curve_" + str(low_threshold) + '_' + str(upper_threshold) + '_' +search_score  + ".pdf")
    fig.savefig(fig_name, bbox_inches="tight", dpi=300)

def calc_least_square_reg(decoy_dic, search_score, inpath):
    ## plots a scatter plot for visualization of score distribution.
    ## Based on the consistency of the score distribution, users have to provide the upper and lower score thresholds for the calculation of least square regression
    
    plot_temp_scatter(decoy_dic, search_score)
    x=[]
    y=[]
    x_filter = []
    y_filter = []
    low_threshold = input('Provide a lower score cutoff point for the calculation of coefficients:')
    upper_threshold = input('Provide a upper score cutoff for the calculation of coefficients:')
    for score in decoy_dic:
        x.append(score)
        frac=float(decoy_dic[score][0]/decoy_dic[score][1])
        y.append(frac)
        if float(low_threshold)<score<float(upper_threshold):
            x_filter.append(score)
            y_filter.append(frac)
            
    coefs = poly.polyfit(np.array(x_filter),np.array(y_filter), 1)
    print("coefficients",coefs)
    
    ## Plots a scattered plot with a slope based on least-square regression calculation for the povided score thresholds
    save_scatter_plot(x, y, coefs, search_score, low_threshold, upper_threshold, inpath)

    return coefs


def _classfdr(inpath, psm_info, search_score, fdr_cutoff, percolator):

    psm_qval = 0
    if fdr_cutoff != None or fdr_cutoff != 0:
        psm_qval = float(fdr_cutoff)
    else:
        psm_qval = 0.01
    
    score_dic = {}
    decoy_dic={}

    novel_targetcount=0
    novel_decoycount=0

    decoycount=0
    targetcount=0

    novpep_dic={}

    output = []
        
    decoy_prefix = "XXX_SAAV@"
    novel_prefix = "SAAV@"
    for line in psm_info:
        pro=line[1]
        score_val = float(line[-2])
    
        if "XXX_" in pro:
            decoycount+=1
            if decoy_prefix in pro:
                novel_decoycount+=1
            decoy_dic[score_val] = [novel_decoycount,decoycount]
        else:
            targetcount+=1
            if novel_prefix in pro:
                novel_targetcount+=1
 
        score_dic[score_val] = [targetcount,decoycount,novel_targetcount,novel_decoycount]

    ## Calculate the regression coefficients from the least-square regression analysis
    coefs = calc_least_square_reg(decoy_dic, search_score, inpath)
    
    new_noveltargetcount = 0
    new_noveldecoycount = 0
    new_score_dic = {}
    int_output = []
    for line in psm_info:
        pep=line[0]
        pro=line[1]
        if novel_prefix not in pro:
            continue;
        
        score_val = float(line[-2])
        counts = score_dic[score_val]   

        target_n=float(counts[0])
        decoy_n=float(counts[1])

        if percolator.upper() == 'PERCOLATOR':
            FDR=float(line[-1])
        else:
            FDR=decoy_n/target_n
        
        novel_targetcount=float(counts[2])
        gamma = poly.polyval(score_val, coefs)
        gamma = (coefs[1] * score_val) + coefs[0]
        if novel_targetcount != 0:
            classFDR = FDR*gamma*(target_n/novel_targetcount)

            if pep not in novpep_dic:
                novpep_dic[pep] = classFDR
        
            if "XXX_SAAV@" not in pro: #write only target PSMs
                int_output.append([pep, pro] + line[2:5] + [score_val, line[-1], FDR, classFDR, novpep_dic[pep]])
                if classFDR < psm_qval:
                    new_noveltargetcount +=1
                elif classFDR > psm_qval:
                    new_noveldecoycount +=1

                new_score_dic[score_val] = [new_noveltargetcount,new_noveldecoycount]

                #output.append([pep.split('_')[0], pep.split('_')[1], pro, str(xcorr), str(line[3]), str(classFDR), str(novpep_dic[pep])])
            
    for line in int_output:
        pep=line[0]
        pro=line[1]
        score_val = float(line[5])
        counts = new_score_dic[score_val]   

        new_target_n=float(counts[0])
        new_decoy_n=float(counts[1])
        if new_target_n != 0:
            FDP = float(new_decoy_n/new_target_n)

            classFDR_error = float(line[-2]) - float(FDP)
            output.append([pep.split('_')[0], pep.split('_')[1], pro] + line[2:5] + [str(score_val), str(line[6]), str(line[7]), str(line[8]),str(FDP), str(classFDR_error)])

    print (f"Hits in novel search space: targe,decoy", novel_targetcount, novel_decoycount)

    return output
