import sys
import os
import getopt
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt



def class_fdr(psm_info):
    score_dic = {}
    decoy_dic={}

    novel_targetcount=0
    novel_decoycount=0

    decoycount=0
    targetcount=0

    novpep_dic={}

    output = []

    psm_qval = 0.01
    decoy_prefix = "XXX_SAAV@"
    novel_prefix = "SAAV@"
    for line in psm_info:
        pro=line[1]
        xcorr = float(line[2])
    
        if "XXX_" in pro:
            decoycount+=1
            if decoy_prefix in pro:
                novel_decoycount+=1
            decoy_dic[xcorr] = [novel_decoycount,decoycount]
        else:
            targetcount+=1
            if novel_prefix in pro:
                novel_targetcount+=1
 
        score_dic[xcorr] = [targetcount,decoycount,novel_targetcount,novel_decoycount]

    x=[]
    y=[]
    x_filter = []
    y_filter = []

    for score in decoy_dic:
        x.append(score)
        frac=float(decoy_dic[score][0]/decoy_dic[score][1])
        y.append(frac)
        if 0<score<2:
            x_filter.append(score)
            y_filter.append(frac)
        
    coefs = poly.polyfit(np.array(x_filter),np.array(y_filter), 1)
    print("coefficients",coefs)
    fit = poly.polyval(x, coefs)

    fig=plt.figure()
    ax = fig.add_subplot(111) 
    ax.scatter(x,y,label = "Real data",s=1)
    ax.plot(x,fit,label = "Polynomial with order=1", color='C1')
    ax.legend()
    plt.xlabel('xcorr')
    plt.ylabel('Proportion (Variant decoy/All decoys)')
    fig.savefig("fitcurve.png")
    
    new_noveltargetcount = 0
    new_noveldecoycount = 0
    new_score_dic = {}
    int_output = []
    for line in psm_info:
        pep=line[0]
        pro=line[1]
        if novel_prefix not in pro:
            continue;
        
        xcorr = float(line[2])
        counts = score_dic[xcorr]   

        target_n=float(counts[0])
        decoy_n=float(counts[1])
        #FDR=decoy_n/target_n
        FDR = float(line[3])

        novel_targetcount=float(counts[2])
        gamma = poly.polyval(xcorr, coefs)
        gamma = (coefs[1] * xcorr) + coefs[0]
        classFDR = FDR*gamma*(target_n/novel_targetcount)

        if pep not in novpep_dic:
            novpep_dic[pep] = classFDR
    
        if "XXX_SAAV@" not in pro: #write only target PSMs
            int_output.append([pep, pro, xcorr, line[3], classFDR, novpep_dic[pep]])
            if classFDR < psm_qval:
                new_noveltargetcount +=1
            elif classFDR > psm_qval:
                new_noveldecoycount +=1

            new_score_dic[xcorr] = [new_noveltargetcount,new_noveldecoycount]

            #output.append([pep.split('_')[0], pep.split('_')[1], pro, str(xcorr), str(line[3]), str(classFDR), str(novpep_dic[pep])])
            
    for line in int_output:
        pep=line[0]
        pro=line[1]
        xcorr = float(line[2])
        counts = new_score_dic[xcorr]   

        new_target_n=float(counts[0])
        new_decoy_n=float(counts[1])
        FDP=new_decoy_n/new_target_n

        classFDR_error = float(line[4]) - float(FDP)
        output.append([pep.split('_')[0], pep.split('_')[1], pro, str(xcorr), str(line[3]), str(line[4]), str(FDP), str(classFDR_error)])

    print ("Hits in novel search space: targe,decoy",novel_targetcount,novel_decoycount)


    return output
