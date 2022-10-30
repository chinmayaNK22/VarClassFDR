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
    #pep_col=header.index("Peptide")
    #prot_col=header.index("Protein")

    #specEval_col=header.index("SpecEValue")
    psm_qval = 0.01
    decoy_prefix = "XXX_SAAV@"
    novel_prefix = "SAAV@"
    for line in psm_info:
        #row=line.strip().split('\t')
        pro=line[1]
        #specEval = -np.log10(float(row[specEval_col]))
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
        if 0.5<score<10:
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
    plt.ylabel('Gamma (class specific decoy hits / total decoy hits)')
    fig.savefig("fitcurve.png")

    #input2=open(input_file,'r')
    
    for line in psm_info:
        #row=line.strip().split('\t')
        pep=line[0]
        pro=line[1]
        if novel_prefix not in pro:
            continue;
        
        #specEval = -np.log10(float(row[specEval_col]))
        xcorr = float(line[2])
        counts = score_dic[xcorr]   
    
        targetcount=float(counts[0])
        decoycount=float(counts[1])
        FDR=decoycount/targetcount
    
        novel_targetcount=float(counts[2])
        gamma = poly.polyval(xcorr, coefs)
        novelFDR = FDR*gamma*(targetcount/novel_targetcount)
    
        if pep not in novpep_dic:
            novpep_dic[pep]= novelFDR
        
        if novelFDR < psm_qval:
            if "XXX_SAAV@" not in pro: #write only target PSMs
                output.append([pep.split('_')[0], pep.split('_')[1], pro, xcorr, line[3], str(novelFDR), str(novpep_dic[pep])])
            #line.append(str(novelFDR))
            #line.append(str(novpep_dic[pep]))
       
            #output.write("\t".join(row)+"\n")
    
    print ("Hits in novel search space: targe,decoy",novel_targetcount,novel_decoycount)


    return output
