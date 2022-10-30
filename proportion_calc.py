from itertools import islice
import matplotlib.pyplot as plt
import numpy as np

decoy_file = "M_avium_Mavium_hominissuis_variant_proteome_search_082021_DecoyPsms.txt"

target_file = "M_avium_Mavium_hominissuis_variant_proteome_search_082021_TargetPsms.txt"


def get_count(infile):
    xcorrs = {}
    c = 0
    #arr = np.array([])
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            xcorr = split_i[-4]
            if len(xcorr) != 0:
                Xcorr = round(float(split_i[-4]), 4)
                if Xcorr not in xcorrs:
                    xcorrs[Xcorr] = [Xcorr]
                else:
                    xcorrs[Xcorr].append(Xcorr)
                    
                if Xcorr >= 2.9:
                    #arr = np.append(arr, xcorr)
                    c+= 1 
                    #xcorrs[xcorr] = 1

    proportions = {k: k/c for k, v in xcorrs.items()}
    for i in sorted(proportions):
        print((i, proportions[i]), end=" ")
        
    return c
##    fig, axs = plt.subplots(1, 1,
##                            figsize =(10, 7),
##                            tight_layout = True)
##     
##    plt.hist(list(xcorrs))
##     
##    # Show plot
##    plt.show()

decoys = get_count(decoy_file)
targets = get_count(target_file)

FDR = decoys/targets

print (decoys, targets, FDR)
