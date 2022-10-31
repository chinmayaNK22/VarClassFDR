from read_fasta_file_v2 import readfasta, decoy_pro
from itertools import islice
from protein_digestor_v2 import ProteinDigestion
from class_FDR import class_fdr
import os
from pdResult_parser import fetch_sqlite3_db

infile = "..\PRIDE\M_avium_Mavium_hominissuis_variant_proteome_search_082021.pdResult"

var_fasta = "Mavium_hominissuis_OCU464_GCF_001865635.4_ASM186563v4_protein_variant_DB_082021.fasta"

fastas = "."

def get_header_index(inlist):
    pro = inlist.index('ParentProteinAccessions')

    pep = inlist.index('Sequence')

    mod_pep = inlist.index('ModifiedSequence')

    mod = inlist.index('Modifications')

    xcorr = inlist.index('XCorr')

    fdr = inlist.index('PercolatorqValue')

    return pro, pep, mod_pep, mod, xcorr, fdr

def digest_fasta(infasta):
    dicts = {}
    for rows in readfasta(infasta).read():
        header = rows[0]
        accession = header.split(' ')[0].split('|')[0]
        seq = rows[1]
        for miss_cleave in range(0, 3):
            peps = ProteinDigestion(seq, miss_cleave, 6, 80).trypsin()
            for pep in peps:
                if pep not in dicts:
                    dicts[pep] = [accession]
                else:
                    dicts[pep].append(accession)

        decoy_accession = 'XXX_' + header.split(' ')[0].split('|')[0]
        decoy_seq = decoy_pro(rows[1]).Reverse()
        for miss_cleave in range(0, 3):
            peps = ProteinDigestion(decoy_seq, miss_cleave, 6, 80).trypsin()
            for pep in peps:
                decoy_pep = 'decoy_' + pep
                if decoy_pep not in dicts:
                    dicts[decoy_pep] = [decoy_accession]
                else:
                    dicts[decoy_pep].append(decoy_accession)
                    
    return dicts

def read_fasta(inpath, var_fasta):
    extension = {'fasta':1, 'faa':2, 'fa':3}
    digested_peps = {}
    for files in os.listdir(os.path.join(inpath)):
        if os.path.isfile(os.path.join(inpath, files)):
            if files.split('.')[-1] in extension:
                if files != os.path.split(var_fasta)[-1]:
                    digested_peps |= digest_fasta(files)
                
    return digested_peps
   
def get_file_info(infile, pep_type):
    peps = read_fasta(fastas, var_fasta)
    mut_peps = digest_fasta(var_fasta)
    if pep_type == 'target':
        header, table = fetch_sqlite3_db(infile, pep_type)
        a = get_header_index(header)
        for i in table:
            pep = i[a[1]]
            mod_pep = i[a[2]]
            pro = i[a[0]]
            mod = i[a[3]]
            xcorr = i[a[4]]
            q_val = i[a[5]]
            if xcorr != None:
                pro = ""
                if pep in peps:
                    pro = ';'.join(peps[pep])
                elif pep in mut_peps:
                    pro = ';'.join('SAAV@' + p for p in mut_peps[pep])

                yield mod_pep + '_' + mod, pro, xcorr, q_val

    elif pep_type == 'decoy':
        header, table = fetch_sqlite3_db(infile, pep_type)
        a = get_header_index(header)
        for i in table:
            pep = i[a[1]]
            mod_pep = i[a[2]]
            pro = i[a[0]]
            mod = i[a[3]]
            xcorr = i[a[4]]
            q_val = i[a[5]]
            if xcorr != None:
                if 'decoy_' + pep in peps:
                    pro = ';'.join(peps['decoy_' + pep])
                elif 'decoy_' + pep in mut_peps:
                    pro = ';'.join('XXX_SAAV@'+ p.lstrip('XXX_') for p in mut_peps['decoy_' + pep])
                    
                yield mod_pep + '_' + mod, pro, xcorr, q_val

output = []
for info in get_file_info(infile, 'target'):
    output.append(info)
    
for info in get_file_info(infile, 'decoy'):
    output.append(info)

scores = {idx:float(o[2]) for idx, o in enumerate(output)}

sort_output = [output[k] for k, v in sorted(scores.items(), key=lambda item: item[1], reverse=True)]

result = class_fdr(sort_output)

outfile = "{0}_classFDR_percolatorFDR.txt".format(infile.rstrip(infile.split('.')[-1]).rstrip('.'))
with open(outfile, 'w') as outf:
   outf.write('Peptide\tModification\tProtein\tXcorr\tPercolator_qvalue\tclass_FDR\tFDP\tclassFDR_error\n')
   outf.writelines('\t'.join(i) + '\n' for i in result)


    
