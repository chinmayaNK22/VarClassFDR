from read_fasta_file_v2 import readfasta, decoy_pro
from itertools import islice
from protein_digestor_v2 import ProteinDigestion
from class_FDR import class_fdr
import os
from pdResult_parser import fetch_sqlite3_db

infile = "F:\\NTNU_M_avium_Proteomics\\PRIDE\\M_avium_Mavium_hominissuis_variant_proteome_search_082021.pdResult"

var_fasta = "F:\\NTNU_M_avium_Proteomics\\classFDR\\Mavium_hominissuis_OCU464_GCF_001865635.4_ASM186563v4_protein_variant_DB_082021.fasta"

fastas = "F:\\NTNU_M_avium_Proteomics\\classFDR"

def get_header_index(inlist):
    pro = inlist.index('ParentProteinAccessions')

    pep = inlist.index('Sequence')

    mod_pep = inlist.index('ModifiedSequence')

    mod = inlist.index('Modifications')

    xcorr = inlist.index('XCorr')

    try:
        msamandascore = inlist.index('AmandaScore')
    except:
        raise ("The search result does not contain AmandaScore")

    fdr = inlist.index('PercolatorqValue')

    return pro, pep, mod_pep, mod, xcorr, fdr, msamandascore

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
                    digested_peps |= digest_fasta(os.path.join(inpath, files))
                
    return digested_peps
   
def get_file_info(infile, pep_type, search_engine, fastas_path, var_fasta):
    peps = read_fasta(fastas_path, var_fasta)
    mut_peps = digest_fasta(var_fasta)
    if pep_type == 'target':
        header, table = fetch_sqlite3_db(infile, pep_type)
        a = get_header_index(header)
        for i in table:
            pep = i[a[1]]
            mod_pep = i[a[2]]
            pro = i[a[0]]
            mod = i[a[3]]
            score = ""
            if search_engine == 'sequest':
                score = i[a[4]]
            elif search_engine == 'amanda':
                score = i[a[-1]]
            q_val = i[a[5]]
            if score != None:
                pro = ""
                if pep.upper() in peps:
                    pro = ';'.join(peps[pep.upper()])
                elif pep.upper() in mut_peps:
                    pro = ';'.join('SAAV@' + p for p in mut_peps[pep.upper()])

                yield mod_pep + '_' + mod, pro, score, q_val

    elif pep_type == 'decoy':
        header, table = fetch_sqlite3_db(infile, pep_type)
        a = get_header_index(header)
        for i in table:
            pep = i[a[1]]
            mod_pep = i[a[2]]
            pro = i[a[0]]
            mod = i[a[3]]
            score = ""
            if search_engine == 'sequest':
                score = i[a[4]]
            elif search_engine == 'amanda':
                score = i[a[-1]]
            q_val = i[a[5]]
            if score != None:
                if 'decoy_' + pep.upper() in peps:
                    pro = ';'.join(peps['decoy_' + pep.upper()])
                elif 'decoy_' + pep.upper() in mut_peps:
                    pro = ';'.join('XXX_SAAV@'+ p.lstrip('XXX_') for p in mut_peps['decoy_' + pep.upper()])
                    
                yield mod_pep + '_' + mod, pro, score, q_val

def run_varclassfdr(infile, search_engine, fastas_path, var_fasta):
    output = []
    for info in get_file_info(infile, 'target', search_engine, fastas_path, var_fasta):
        output.append(info)
    
    for info in get_file_info(infile, 'decoy', search_engine, fastas_path, var_fasta):
        output.append(info)

    scores = {idx:float(o[2]) for idx, o in enumerate(output)}

    sort_output = [output[k] for k, v in sorted(scores.items(), key=lambda item: item[1], reverse=True)]

    inpath = os.path.split(infile)[0]
    
    result = class_fdr(inpath, sort_output, search_engine)

    outfile = "{0}_SequestHT_classFDR.txt".format(infile.rstrip(infile.split('.')[-1]).rstrip('.'))
    with open(outfile, 'w') as outf:
       outf.write('Peptide\tModification\tProtein\tXcorr\tPercolator qvalue\tclass-FDR\tFalse Discovery Proportion (FDP)\tclass-FDR error\n')
       outf.writelines('\t'.join(i) + '\n' for i in result)


run_varclassfdr(infile, 'amanda', fastas, var_fasta)    
