from read_fasta_file_v2 import readfasta, decoy_pro
from itertools import islice
from protein_digestor_v2 import ProteinDigestion
from class_FDR import class_fdr
import os

var_fasta = "Mavium_hominissuis_OCU464_GCF_001865635.4_ASM186563v4_protein_variant_DB_082021.fasta"

target_file = "M_avium_Mavium_hominissuis_variant_proteome_search_082021_TargetPsms.txt"
decoy_file = "M_avium_Mavium_hominissuis_variant_proteome_search_082021_DecoyPsms.txt"

fastas = "."

def get_header_index(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')

            pro = split_i.index('ParentProteinAccessions')

            pep = split_i.index('Sequence')

            mod_pep = split_i.index('ModifiedSequence')

            mod = split_i.index('Modifications')

            xcorr = split_i.index('XCorr')

            fdr = split_i.index('PercolatorqValue')

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
    a = get_header_index(infile)
    mut_peps = digest_fasta(var_fasta)
    if pep_type == 'target':
        with open(infile) as file:
            for i in islice(file, 1, None):
                split_i = i.rstrip().split('\t')
                pep = split_i[a[1]]
                mod_pep = split_i[a[2]]
                pro = split_i[a[0]]
                mod = split_i[a[3]]
                xcorr = split_i[a[4]]
                q_val = split_i[a[5]]
                if len(xcorr) != 0:
                    pro = ""
                    if pep in peps:
                        pro = ';'.join(peps[pep])
                    elif pep in mut_peps:
                        pro = ';'.join('SAAV@' + p for p in mut_peps[pep])

                    yield mod_pep + '_' + mod, pro, xcorr, q_val

    elif pep_type == 'decoy':
        with open(infile) as file:
            for i in islice(file, 1, None):
                split_i = i.rstrip().split('\t')
                pep = split_i[a[1]]
                mod_pep = split_i[a[2]]
                pro = split_i[a[0]]
                mod = split_i[a[3]]
                xcorr = split_i[a[4]]
                q_val = split_i[a[5]]
                if len(xcorr) != 0:
                    if 'decoy_' + pep in peps:
                        pro = ';'.join(peps['decoy_' + pep])
                    elif 'decoy_' + pep in mut_peps:
                        pro = ';'.join('XXX_SAAV@'+ p.lstrip('XXX_') for p in mut_peps['decoy_' + pep])
                        
                    yield mod_pep + '_' + mod, pro, xcorr, q_val

output = []
for info in get_file_info(target_file, 'target'):
    output.append(info)
    
for info in get_file_info(decoy_file, 'decoy'):
    output.append(info)

result = class_fdr(output)
outfile = "{0}_classFDR.txt".format(target_file.rstrip('txt').rstrip('.'))
with open(outfile, 'w') as outf:
    outf.write('Peptide\tModification\tProtein\tXcorr\tPercolator_qvalue\tclass_FDR\n')
    outf.writelines('\t'.join(i) + '\n' for i in output)


    
