from read_fasta_file_v2 import readfasta, decoy_pro
from itertools import islice
from protein_digestor_v2 import ProteinDigestion
from class_FDR import class_fdr

infasta = "Mavium_hominissuis_OCU464_GCF_001865635.4_ASM186563v4_protein_variant_DB_082021.fasta"

target_file = "M_avium_Mavium_hominissuis_variant_proteome_search_082021_TargetPsms.txt"
decoy_file = "M_avium_Mavium_hominissuis_variant_proteome_search_082021_DecoyPsms.txt"

def get_header_index(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')

            pro = split_i.index('ParentProteinAccessions')

            pep = split_i.index('ModifiedSequence')

            mod = split_i.index('Modifications')

            xcorr = split_i.index('XCorr')

            fdr = split_i.index('PercolatorqValue')

            return pro, pep, mod, xcorr, fdr


def digest_fasta(infasta, peps_type):
    dicts = {}
    for rows in readfasta(infasta).read():
        header = rows[0]

        if peps_type == 'target':
            accession = 'SAAV@' + header.split(' ')[0].split('|')[0]
            seq = rows[1]
            for miss_cleave in range(0, 3):
                peps = ProteinDigestion(seq, miss_cleave, 6, 80).trypsin()
                for pep in peps:
                    if pep not in dicts:
                        dicts[pep] = [accession]
                    else:
                        dicts[pep].append(accession)

        elif peps_type == 'decoy':
            decoy_accession = 'XXX_SAAV@' + header.split(' ')[0].split('|')[0]
            decoy_seq = decoy_pro(rows[1]).Reverse()
        
            for miss_cleave in range(0, 3):
                peps = ProteinDigestion(decoy_seq, miss_cleave, 6, 80).trypsin()
                for pep in peps:
                    if pep not in dicts:
                        dicts[pep] = [decoy_accession]
                    else:
                        dicts[pep].append(decoy_accession)
                    
    return dicts
   
def get_file_info(infile, pep_type):
    a = get_header_index(infile)
    if pep_type == 'target':
        peps = digest_fasta(infasta, pep_type)
        with open(infile) as file:
            for i in islice(file, 1, None):
                split_i = i.rstrip().split('\t')
                pep = split_i[a[1]]
                pro = split_i[a[0]]
                mod = split_i[a[2]]
                xcorr = split_i[a[3]]
                q_val = split_i[a[4]]
                if len(xcorr) != 0:
                    if pep in peps:
                        pro = ';'.join(peps[pep])
                        
                        yield pep + '_' + mod, pro, xcorr, q_val

    elif pep_type == 'decoy':
        peps = digest_fasta(infasta, pep_type)
        with open(infile) as file:
            for i in islice(file, 1, None):
                split_i = i.rstrip().split('\t')
                pep = split_i[a[1]]
                pro = split_i[a[0]]
                mod = split_i[a[2]]
                xcorr = split_i[a[3]]
                q_val = split_i[a[4]]
                if len(xcorr) != 0:
                    if pep in peps:
                        pro = ';'.join(peps[pep])
                        
                        yield pep + '_' + mod, pro, xcorr, q_val

output = []
for info in get_file_info(target_file, 'target'):
    output.append(info)
    
for info in get_file_info(decoy_file, 'decoy'):
    output.append(info)

for result in class_fdr(output):
    print (result)


    
