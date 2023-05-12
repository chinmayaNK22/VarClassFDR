from itertools import islice
import os
from scripts.read_fasta_file_v2 import readfasta, decoy_pro
from scripts.protein_digestor_v2 import ProteinDigestion
from scripts import pdResult_parser
from scripts import class_FDR
import argparse

parser = argparse.ArgumentParser(description='''Calculate False Discovery Rate (FDR) for variant peptides identified from Proteome Discoverer tool using various search engines (SequestHT//Mascot//MSAmanda)''')

parser.add_argument('infile', metavar='-ip', type=str, nargs='+', help='Proteome Discoverer output in msf/pdResult format')
parser.add_argument('fastas', metavar='-vf', type=str, nargs='+', help='Path to a reference proteome FASTA file. Provide path to a folder if more than single FASTA files were used for the search apart from variant FASTA file')
parser.add_argument('var_fasta', metavar='-rf', type=str, nargs='+', help='Variant protein sequences used for the search in FASTA format')
parser.add_argument('searchEngine', metavar='-SE', type=str, nargs='+', help='''Currently, database search output from SequestHT (1), Mascot (2)and MSAmanda (3)can be used.
                                                                                Each search engines is annotated with a number, which has to be provide when defining the type of search engine used.''')
parser.add_argument('q-value_cutoff', metavar='-fdr', type=str, default=0.01, nargs='+', help='Set the Class-FDR (q-value) cutoff to be applied. Ex: 1%% class-FDR will be 0.01 and 5%% will be 0.05')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args = parser.parse_args()

algorithm_score_type = {'Sequest HT':'XCorr', 'MS Amanda 2.0':'AmandaScore', 'Mascot':'IonScore'}

def header_index(inlist, column_name):
    return inlist.index(column_name)


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

   
def get_file_info(infile, fastas_path, var_fasta, search_eng):
    peps = read_fasta(fastas_path, var_fasta)
    mut_peps = digest_fasta(var_fasta)
    t_header, t_table = pdResult_parser.fetch_sqlite3_db(infile, 'target')
    d_header, d_table = pdResult_parser.fetch_sqlite3_db(infile, 'decoy')
    
    search_engines = {i[header_index(t_header, 'IdentifyingNodeTypeName')]:1 for i in t_table}
    print (f'The search was performed using ', ' and '.join(list(search_engines.keys())), 'alogorithms.')

    for k in search_engines.keys():
        if search_eng == k:
            spectra_hits = []
            for i in t_table:
                pep = i[header_index(t_header,'Sequence')] 
                mod_pep = i[header_index(t_header,'ModifiedSequence')] 
                pro = i[header_index(t_header,'ParentProteinAccessions')]
                mod = i[header_index(t_header,'Modifications')] 
                score = i[header_index(t_header, algorithm_score_type[k])]
                q_val = i[header_index(t_header,'PercolatorqValue')]
                rawfile = i[header_index(t_header, 'SpectrumFileName')]
                scan = str(i[header_index(t_header, 'FirstScan')])
                rt = str(round(float(i[header_index(t_header,'RetentionTime')]),4))
                #pep_id = i[header_index(t_header,'PeptideID')]
                if score != None:
                    pro = ""
                    if pep.upper() in peps:
                        pro = ';'.join(peps[pep.upper()])
                    elif pep.upper() in mut_peps:
                        pro = ';'.join('SAAV@' + p for p in mut_peps[pep.upper()])

                    spectra_hits.append([mod_pep + '_' + mod, pro, rt, scan, rawfile, score, q_val])
            
            for d_i in d_table:
                pep = d_i[header_index(d_header,'Sequence')]
                mod_pep = d_i[header_index(d_header,'ModifiedSequence')] 
                pro = d_i[header_index(d_header,'ParentProteinAccessions')]
                mod = d_i[header_index(d_header,'Modifications')] 
                score = d_i[header_index(d_header, algorithm_score_type[k])]
                q_val = d_i[header_index(d_header,'PercolatorqValue')]
                rawfile = d_i[header_index(d_header, 'SpectrumFileName')]
                
                if score != None:
                    if 'decoy_' + pep.upper() in peps:
                        pro = ';'.join(peps['decoy_' + pep.upper()])
                    elif 'decoy_' + pep.upper() in mut_peps:
                        pro = ';'.join('XXX_SAAV@'+ p.lstrip('XXX_') for p in mut_peps['decoy_' + pep.upper()])

                    spectra_hits.append([mod_pep + '_' + mod, pro, '','', rawfile, score, q_val])

            return k, spectra_hits

def run_classfdr(infile, fastas_path, var_fasta, search_id, fdr_cutoff, percolator):
    searchEngines = {'1':'Sequest HT', '2':'Mascot', '3':'MS Amanda 2.0'}

    searchEngine, output = get_file_info(infile, fastas_path, var_fasta, searchEngines[str(search_id)])
    
    print (f'There were ', len(output), ' PSMs stored from ', searchEngine, ' search algorithm hits.')
    
    if len(output) != 0:

        scores = {idx:float(o[-2]) for idx, o in enumerate(output)}

        sort_output = [output[k] for k, v in sorted(scores.items(), key=lambda item: item[1], reverse=True)]

        ### Perform class-FDR calculation
        inpath = os.path.split(infile)[0]
        
        result = class_FDR._classfdr(inpath, sort_output, algorithm_score_type[searchEngine], fdr_cutoff, percolator)

        ### Write the result
        
        if percolator.upper() == "PERCOLATOR":
            outfile = infile.rstrip(infile.split('.')[-1]).rstrip('.') + "_" + searchEngines[str(search_id)] + '_Percolator_classFDR.txt'
        else:
            outfile = infile.rstrip(infile.split('.')[-1]).rstrip('.') + "_" + searchEngines[str(search_id)] + '_classFDR.txt'
            
        print (f'The class-FDR calculation was performed for ', len(result), 'PSMs and stored in', outfile)
        
        with open(outfile, 'w') as outf:
           outf.write('Peptide\tModification\tProtein\tRT (min)\tScan\tRaw File\tXcorr\tPercolator qvalue\tGlobal-FDR\tclass-FDR\tFalse Discovery Proportion (FDP)\tclass-FDR error\n')
           outf.writelines('\t'.join(i) + '\n' for i in result)
    else:
        raise Exception
    
if __name__ == '__main__':
    run_classfdr(args.infile[0], args.fastas[0], args.var_fasta[0], args.searchEngine[0], args.q-value_cutoff[0], '')
