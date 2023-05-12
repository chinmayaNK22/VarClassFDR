# VarClassFDR -Variant peptide class FDR calculator

Measurement of False Discovery Rate (FDR) specific to variant peptides (class-FDR) for is crucial their accurate identification and removes any false hits being considered as true hits based on Global-FDR calculation. Therefore, we developed a Python-based tool for the calculation of class-FDR for variant peptides identified using database search in Proteome Discoverer (Thermo Scientific).

## Command line usage of VarClassFDR.py
```
>python VarClassFDR.py input.msf variant_proteome.fasta reference_proteome.fasta 1 0.01
```

## Details of required inputs
```
usage: VarClassFDR.py [-h] [--version] -ip [-ip ...] -vf [-vf ...] -rf [-rf ...] -SE [-SE ...] -fdr [-fdr ...]

Calculate False Discovery Rate (FDR) for variant peptides identified from Proteome Discoverer tool using various search engines (SequestHT//Mascot//MSAmanda)

positional arguments:
  -ip         Proteome Discoverer output in msf/pdResult format
  -vf         Path to a reference proteome FASTA file. Provide path to a folder if more than single FASTA files were used for the search apart from variant FASTA file
  -rf         Variant protein sequences used for the search in FASTA format
  -SE         Currently, database search output from SequestHT (1), Mascot (2)and MSAmanda (3)can be used. Each search engines is annotated with a number, which has to be provide when defining the type of search engine used.
  -fdr        Set the Class-FDR (q-value) cutoff to be applied. Ex: 1% class-FDR will be 0.01 and 5% will be 0.05

options:
  -h, --help  show this help message and exit
  --version   show program's version number and exit
  ```
