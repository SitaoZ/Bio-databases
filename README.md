# Bio-databases
Highly used databases in bioinformatics

## Reference database

To download reference data, there are some different sources available:
General biological databases: Ensembl, NCBI, and UCSC
Organism-specific biological databases: Wormbase, Flybase, TAIR etc. (often updated more frequently, so may be more comprehensive)
Reference data collections: Illumina's iGenomes, one location to access genome reference data from Ensembl, UCSC and NCBI

1. Ensembl
```bash 
Ensembl identifiers
When using Ensembl, note that it uses the following format for biological identifiers:
bold=$(tput bold)
normal=$(tput sgr0)
ENS${bold}G${normal}###########: Ensembl Gene ID
ENS**T**###########: Ensembl Transcript ID
ENS**P**###########: Ensembl Peptide ID
ENS**E**###########: Ensembl Exon ID
For non-human species a suffix is added:

ENSMUSG###: MUS (Mus musculus) for mouse
ENSDARG###: DAR (Danio rerio) for zebrafish

```
2. NCBI
https://www.ncbi.nlm.nih.gov/home/download/

  2.1 GeneBank   
  
**GCA** for GenBank assemblies   
**GCF** for RefSeq assemblies   
GCA(or GCF) is followed by an underscore and 9 digits. GRCh38.p11 is GCA_000001405.26.      

  2.2 BioProject

```bash
$ wget https://ftp.ncbi.nlm.nih.gov/bioproject/summary.txt
$ cat summary.txt | awk -F "\t" '{print $3}' | sed -E 's/[0-9]+//g' | sort | uniq -c
292 PRJDA
12933 PRJDB
497 PRJEA
57137 PRJEB
658381 PRJNA
1 Project Accession
```
**PRJNA** stands for "Project accession number" in NCBI.  
**PRJEB** stands for "European Bioinformatics Institute (EBI) Project Accession" and is specifically used to identify projects hosted by the EBI.   
**PRJEA** stands for "European Nucleotide Archive (ENA) Project Accession".    
**PRJDB** stands for "DNA Data Bank of Japan (DDBJ) Project Accession" and is specific to projects hosted by the DDBJ.    
**PRJDA** stands for "DNA Data Archive (DRA) Project Accession" and is used to identify projects hosted by the DRA.The DRA is an archive maintained by the National Bioscience Database Center (NBDC) in Japan and is part of the SRA consortium.    

4. UCSC
https://genome.ucsc.edu/goldenpath/help/ftp.html

## Gene Expression Omnibus (GEO)
ftp://ftp.ncbi.nlm.nih.gov/geo/
```bash
$ wget --recursive --no-parent -nd ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50499/suppl/
$ wget -r -np -nd -R "index.html*" ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50499/suppl/
```

## Sequence Read Archive

```bash
# There are four hierarchical levels of SRA entities and their accessions:
STUDY with accessions in the form of SRP, ERP, or DRP
SAMPLE with accessions in the form of SRS, ERS, or DRS
EXPERIMENT with accessions in the form of SRX, ERX, or DRX
RUN with accessions in the form of SRR, ERR, or DRR

## The fastq-dump command will only download the fastq version of the SRR, given the SRR number and an internet connection
$ fastq-dump SRR1013512
```

[HBCtraining](https://github.com/hbctraining/Accessing_public_genomic_data/blob/master/lessons/downloading_from_SRA.md)
