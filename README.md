# Bio-databases
Highly used databases in bioinformatics

## 1 Reference database

To download reference data, there are some different sources available:
General biological databases: Ensembl, NCBI, and UCSC
Organism-specific biological databases: Wormbase, Flybase, TAIR etc. (often updated more frequently, so may be more comprehensive)
Reference data collections: Illumina's iGenomes, one location to access genome reference data from Ensembl, UCSC and NCBI

### 1.1 Ensembl

```bash 
Ensembl identifiers
When using Ensembl, note that it uses the following format for biological identifiers:

ENSG###########: Ensembl Gene ID
ENST###########: Ensembl Transcript ID
ENSP###########: Ensembl Peptide ID
ENSE###########: Ensembl Exon ID
For non-human species a suffix is added:

ENSMUSG###: MUS (Mus musculus) for mouse
ENSDARG###: DAR (Danio rerio) for zebrafish

```

```bash
human genome assemblies
GRCh38 (aka hg38)
many rare/private alleles replaced
www.ensembl.org
Most up-to-date and supported

GRCh37 (aka hg19)
Some large gaps
grch37.ensembl.org
Limited data and software updates
Still the preferred genome of the clinical community

NCB136 (aka hg18)
many gaps
ncbi36.ensembl.org
No longer updated

```
### 1.2 NCBI
https://www.ncbi.nlm.nih.gov/home/download/

####  1.2.1 GeneBank   
  
**GCA** for GenBank assemblies   
**GCF** for RefSeq assemblies   
GCA(or GCF) is followed by an underscore and 9 digits. GRCh38.p11 is GCA_000001405.26.      

####  1.2.2 BioProject

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


####  1.2.3 Blast

```bash
$ nr.gz 非冗余蛋白数据库,来源于Swisprot,PIR,PDF,PDB和RefSeq
$ nt.gz 核酸数据库,来源于GeneBank, EMBL, DDBJ.部分非冗余.
$ # BLAST+程序包，提供了一个脚本update_blastdb.pl可以很方便的进行本地化下载blast数据库
$ perl update_blastdb.pl # 可以很方便的查看有哪些数据库支持下载

$ nohup perl update_blastdb.pl --decompress nt &> update.log & # 后台下载，支持断点续传
```


### 1.3 UCSC
加州大学圣克鲁兹分校基因组数据库
https://genome.ucsc.edu/goldenpath/help/ftp.html

### 1.4 CNGB
国家基因库
```bash
$ wget -c -nH -np -r -R "index.html*" --cut-dirs 4 ftp://ftp.cngb.org/pub/CNSA/data1/CNP0000405/CNS0064392/
```

### 1.5 GATK broad institute
- hg19
```bash
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/*.gz # 一次性全部下载
$ 
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3_hg19_pop_stratified_af.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.sites.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.sites.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.knowledgebase.snapshot.20131119.hg19.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.knowledgebase.snapshot.20131119.hg19.vcf.idx.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.dict.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.fai.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz
```

- hg38
```bash
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/*.gz # 全部下载
$
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_144.hg38.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3_grch38_pop_stratified_af.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

### 1.6 Human reference data

#### CHM13
目前最新的T2T基因组
- URL: https://github.com/marbl/CHM13
- CONTENTS: chr1-22(CHM13),chrX(CHM13),chrY(NA24385),chrM(CHM13)
- CITATION : Nurk, S., Koren, S., Rhie, A., Rautiainen, M., Bzikadze, A. V., Mikheenko, A., ... & Phillippy, A. M. (2022). The complete sequence of a human genome. Science, 376(6588), 44-53.

#### GRCh39
```bash

```
#### GRCh38

#### hg38

#### CRCh37

#### hg19

#### hs37d5

#### humanG1Kv37

#### b37

#### Human pangenome reference consortium

#### 1000 Genome Project 30X high coverage (hg38)

#### 1000 Genome Project Phase 3 (hg19)

#### HAPMAP3


## 2 Gene Expression Omnibus (GEO)
ftp://ftp.ncbi.nlm.nih.gov/geo/
```bash
$ wget --recursive --no-parent -nd ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50499/suppl/
$ wget -r -np -nd -R "index.html*" ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50499/suppl/
```

## 3 Sequence Read Archive (SRA)

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


## 4 dbSNP
[dbSNP of NCBI](https://www.ncbi.nlm.nih.gov/snp/)

```bash
rsid (reference SNP)
dbSNP Reference SNP (rs or RefSNP) number is a locus accession for a variant type assigned by dbSNP
```
download path : https://ftp.ncbi.nih.gov/snp/latest_release/VCF   
web search: https://www.ncbi.nlm.nih.gov/snp/rs328   

## 5 lncRNA 
长链非编码RNA [LNCipedia](https://lncipedia.org/)
LNCipedia数据库注释的是human基因组，目前版本包括 127,802 transcripts 和 56,946 genes。
```bash
$ less -S lncipedia_5_2_hg38.gtf | cut -f 9 | grep -v "^#" | awk -F ";" '{print $1}' | grep gene_id | awk -F " " '{print $2}' | sort | uniq | wc
   56946   56946  704558
(base) [09:34:00] zhusitao zhusitaodeMacBook-Air ~/Downloads 
$ less -S lncipedia_5_2_hg38.gtf
(base) [09:34:12] zhusitao zhusitaodeMacBook-Air ~/Downloads 
$ less -S lncipedia_5_2_hg38.gtf | cut -f 9 | grep -v "^#" | awk -F ";" '{print $2}' | grep transcript_id | awk -F " " '{print $2}' | sort | uniq | wc
  127802  127802 1800082
```

## 6 Single Cell Portal
Broad Institute单细胞测序数据，目前收录645 studies，包括 40,522,660 cells。[single cell](https://singlecell.broadinstitute.org/single_cell)。


