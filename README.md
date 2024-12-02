# Final-Respository (BIO 312 - Extra Credit)

## A Tutorial on Gene History Analysis of EEF2K
  ### Contents  
  1. [Introduction](#1-introduction)
  2. [BLAST](#2-BLAST)
  3. [Alignment](#3-Alignment) 
  4. [IQ-Tree](#4-IQ-Tree) 
  5. [Reconciling](#5-Reconciling)
  6. [Protein-Domain](#6-Protein-Domain)
  7. [Results](#7-Results)
  8. [Conclusion](#8-Conclusion)

# 1. Introduction

This repository consolidates analyses from Labs 3 to 8, focused on exploring the EEF2K gene family across various species. The primary goal is to investigate the evolutionary relationships, homolog distribution, and functional divergence of the EEF2K gene. Using bioinformatics tools like BLAST, sequence alignment, and phylogenetics, we identify high-confidence homologs, analyze their distribution across species, and infer evolutionary trends. This work aims to enhance our understanding of the EEF2K gene family's role in regulating protein synthesis and its functional adaptations in different vertebrate lineages.

**Before each lab** 

 Use this command to clone lab repository to Github

```
git clone https://github.com/Bio312/lab(insert lab)-$MYGIT
```

Within terminal (VS Code in this case), use this command to move to specific repository 

```
cd lab(insert lab)-$MYGIT
```
This is the same first step for all labs 3-8

Use this command to ensure at any time that you are in the proper working directory

```
pwd
```

# 2. BLAST

This section focuses on identifying homologous sequences of the EEF2K gene family across various species using BLAST. The goal was to locate high-confidence homologs to understand the gene family's evolutionary distribution. By creating a BLAST database from single-isoform proteomes, querying it with a representative EEF2K sequence, and filtering for significant hits, this analysis uncovered homologs across multiple species, laying the foundation for deeper studies.

## Create BLAST Database

Go to your lab 3 folder

```
cd ~/lab03-$MYGIT
```

Now, extraxt all compressed proteome files for analysis by decompressing

```
gunzip proteomes/*.gz
```

## Prepare BLAST Database

Now that we have extracted all proteome files for analysis, let us combine all individual proteome files into single FASTA file

```
cat proteomes/*.faa > allprotein.fas
```

We can now create the BLAST database from the combines protein sequences
```
makeblastdb -in allprotein.fas -dbtype prot
```

## Perform BLAST Search

Let us download the query protein sequence (EEF2K in this case) in FASTA format

```
ncbi-acc-download -F fasta -m protein "NP_037434.2"
```

Perform BLAST search using query protein against database to generate a typical output

```
blastp -db allprotein.fas -query NP_037434.2.fa -outfmt 0 -max_hsps 1 -out EEF2K.blastp.typical.out
```

Now, we can also perform another BLAST search to generate detailed tabular output (for easier analysis)

```
blastp -db ../allprotein.fas -query NP_037434.2.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out EEF2K.blastp.detail.out
```

## Filter BLAST Results

Now we filter our BLAST results to keep only high-confidence hits with an e-value less than 1e-30 using this command

```
awk '{if ($6 < 1e-30) print $1}' EEF2K.blastp.detail.out > EEF2K.blastp.detail.filtered.out
```

Count the number of high-confidence homologs identified after filtering

```
wc -l EEF2K.blastp.detail.filtered.out
```

## Analyze BLAST Results

Count how many paralogs found in each species 

```
grep -o -E "^[A-Z]\.[a-z]+" EEF2K.blastp.detail.filtered.out | sort | uniq -c
```

### See Results for further information on expected output of these commands








# 7. Results

## Lab 03

The BLAST analysis for EEF2K - after applying E-value threshold of 1e-30 - found 14 high-confidence homologs across the species in the database. A 11-species table should appear as below 

| Species      |   Count   |
|--------------|-----------|
| C.carcharias |    1       |
| C.mydas      |     1      |
| D.rerio      |     1       |
| E.caballus   |      1    |
| F.catus      |      1     |
| G.aculeatus  |      1     |
| G.gallus     |       1    |
| H.sapiens    |      2     |
| S.salar      |      2    |
| S.townsendi  |       1    |
| X.laevis     |      2     |

### Highest scoring hit identified as query sequence itself (_Homo sapeins_ EEF2K)

### No species found to have zero homologs

**IMPORTANT NOTE:** For your project, it is desirable to work with between 20 and 85 homologs. In addition, you should have between 2 and 15 copies from each species. Rest has pre-selected gene families that should be in these ranges. If you are not in that range, you will need to change the e-value threshold to increase or decrease the number of hits  - but you **must** talk to your TA when making this decision.  

























