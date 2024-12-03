# Final-Respository (BIO 312 - Extra Credit)

## A Tutorial on Gene History Analysis of EEF2K
  ### Contents  
  1. [Introduction](#1-introduction)
  2. [BLAST](#2-BLAST)
  3. [Multiple Sequence Alignment](#2-Multiple-Sequence-Alignment)
  4. [IQ-Tree](#4-IQ-Tree) 
  5. [Reconciliation](#5-Reconciliation)
  6. [Protein Domain Prediction](#6-Protein-Domain-Prediction)
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

### See [Results](#7-Results) for further information on expected output of these commands

# 3. Multiple Sequence Alignment 
This repository focuses on performing multiple sequence alignment (MSA) for the EEF2K gene family, a continuation from Lab 3's homolog identification. The aim is to align sequences with tools such as MUSCLE, AlignBuddy, and T-Coffee, assess sequence conservation, and calculate alignment statistics. Key analyses include identifying conserved regions, calculating average percent identity, and visualizing alignments to ensure their biological validity.

## Set Up Environment

Make sure you create a folder for the EEF2K family for this lab

```
mkdir ~/lab04-$MYGIT/EEF2K
```

Now, go into the directory

```
cd ~/lab04-$MYGIT/EEF2K
```

Use the ```pwd``` command if you are unsure whether you are in the right working directory

## Extract Homolog Sequences
Extract sequences for homologs identified in lab 03 and save them as ```EEF2K.homologs.fas``` using this command 

```
seqkit grep --pattern-file ~/lab03-$MYGIT/EEF2K/EEF2K.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.fas
```

## Perform Multiple Sequence Alignment
Align homolog sequence using MUSCLE to create hypothesis of homology among positions in sequence using this command

```
muscle -align ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.fas -output ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.al.fas
```

## Generate Alignment PDF
Use R script to create high-resolution PDF of alignment 

```
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.al.fas
```
PDF should have the name ```EEF2K.homologs.al.fas.pdf```


## Obtain Alignment Statistics 
Using _alignbuddy_, we can calculate the total width (number of columns) in the alignment, the number of columns containing gaps, and the the number of columns in the alignment that are invariant. Here are the three separate codes for each command

```
alignbuddy -al ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.al.fas
```

```
alignbuddy -trm all ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.al.fas | alignbuddy -al
```

```
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.al.fas | alignbuddy -al
```

We can use _t_coffee_ to calculate average percent identity for alignment

```
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.al.fas -output sim
```

Once again, _alignbuddy_ can also be used to calculate average percent identity. Here is the command 

```
alignbuddy -pi ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.al.fas | awk '(NR>2) { for (i=2;i<=NF;i++){ sum+=$i;num++} } END{ print(100*sum/num) }'
```

### See [Results](#7-Results) for further information on expected output of these commands

# 4. IQ-TREE
This repository focuses on constructing a maximum-likelihood phylogenetic tree for the EEF2K gene family using sequence alignments derived from Lab 4. Key steps include refining the sequence alignment by removing duplicates, selecting an appropriate substitution model, and estimating the optimal phylogenetic tree with bootstrap support using IQ-TREE. Additionally, the phylogeny was visualized both as an unrooted and midpoint-rooted tree, with graphical outputs saved as PDFs for detailed examination. The analysis provides insights into the evolutionary relationships among EEF2K homologs across species, emphasizing branch lengths, clade support, and overall topology.

## Set Up Environment

Make sure you create a folder for the EEF2K family for this lab

```
mkdir ~/lab05-$MYGIT/EEF2K
```

Now, go into the directory

```
cd ~/lab05-$MYGIT/EEF2K
```

Use the ```pwd``` command if you are unsure whether you are in the right working directory

## Refine Sequence Alignment
First, remove duplicate sequence labels and save refined file in directory using this command

```
sed 's/ /_/g' ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.al.fas | seqkit grep -v -r -p "dupelabel" > ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.fas
```
## Generate Maximum-Likelihood Tree
Run IQ-Tree to estimate the maximum-likelihood tree for EEF2K homologs using refined alignment 

```
iqtree -s ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.fas -bb 1000 -nt 2
```

Use R script to generate a graphical PDF of unrooted tree with adjusted label sizes and lengths

```
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.fas.treefile ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.fas.treefile.pdf 0.4 15
```

PDF name should be ```EEF2K.homologsf.al.fas.treefile.pdf```

## Generate Midpoint-Rooted Tree
Perform the midpoint rooting by rooting at the midpoint of the longest branch using this command

```
gotree reroot midpoint -i ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile
```

Generate an SVG image of the midpoint-rooted tree

```
nw_order -c n ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s > ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.svg
```

Convert SVG image into PDF for easier sharing and analysis

```
convert ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.svg ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.pdf
```

PDF name should be ```EEF2K.homologsf.al.mid.treefile.pdf```

## Generate Cladogram from Rooted Tree
Create cladogram from rooted tree and save it as SVG

```
nw_order -c n ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.midCl.treefile.svg
```

Convert SVG image into PDF 

```
convert ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.midCl.treefile.pdf
```

PDF name should be ```EEF2K.homologsf.al.midCl.treefile.pdf```


### See [Results](#7-Results) for further information on expected output of these commands

# 5. Reconciliation

The goal of this lab was to reconcile the evolutionary history of the EEF2K gene family by comparing its gene tree to the species tree of the studied organisms. This process involved identifying gene duplication, loss, and speciation events, which provide insight into the evolutionary dynamics of the EEF2K gene family. 

## Set Up Environment

Make sure you create a folder for the EEF2K family for this lab

```
mkdir ~/lab05-$MYGIT/EEF2K
```

Now, go into the directory

```
cd ~/lab05-$MYGIT/EEF2K
```

Use the ```pwd``` command if you are unsure whether you are in the right working directory

## Copy Midpoint Tree from Lab 05 to Lab 06

```
cp ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile ~/lab06-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile
```

## Run Reconciliation using NOTUNG
Use NOTUNG to identify duplication, losses, and speciation events through this command

```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab06-$MYGIT/EEF2K/species.tre -g ~/lab06-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/EEF2K/
```


## Convert Reconciliation
Converts Notung reconciliation into RecPhyloXML format for visualization

```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.rec.ntg --include.species
```

## Visualize Reconciliation
Generate a visual representation of the gene-within-species reconciliation tree

```
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.rec.ntg.xml -o ~/lab06-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.rec.svg
```

## SVG to PDF conversion
Format to PDF for easier viewing

```
convert -density 150 ~/lab06-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.rec.pdf
```
see both png and pdf versions ```EEF2K.homologsf.al.mid.treefile.rec..png``` and ```EEF2K.homologsf.al.mid.treefile.rec.pdf```

## Examine Reconciliation
Displays summary table of duplication, loss, and speciation events inferred during reconciliation using this command

```
cat ~/lab06-$MYGIT/EEF2K/EEF2K.homologsf.al.mid.treefile.rec.events.txt
```
should find a text file ```EEF2K.homologsf.al.mid.treefile.rec.events.txt```

### See [Results](#7-Results) for further information on expected output of these commands

# 6. Protein Domain Prediction

This lab focuses on identifying protein domains in the EEF2K gene family using RPS-BLAST with the Pfam database to predict functional domains. The raw protein sequences are processed to remove stop codons and then analyzed to identify conserved and lineage-specific domains. The predicted domains are plotted alongside a phylogenetic tree generated in previous labs to visualize their distribution across species. This analysis provides insights into the functional conservation and potential evolutionary adaptations of the EEF2K gene family. The lab also involves comparing domain annotations across species and assessing their significance by adjusting e-values. These results, combined with previous phylogenetic findings, contribute to the results section of the final research paper, offering a deeper understanding of the evolutionary dynamics and functional roles of EEF2K across lineages.

## Set Up Environment

Make sure you create a folder for the EEF2K family for this lab

```
mkdir ~/lab05-$MYGIT/EEF2K
```

Now, go into the directory

```
cd ~/lab05-$MYGIT/EEF2K
```

Use the ```pwd``` command if you are unsure whether you are in the right working directory

## Input Protein Sequences
Make a copy of raw unaligned sequence from lab04 to remove any stop codons (asterisks) from raw sequences and save clean file with this command

```
sed 's/*//' ~/lab04-$MYGIT/EEF2K/EEF2K.homologs.fas > ~/lab08-$MYGIT/EEF2K/EEF2K.homologs.fas
```

## Run RPS-BLAST
With clean sequence, predict Pfam domains using RPS-BLAST with stringent e-value threshold for accuracy 

```
rpsblast -query ~/lab08-$MYGIT/EEF2K/EEF2K.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/EEF2K/EEF2K.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```

## Copy Phylogenetic Tree File and Plot Pfam Domains

copy the phylogenetic tree file from lab05 to current directory

```
cp ~/lab05-$MYGIT/EEF2K/EEF2K.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/EEF2K/
```

Use Rscript to generate a PDF that visualizes Pfam domains alongside the phylogenetic tree for EEF2K

```
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/EEF2K/EEF2K.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/EEF2K/EEF2K.rps-blast.out ~/lab08-$MYGIT/EEF2K/EEF2K.tree.rps.pdf
```
The file created here is ```EEF2K.tree.rps.pdf```

### See [Results](#7-Results) for further information on expected output of these commands


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

### Highest scoring hit
identified as query sequence itself (_Homo sapiens_ EEF2K). 

Note: _Xenopus laevis_ and _Gallus gallus_ had very high scores highlighting strong conservation of EEF2K gene across mammals, amphibians, and birds

### Zero Homolog Test
No species found to have zero homologs

### Important Note
As discussed with TA, this table was approved due to high degree of conservation of the EEF2K gene family. Decreases in stringency found hits with significantly lower percent identity (~30%)

## Lab 04

| Metric      |   Value   |
|--------------|-----------|
| Alignment length |    853      |
| Columns in alignment containing gaps      |     178   |
|Columns in alignment that are invariant    |     292    |
|Average percent identity (_t_coffee_)  |     72.46   |
|Average percent identity ( _alignbuddy_)   |     67.8115    |

The alignment results demonstrate a high degree of conservation in the EEF2K gene family, emphasizing its critical role in cellular regulation. The high percent identity and the significant number of invariant columns provide strong evidence for the functional conservation of EEF2K across species. The presence of gaps and variable regions suggests that while the core functional regions of EEF2K are conserved, some evolutionary divergence has occurred in non-essential regions.

### PDF of Rscript command can be found within its respective folder (Lab 04) in final repository

## Lab 05

### Phylogenetic Tree Estimation Parameters

| **Parameter**                | **Value**                       |
|-------------------------------|---------------------------------|
| **Model of Substitution**     | JTT+G4                         |
| **Model of Rate Heterogeneity** | Gamma with 4 categories         |
| **Gamma Shape (α)**            | 0.6404                          |

| **Category** | **Relative Rate** | **Proportion** |
|--------------|-------------------|----------------|
| 1            | 0.0615            | 0.2500         |
| 2            | 0.3340            | 0.2500         |
| 3            | 0.9006            | 0.2500         |
| 4            | 2.7039            | 0.2500         |

*Relative rates are computed as the mean of the portion of the Gamma distribution falling in each category.*

### Analysis of Bootstrap values
Generally high bootstrap values (lowest being 72) which shows strong general conservation across EEF2K homologs. Mammals form distinct clade with high bootstrap value, consistent with evolutionary relationship. Reptiles and amphibiand more distantly related to mammals, as evidence by branching pattern. Relatively short branch lengths within clades imply slower rate of divergence, supporting functional conservation. Some lower bootstrap values (e.g., 72 between Gasterosteus aculeatus and Salmo salar) highlight regions of uncertainty, potentially reflecting rapid evolutionary events or incomplete lineage sorting.

### PDF of Rscript command can be found within its respective folder (Lab 05) in final repository

## Lab 06

| **Parameter**                                    | **Value**                                  |
|-------------------------------------------------|--------------------------------------------|
| **Cost of Reconciled Tree**                     | 12.0                                      |
| **Copies in Common Ancestor of Gnathostomes**   | 1                                         |
| **Events Leading to Common Ancestor of Tetrapods** | 1 speciation, 1 gene loss                 |
| **Total Duplications**                          | 4                                         |
| **Total Losses**                                | 6                                         |

Four gene duplication events are suggestive of functional diversification in the EEF2K gene family. Six gene losses suggest that EEF2K may not be essential in all lineages or has been replaced by alternative mechanisms. However, the single copy in Gnathostomata is suggestive that EEF2K originated in a common ancestor and was conserved. 

### PDF of Rscript command can be found within its respective folder (Lab 06) in final repository

## Lab 08

| **Parameter**                                      | **Value**                                                                                     |
|----------------------------------------------------|-----------------------------------------------------------------------------------------------|
| **Gene with Most Pfam Domains**                   | *H.sapiens_EEF2K_eukaryotic_elongation_factor_2_kinase* (14 domains)                          |
| **Pfam Domains Found in Most Genes**              | *pfam02816 (Alpha_kinase, Alpha-kinase family)*                                               |
| **Pfam Domains Found Multiple Times in a Gene**   | *H.sapiens_EEF2K_eukaryotic_elongation_factor_2_kinase*: *pfam02816*                          |
| **Gene with Longest Pfam Domain Annotation**      | *S.townsendi_EEF2K_LOW_QUALITY_PROTEIN_eukaryotic_elongation_factor* (183 residues)           |
| **Domain with Best E-value**                      | *D.rerio_eef2k_eukaryotic_elongation_factor_2_kinase* (*pfam02816*, E-value: 2.52e-74)        |
| **Domains Present Across All Lineages**           | *pfam02816 (Alpha_kinase)*                                                                    |
| **Lineage-Specific Domain**                       | *pfam01270 (Glyco_hydro_8)* found in *X.laevis*, indicating potential functional adaptation.   |
| **Evidence of Functional Divergence**             | Some gene copies lack *pfam01270*, suggesting divergence in functional roles across lineages. |
| **Evidence of Functional Gain**                   | Acquisition of *pfam01270* in *X.laevis* indicates possible lineage-specific functional gain.  |

The Alpha_kinase (pfam02816) domain is highly conserved and ubiquitous across the analyzed gene copies, reaffirming its critical role in EEF2K functionality. Lineage-specific differences, such as the presence of Glyco_hydro_8 (pfam01270) in X.laevis, highlight potential adaptations that may be explored further in functional and evolutionary contexts.


### PDF of Rscript command can be found within its respective folder (Lab 08) in final repository

# 8. Conclusion 

The work conducted from Lab 3 to Lab 8 provides a detailed evolutionary and functional analysis of the EEF2K gene family. Lab 3 focused on aligning homologous sequences, identifying conserved and variable regions. Lab 4 involved constructing a maximum-likelihood phylogenetic tree, revealing relationships between homologs. Lab 5 refined the phylogeny through rooting and bootstrap analysis, establishing divergence patterns. Lab 6 reconciled the gene and species trees, identifying key duplication and loss events, shedding light on lineage-specific adaptations. Finally, Lab 8 integrated domain predictions with the phylogeny, confirming conserved Alpha_kinase domains across species and highlighting lineage-specific gains, such as the Glyco_hydro_8 domain in X. laevis. Together, these analyses reveal how evolutionary pressures shaped the EEF2K gene family, providing insights into its conserved roles and lineage-specific adaptations.


