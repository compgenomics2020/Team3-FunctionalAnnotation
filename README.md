# Functional Annotation
This pipeline combine several tools together to make functional annotation on gene prediction results.
***
## Members 
* Allison Rozanski <br />
* Pallavi Misra <br />
* Gulay Bengu Ulukaya <br />
* Cheng Shen-Yi <br />
***
## USEARCH <br />
We chose to cluster the genes based upon similarity in order to reduce the amount of overlap when annotating these genes. This is executed through the UCLUST algorithm. UCLUST preforms this by creating clusters that contain a single centroid sequence upon which the other sequences must have a certain sequence similarity to be considered apart of the cluster. We can set an identity threshold which can be thought of as the radius of the cluster.
***
## Homology Tools

### eggNOG mapper V2 <br /> 
eggNOG-mapper is a Tool for functional assignments based on precomputed orthologous clusters present it the eggNOG database. This is performed in the steps as follows

1) Sequence Mapping using HMMER3 or DIAMOND, for our analyses we use DIAMOND as a result of size of input and as because it is recommended over HMMER3 when annotating organisms with close relatives among the species covered by eggNOG

2) Orthology assignment

3) Functional Annotation, which is restricted to closest orthologs for reduction of false positives

### CARD-rgi <br />
Comprehensive Antibiotic Resistance Database (CARD) is a rigorously curated collection of characterized, peer-reviewed Antibiotic Resistance Genes which is monthly updated. Resistance Gene Identifier(RGI) is a toolkit based on CARD for annotating Antimicrobial genes.

### VFDB <br />
Virulence Factor Database (VFDB) is an integrated and comprehensive online resource for curating information about virulence factors of bacterial pathogens (recently updated in 2019). The database contains information such as structure features of the virulence factors, functions and mechanisms used by the pathogens for circumventing host defense mechanisms and causing pathogenicity. Core dataset of DNA sequences was downloaded from VFDB website, which include genes associated with experimentally verified Virulence Factors only. BLAST database was build based on the downloaded dataset from VFDB and BLASTN was used.
***
## Ab-initio Tools <br />

### SignalP <br />
Signal peptides are short amino acid sequences of newly synthesized proteins that target proteins into, or across, membranes. SignalP 5.0 distinguishes three types of signal peptides in prokaryotes: Sec substrates cleaved by SPase I (Sec/SPI), Sec substrates cleaved by SPase II (Sec/SPII), and Tat substrates cleaved by SPase I (Tat/SPI). SignalP consists of two different predictors based on neural network and hidden Markov model algorithms. In order to predict potential signal peptides of proteins, the D-score from the SignalP output is used for discrimination of signal peptide versus non-signal peptide.

### HMMTOP <br />
HMMTOP (Hidden Markov Model for Topology Prediction) transmembrane topology prediction tool predicts both the localization of helical transmembrane segments, their start and end positions in the sequence, and the topology of transmembrane proteins. HMMTOP method is based on the transmembrane proteins determined by the maximum divergence of amino acid composition of sequence segments.

### Piler-CR <br />
CRISPR are family of DNA sequences found in the genomes of prokaryotic organisms- bacteria and archaea. They are derived from DNA fragments of viruses that had previously infected the prokaryote and provides protection from viruses and plays a major role in antiviral defense system. PILERCR identifies CRISPR repeats by using BLAST to find their fragmented/ degraded copies. A CRISPR array is found when it fulfills the criteria of having a set of CRISPR repeats with intervening unique sequences known as spacers. This program provides fast identification and classification of CRISPR genes and also has both high sensitivity and high specificity.
***

## Environment Instructions
[Environment_Instructions](https://github.gatech.edu/compgenomics2020/Team3-FunctionalAnnotation/blob/master/Environment_Instructions)
***

## Script Execution
`pipeline.py -e <path/to/eggNOG> -p <path/to/Piler-CR> -s <path/to/signalP> -i <path/to/input>`<br />
***
