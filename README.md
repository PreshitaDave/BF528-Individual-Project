# BF528-Individual-Project


In this project, I aim to attempt to replicate the Baron et al. results using a 51 year old female human donor’s cells. The reads and barcodes would be processed and then aligned to a reference transcriptome to generate read counts for unique molecular identifiers (UMIs). Clustering will then be performed to determine the number of cell types found. 

Reference: Baron, Maayan, Adrian Veres, Samuel L. Wolock, Aubrey L. Faust, Renaud Gaujoux, Amedeo Vetere, Jennifer Hyoje Ryu, et al. 2016. “A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-Cell Population Structure.” Cell Systems 3 (4): 346–60.e4.

# Repository contents
The data curator folder contains the following scripts:

1. counts.qsub: Extracts barcodes from each run, and counts the number of reads per distinct barcode. 


2. hist.R: For each run, plots the cumulative distribution of reads and retrieves whitelist.


3. index.qsub: Creates the index of the reference genome file for Human Release 40 (GRCh38.p13).


4. mapping.qsub: Maps the transcript to gene. 


5. salmon.qsub: Runs salmon alevin to create UMI matrix. 


The programmer folder contains the script:

1. program.R: Preprocesses the UMI Counts Matrix and saves the rds object with cluster information. 


Analyst:

analyst.R: This script analyses and identifies marker genes for clusters that were created.
