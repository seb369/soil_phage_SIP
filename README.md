# Metagenomic stable isotope probing reveals bacteriophages involved in soil carbon cycling
Analysis and data from a study using metagenomic-SIP to examine soil phages and the soil carbon cycle.

Barnett, S.E. & Buckley, D.H. (2023) Metagenomic stable isotope probing reveals bacteriophage participation in soil carbon cycling. *Environmental Microbiology*, 25:1785–1795. https://doi.org/10.1111/1462-2920.16395

## Files

### R analyses
This directory contians all R based analyses, all writen in rmarkdown (.rmd) with GitHub formatted markdown output (.md) for easy reading and figure visualization.

* **Phage_metagenomicSIP:** metagenomic-SIP analysis and identification and examination of vOTUs.
* **Phage_abundance_dynamics:** analysis of vOTU and Streptomyces rpoB qPCR and amplicon sequencing.

### Data
This directory contains important data not included in the manuscript supplemental, specifically the rpoB OTU data. The raw sequencing data can be found on the NCBI SRA under BioProject PRJNA746975.

* **final.otu.fasta:** rpoB OTU consensus sequences.
* **final.otutab.biom:** rpoB OTU table in biom format.
* **final.otu.fasta.blca.out:** rpoB OTU taxonomic classifications using BLCA and a custom rpoB database.
