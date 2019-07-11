# Immune surveillance in clinical regression of pre-invasive squamous cell lung cancer

This repository provides code used to generate figures and results quoted in "Immune surveillance in clinical regression of pre-invasive squamous cell lung cancer" by Pennycuick, Teixeira et al. It is shared to assist researchers who wish to replicate our results or apply our analyses to their own datasets. Due to the complexity of pre-processing steps it is not intended to be a simple 'one-click' analysis, rather it is provided as a reference to provide full transparency for the reader to understand our analyses.

## Data

To run the code in this repository you will need to download raw data. This is shared in the following repositories:

* CIS gene expression data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108124
* Matched stromal expression data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133690
* Methylation data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108124
* Whole genome sequencing data: https://www.ebi.ac.uk/ega/datasets/EGAD00001003883
* IHC and imaging data: raw data not publicly available

The file load.data.R in this repository parses these raw data into the format required by the main analysis function, plot.figures.R. However, due to computational requirements and limitations in sharing potentially identifiable genomic data, we are not able to share full pipelines to do this. Full replication of this analysis does therefore require bioinformatic expertise to pre-process the above raw data.

Pre-processing steps are summarised as follows:

* Calling of mutations, indels, rearrangements and copy number alterations, as well as normalisation of gene expression and methylation data was performed as detailed here https://doi.org/10.1038/s41591-018-0323-0
* Detection of LOH in HLA regions was performed using LOHHLA, with manual curation of LOH calls: https://bitbucket.org/mcgranahanlab/lohhla/src/master/ ; https://doi.org/10.1016/j.cell.2017.10.001. The same methodology was applied to TCGA data for comparison.
* Neoantigens were predicted using netMHCpan as described here https://doi.org/10.1038/s41586-019-1032-7

## Analysis

The file `plot.figures.R` generates all figures and supplementary tables presented in our paper. In addition, it generates a file `results.in.text.txt` containing additional statistics referenced in the paper. The file should contain sufficient comments for a reader to follow.

## Software dependencies

Figures were produced using R version 3.5.0 and Bioconductor version 3.7. Required packages are listed in `plot.figures.R`. Full package version information is cached on running this file; the full package version list used to generate publication figures is included as `package.versions.csv`.

Excluded from this file is CIBERSORT.R and the MethylCIBERSORT R package, which is not available on Bioconductor. Code for CIBERSORT and MethylCIBERSORT can be downloaded from https://cibersort.stanford.edu and https://figshare.com/articles/MethylCIBERSORT_-_an_R_package_for_methylation_based_deconvolution_of_tumour_data_/6453650. The published analysis uses CIBERSORT version 1.0.4 and MethylCIBERSORT version 0.2.0.
