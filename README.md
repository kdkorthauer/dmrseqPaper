# Scripts for dmrseq Paper

This repository contains the code to reproduce the analyses and figures in "Detection and accurate False Discovery Rate control of differentially methylated regions from Whole Genome Bisulfite Sequencing", beginning with preprocessing steps, through the analyses with `dmrseq` and other tools used in the comparisons, and generation of figures in the main text and supplement.

The directory contains three folders:

1. **PREPROCESSING** - scripts to preprocess data from human tissue and mouse leukemia model datasets
2. **ANALYSIS** - code to reproduce the analysis of the human tissue and mouse leukemia model datasets as well as generate and analyze simulated data based on human dendritic cell data
3. **FIGURES** - code to generate the figures in the main text and supplement

More detail on the contents of each of these folders is provided below.

## Preprocessing

The scripts in this folder assume you have downloaded the necessary data files. The datasets used in the paper are all publicly available. Here is where you can find them for download:

* Roadmap (Human tissue samples):
* DNMT3A (Murine leukemia models):
* Dendritic cells (used to generate simulated datasets): - Note that there is no preprocessing script here since the processed data was downloaded from GEO and converted to bsseq objects in the R scripts in the **ANALYSIS** directory.

## Analysis

The scripts in this folder assume the preprocessing steps have been carried out.

## Figures

The scripts in this folder assume the analysis steps have been carried out, and read in the saved output files.
