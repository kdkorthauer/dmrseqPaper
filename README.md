# Scripts for dmrseq Paper

This repository contains the code to reproduce the analyses and figures in "Detection and accurate False Discovery Rate control of differentially methylated regions from Whole Genome Bisulfite Sequencing", beginning with preprocessing steps, through the analyses with `dmrseq` and other tools used in the comparisons, and generation of figures in the main text and supplement. These scripts rely on the package dmrseq from [github](https://github.com/kdkorthauer/dmrseq/tree/01b82b69654755e5b991041ac6279fb68757fa75) (commit `01b82b69654755e5b991041ac6279fb68757fa75` on December 12, 2017).

The directory contains three folders:

1. **PREPROCESSING** - scripts to download and preprocess data from human tissue, dendritic cell line, and mouse leukemia model datasets
2. **ANALYSIS** - code to reproduce the analysis of the human tissue and mouse leukemia model datasets as well as generate and analyze simulated data based on human dendritic cell data
3. **FIGURES** - code to generate the figures in the main text and supplement
4. **BENCHMARK** - code to download dendritic cell data and generate simulated benchmark datasets for two sample sizes

To reproduce all results from the paper, run the code in each folder in the order above (i.e. the scripts in **FIGURES** assume that the code in **ANALYSIS** has been run. The scripts in **PREPROCESSING** and **ANALYSIS** are set up to use the batch processing system _Slurm_ with `.sh` files that call the `.R` analysis files with various arguments you need to specify in the headers. You may the `.sh` scripts them to use a different batch processing system, or change the `SLURM` argument in the `.sh` files to `F` to run each step separately (in serial). More detail on the contents of each of these folders is provided below.

## Preprocessing

The scripts in this folder will downloaded and preprocess the necessary data files. The datasets used in the paper are all publicly available. Details for each subfolder are described in the next subsections. The scripts used for downloading and preprocessing are a mix of bash scripts (`.sh`) and R scripts (`.R`). In the header of each script, there are variables that you should change to reflect your system (e.g. path to where SRATools is installed, path to where you’d like to download each dataset, etc). It assumes that you already have SRATools, Bismark, cutadapt, and the various necessary R packages installed. 

### `ROADMAP`: Human tissue samples

#### WGBS data 

The list of SRA accession numbers in the file `SRR_Acc_List_public_heart_intestine.txt` was used with the SRATools `prefect` command to obtain `.sra` files for the selected runs from SRP000941.

1. Download SRA files: `prefetch.public` downloads the SRA files in the list `SRR_Acc_List_public_heart_intestine.txt` using prefetch from the SRAToolkit. This list was generated using the SRA run selector to include all runs from the selected tissues (heart right ventricle, heart left ventricle, small intestine, and sigmoid colon). 

2. Prepare Bismark reference genome: `assembleRefGenome_public.sh` will download the reference assembly for GRCh38 from the NCBI ftp and then build the reference genome using the bismark_genome_preparation function. It only needs to be done once before mapping.

3. Map and extract cytosine counts: `fastqToBismarkExtractCytosine_public.slurm` calls the `fastqToBismarkExtractCytosine_public.sh` script, which is the major workhorse from fastq-dump to obtain fastq from sra, to mapping the reads and extracting cytosine counts using Bismark, and then converting to bsseq Data objects by calling the `cytosineToBSseq_public.R script`.

#### RNA-seq data 

Download and process RNA-seq data: the `RNAseq_betToBam.sh` script will download the mapped matched RNA-seq bed files from NCBI FTP (specifically, the bed.gz files for mRNA-Seq for the tissues of interest: ‘heart_left_ventricle`, `heart_right_ventricle`, `small_intestine`, and `sigmoid_colon`), as well as the hg19 human gif annotation file from Gencode version 19. Then it will use bedrolls to convert bed files to bam files, which will then be counted using featureCounts. The featureCounts output is then combined into one table for all samples, which is written to a text file for later use.

### `DNMT3A`: Murine leukemia models

#### WGBS data 

1. Download raw WGBS data and prepare genome assembly: The `getSRA_prepccessBismark.sh` script downloads the SRA files in experiment id SRP074072 then prepares the mm10 genome assembly for Bismark. This script uses the `SRR_Acc_List_Url.txt` file which was downloaded from the SRA website. 

2. Map reads and extract cytosine information: The next processing steps are carried out with the `runBismark.sh` bash script which 
(1) utilizes sratoolkit to convert the SRA files to fastq format and perform QC
(2) utilizes trim galore to trim adapter sequences from reads
(3) utilizes bismark to map reads onto reference genome and create bam files
(3) utilizes samtools to merge bam files
(4) utilizes bismark to extract the count of methylated and unmethylated reads for each CpG in the genome
(5) utilizes bismark to extract the cytosine reports, which are the same as in (4) except that strand info is preserved, which allows for CpGs on opposite strands to be collapsed into one site.
(6) Finally, the R script `cytosineToBSseq.R` is called, which uses modified versions of functions in the bsseq Bioconductor package to read in the Bismark output and create bsseq objects, which are saved as RData objects in a specified data directory.

This script can be processed in batch using the Slurm scheduler with the `runBismark.slurm` script. Alternatively, this file can be ignored and `runBismark.sh` can be run in serial.

#### RNA-seq data

Download processed RNA-seq data: the `download_RNAseq.sh` script will download the processed matched RNA-seq count table from GEO accession number GSE61969 (specifically, the cufflinks gene count table `GSE61969_cufflinks.gene.counts.txt`).


### `DENDRITIC`: Dendritic cells (used to generate simulated datasets)

#### WGBS data  

Download the data with the `download.sh` script. Note that there are no further preprocessing steps analogous to the other two datasets since the processed data in the form of text files was downloaded from GEO (GSE64177). The text files are read in, converted to bsseq objects, and used to generate simulated data in the R scripts in the **ANALYSIS** directory. 

#### RNA-seq data 

Since the dendritic cell line data was used to generate null comparisons and simulated data, no RNA-seq data was used for this set.

## Analysis

The scripts in this folder assume the preprocessing steps have been carried out. The `.sh` scripts are set up to use the batch processing system _Slurm_ with `.sh` files that call the `.R` analysis files with various arguments you need to specify in the headers. These arguments include: (1) file paths to where data from the previous step is stored and results are output, (2) path to where metilene is installed on your system, (3) whether the alternative methods should be run in default or the optimized settings as described in the paper (to reproduce all results, you’ll need to run both versions), and (4) whether to use a naive test statistic (`area` or `avg`) for dmrseq instead of the default as described in the paper (to reproduce all results, you’ll need to run both versions). You may the `.sh` scripts them to use a different batch processing system, or change the `SLURM` argument in the `.sh` files to `F` to run each method and dataset separately (in serial). 

## Figures

The scripts in this folder assume the analysis steps have been carried out, and read in the saved output files in order to generate summary figures. 

- `SimulationResults.R` generates Figure 1b, Figure 2, and Supplementary Figures S3-S5
- `EmpiricalResults_byWindow.R` generates Figure 3, and Supplementary Figures S6-S9
- `fig4.Rmd` calls the script `rankComparisonFigure.R` to generate Figure 4
- `SupplementaryFigureS1.R` generates Supplementary Figure S1
- `DataDescription.R` generates Supplementary Figure S2
- `fig1a.pptx` generates the schematic illustration for Figure 1a
- `SimulationResults_subset.R` - generates Supplementary Figures 10-12 

## Benchmark

This directory contains a separate stand-alone R script that will generate the simulated benchmark data for two sample sizes (2 per condition and 3 per condition). This simulated data is based on the dendritic cell data, which the script downloads from GEO. After reading in the text files from GEO and converting them into a bsseq object, the `simDMRs` function in the `dmrseq` package is used to add 3000 DMRs. Differences are added between the first half of samples and the second half of samples. The simulated data is then saved as an `.rds` file which contains a bsseq object as well as a Granges object denoting the true DMRs. Output is also saved in the form of compressed tab-delimited text files for each sample containing methylated and unmethylated counts per CpG, in addition to a `.bed` file containing the DMRs. These outputs are also provided in a [FigShare repository for direct download](https://figshare.com/projects/Whole_Genome_Bisulfite_Sequencing_WGBS_Benchmark_Data/27532). The `simDMRs` function in the dmrseq package may also be used to generate additional benchmarking datasets starting from other sets of biological replicates by following a similar approach.

