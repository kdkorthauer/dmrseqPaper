#!/bin/bash

# change the following to where you want to download the data
# cd path/to/data/download

# this will download the DNMT3a expression data from GEO (GSE61969)
# the cufflinks gene counts table was used

wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61969/suppl/GSE61969_cufflinks.gene.counts.txt.gz'

# unzip the files

gunzip GSE61969_cufflinks.gene.counts.txt.gz


