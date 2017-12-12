#!/bin/bash

# change the following to where you want to download the data
# cd path/to/data/download

# this will download the Dendritic data from GEO (GSE64177)

wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64177/suppl/GSE64177_RAW.tar'

# untar and unzip the files

tar -xf GSE64177_RAW.tar
gunzip *gz 

