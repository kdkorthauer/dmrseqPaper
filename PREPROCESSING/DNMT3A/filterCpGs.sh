 #!/bin/bash       

# bash script to extract the loci corresponding to CpG loci
# pattern matching on the context field to get both (1) NNCGN and (2) NCGNN
 
cd /n/irizarryfs01_backed_up/kkorthauer/WGBS/DNMT3A/DATA
 
for f in *.txt; 
 	do echo "Processing $f file.."; 
 	awk '$4~/^([AGCT][AGCT]CG[ACGT]|[ACGT]CG[ACGT][AGCT])$/' $f > CpG/${f}_CpG.txt
done;
