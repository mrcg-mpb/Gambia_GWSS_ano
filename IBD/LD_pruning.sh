#!/bin/bash

# Generate a genetic plink file from the VCF file
plink --vcf 3R_chrom_renamed.vcf.gz --make-bed --out ./LD_process/3R_data

# Perform LD pruning
plink --bfile ./LD_process/3R_data --indep-pairwise 50 5 0.2 --out ./LD_process/pruned

# Extract the pruned SNPs
plink --bfile ./LD_process/3R_data --extract ./LD_process/pruned.prune.in --make-bed --out ./LD_process/3R_data_pruned

# Convert the pruned data back to VCF format
plink --bfile ./LD_process/3R_data_pruned --recode vcf --out ./LD_process/3R_chrom_pruned

# Compress the VCF file
bgzip ./LD_process/3R_chrom_pruned.vcf
