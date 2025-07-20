#!/bin/bash

##### Quality check with FastQC and MultiQC


### How to install FastQC and MultiQC with conda

# Create new conda environment: conda create --name fastqc
# Activate fastqc environment: conda activate fastqc
# Install fastqc: conda install -c bioconda -c conda-forge fastqc
# Install multiqc: conda install -c bioconda -c conda-forge multiqc


### Running the script

# Activate fasqc environment: conda activate fastqc

# Create folders for results
mkdir ../Output/fastqc_report_testdata


# Analyse ITS reads (1 runs)
# *R1.fastq = files that end with this


for f in ../Data/fastq/*R1.fastq.gz ../Data/fastq/*R2.fastq.gz;
do fastqc $f --outdir=../Output/fastq;
done


# Visualize the quality analyses for each locus and direction with MultiQC
multiqc ../Output/fastq/*R1_fastqc.zip -n multiqc_report_R1 -o ../Output/fastq
multiqc ../Output/fastq/*R2_fastqc.zip -n multiqc_report_R2 -o ../Output/fastq


