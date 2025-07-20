
### Script to run Lotus2 to analyse ITS2 data from demultiplexed NextSeq sequences


# Intstall and activate Lotus2

#conda create -n lotus2 lotus2
#conda activate lotus2


# run lotus 2, using UNITE V10 as the custom reference DB


lotus2 -i ../Data/fastq \
-m ../Data/Mapping_file.txt \
-o ../Output/lotus2_scats \
-refDB ../Data/UNITEv10/sh_general_release_dynamic_19.02.2025.fasta \
-tax4refDB ../Data/UNITEv10/sh_general_release_dynamic_19.02.2025.tax \
-amplicon_type ITS2 \
-LCA_idthresh 97,95,93,91,88,78,0 \
-tax_group fungi \
-taxAligner blast \
-forwardPrimer GTGARTCATCRARTYTTTG \
-reversePrimer CCTSCSCTTANTDATATGC \
-clustering vsearch \
-derepMin 8:1,4:2,3:3 \
-id 0.97 \
/



~               