#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - LARGE queue      |
#   |    - 20 CPU + 20 Gb   |
#   +-----------------------+
#
#  Replace the USER name in this script with your username and
#  call your project whatever you want
#
#  This script must be made executable like this
#    chmod +x my_script
#
#  Submit this script to the queue with a command like this
#    run_script my_script.sh
#
#  My preferred setup before running:
#    -- script to be run in /home/scripts
#    -- project directory (of same name as script) in /home/
#    -- /input/ and /output/ subdirs within project dir

##  Set username
USER=aubaxh002

## Set project name
PROJ=psmc

## Create a directory on /scratch
mkdir /scratch/${USER}_${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}_${PROJ}/

##  Copy input files to scratch
cp /home/$USER/$PROJ/input/* /scratch/${USER}_${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/


## --------------------------------
## Load modules
# module load bwa/0.7.12
module load samtools/1.11


## --------------------------------
## Index reference genomes
# bwa index ./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa
# bwa index ./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa
	

## --------------------------------
## Align reads to reference genomes
## D genome - ordii
# while read -a line
# 	do
# 	bwa mem -t 20 -M \
# 	./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa \
# 	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
# 	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
# 	> ./align_files/${line[0]}_Dgenome.bam

## C genome	- ordii 
# 	bwa mem -t 20 -M \
# 	./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa \
# 	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
# 	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
# 	> ./align_files/${line[0]}_Cgenome.bam
# 	done < d_ordii_sra_list.txt

	
## D genome - stephensi
# while read -a line
# 	do
# 	bwa mem -t 20 -M \
# 	./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa \
# 	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
# 	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
# 	> ./align_files/${line[0]}_Dgenome.bam

## C genome	- stephensi 
# 	bwa mem -t 20 -M \
# 	./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa \
# 	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
# 	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
# 	> ./align_files/${line[0]}_Cgenome.bam
# 	done < d_stephensi_sra_list.txt
	
## --------------------------------
## Sort, merge, and sort BAM files
cd ./align_files/

## D genome - ordii
# while read -a line
# 	do
# 	samtools sort -@ 19 -o sorted_${line[0]}_Dgenome.bam ${line[0]}_Dgenome.bam
# 	done < ../d_ordii_sra_list.txt
# 
# samtools merge -@ 19 \
# d_ordii_Dgenome.bam \
# sorted_SRR1646412_Dgenome.bam \
# sorted_SRR1646413_Dgenome.bam \
# sorted_SRR1646414_Dgenome.bam \
# sorted_SRR1646416_Dgenome.bam \
# sorted_SRR1646417_Dgenome.bam \
# sorted_SRR1646418_Dgenome.bam \
# sorted_SRR1646419_Dgenome.bam \
# sorted_SRR1646420_Dgenome.bam \
# sorted_SRR1646421_Dgenome.bam \
# sorted_SRR1646422_Dgenome.bam \
# sorted_SRR1646423_Dgenome.bam
# 
# samtools sort -@ 19 -o sorted_d_ordii_Dgenome.bam d_ordii_Dgenome.bam

## C genome - ordii
# while read -a line
# 	do
# 	samtools sort -@ 19 -o sorted_${line[0]}_Cgenome.bam ${line[0]}_Cgenome.bam
# 	done < ../d_ordii_sra_list.txt
# 
# samtools merge -@ 19 d_ordii_Cgenome.bam \
# SRR1646412_Cgenome.bam \
# SRR1646413_Cgenome.bam \
# SRR1646414_Cgenome.bam \
# SRR1646415_Cgenome.bam \
# SRR1646416_Cgenome.bam \
# SRR1646417_Cgenome.bam \
# SRR1646418_Cgenome.bam \
# SRR1646419_Cgenome.bam \
# SRR1646420_Cgenome.bam \
# SRR1646421_Cgenome.bam \
# SRR1646422_Cgenome.bam \
# SRR1646423_Cgenome.bam
# 
# samtools sort -@ 19 -o sorted_d_ordii_Cgenome.bam d_ordii_Cgenome.bam

## D genome - stephensi
samtools sort -@ 19 -o sorted_d_stephensi_Dgenome.bam SRR14572526_Dgenome.bam 

## C genome	- stephensi
samtools sort -@ 19 -o sorted_d_stephensi_Cgenome.bam SRR14572526_Cgenome.bam

## --------------------------------
## Generate VCF for each alignment
cd /scratch/aubaxh002_psmc/variant_files
module purge
module load bcftools

## D genome - ordii
bcftools mpileup --threads 20 \
-f ../sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa \
../align_files/sorted_d_ordii_Dgenome.bam | bcftools call --threads 20 -c \
--output-type v \
--output unfilt_d_ordii_Dgenome.vcf

## C genome - ordii
bcftools mpileup --threads 20 \
-f ../sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa \
../align_files/sorted_d_ordii_Cgenome.bam | bcftools call --threads 20 -c \
--output-type v \
--output unfilt_d_ordii_Cgenome.vcf

## D genome - stephensi
bcftools mpileup --threads 20 \
-f ../sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa \
../align_files/sorted_d_stephensi_Dgenome.bam | bcftools call --threads 20 -c \
--output-type v \
--output unfilt_d_stephensi_Dgenome.vcf

## C genome  - stephensi
bcftools mpileup  --threads 20 \
-f ../sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa \
../align_files/sorted_d_stephensi_Cgenome.bam | bcftools call --threads 20 -c \
--output-type v \
--output unfilt_d_stephensi_Cgenome.vcf

