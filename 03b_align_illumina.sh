#!/bin/bash
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
module load bwa/0.7.12
module load samtools/1.11


## --------------------------------
## Index reference genomes
samtools faidx ./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fna.gz
samtools faidx ./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fna.gz
	


## --------------------------------
## Align reads to reference genomes
## D genome - ordii
while read -a line
	do
	bwa mem -10 -M \
	./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fna.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	> ./align_files/${line[0]}_Dgenome.bam

## C genome	- ordii 
	bwa mem -10 -M \
	./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fna.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	> ./align_files/${line[0]}_Cgenome.bam
	done < d_ordii_sra_list.txt
	
## D genome - stephensi
while read -a line
	do
	bwa mem -10 -M \
	./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fna.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	> ./align_files/${line[0]}_Dgenome.bam

## C genome	- stephensi 
	bwa mem -10 -M \
	./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fna.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	> ./align_files/${line[0]}_Cgenome.bam
	done < d_stephensi_sra_list.txt
	
## --------------------------------
## Merge and sort BAM files
cd ./align_files/

## D genome - ordii
samtools merge -o d_ordii_Dgenome.bam \
SRR1646412_Dgenome.bam \
SRR1646413_Dgenome.bam \
SRR1646414_Dgenome.bam \
SRR1646415_Dgenome.bam \
SRR1646416_Dgenome.bam \
SRR1646417_Dgenome.bam \
SRR1646418_Dgenome.bam \
SRR1646419_Dgenome.bam \
SRR1646420_Dgenome.bam \
SRR1646421_Dgenome.bam \
SRR1646422_Dgenome.bam \
SRR1646423_Dgenome.bam

samtools sort -o sorted_d_ordii_Dgenome.bam d_ordii_Dgenome.bam

## C genome - ordii
samtools merge -o d_ordii_Cgenome.bam \
SRR1646412_Cgenome.bam \
SRR1646413_Cgenome.bam \
SRR1646414_Cgenome.bam \
SRR1646415_Cgenome.bam \
SRR1646416_Cgenome.bam \
SRR1646417_Cgenome.bam \
SRR1646418_Cgenome.bam \
SRR1646419_Cgenome.bam \
SRR1646420_Cgenome.bam \
SRR1646421_Cgenome.bam \
SRR1646422_Cgenome.bam \
SRR1646423_Cgenome.bam

samtools sort -o sorted_d_ordii_Cgenome.bam d_ordii_Cgenome.bam

## D genome - stephensi
samtools sort sorted_d_stephensi_Dgenome.bam SRR14572526_Dgenome.bam 

## C genome	- stephensi
samtools sort sorted_d_stephensi_Cgenome.bam SRR14572526_Cgenome.bam

## --------------------------------
## Generate VCF for each alignment
## D genome - ordii
samtools mpileup -uf ../sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fna.gz \
sorted_d_ordii_Dgenome.bam | bcftools call -c \
--output-type v \
--output ../variant_files/unfilt_d_ordii_Dgenome.vcf

## C genome - ordii
samtools mpileup -uf ../sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fna.gz \
sorted_d_ordii_Cgenome.bam | bcftools call -c \
--output-type v \
--output ../variant_files/unfilt_d_ordii_Cgenome.vcf

## D genome - stephensi
samtools mpileup -uf ../sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fna.gz \
sorted_d_stephensi_Dgenome.bam | bcftools call -c \
--output-type v \
--output ../variant_files/unfilt_d_stephensi_Dgenome.vcf

## C genome  - stephensi
samtools mpileup -uf ../sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fna.gz \
sorted_d_stephensi_Cgenome.bam | bcftools call -c \
--output-type v \
--output ../variant_files/unfilt_d_stephensi_Cgenome.vcf

