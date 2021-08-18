#!/bin/bash
#
#   +-------------+
#   |  USE 5 CPU  |
#   +-------------+
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

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/


## --------------------------------
## Load modules
module load fastqc/0.11.9
module load trimmomatic/0.39

## --------------------------------
## Run FASTQC, see where we're at for the D. ordii and D. stephensi reads
## D. ordii
while read -a line
	do
	fastqc -t 5 -o fastqc_output/ \
	./sra_downloads/${line[0]}_1.fastq.gz \
	./sra_downloads/${line[0]}_2.fastq.gz
	done < d_ordii_sra_list.txt
	
## D. stephensi
while read -a line
	do
	fastqc -t 5 -o fastqc_output/ \
	./sra_downloads/${line[0]}_1.fastq.gz \
	./sra_downloads/${line[0]}_2.fastq.gz
	done < d_stephensi_sra_list.txt
	
exit
	

## --------------------------------
## Trimmomatic to clean and trim adapters from reads
## D. ordii
while read -a line
	do
	trimmomatic PE -phred33 -threads 5 \
	./sra_downloads/${line[0]}_1.fastq.gz \
	./sra_downloads/${line[0]}_2.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_unpaired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	./cleaned_reads/trimmed_unpaired_${line[0]}_2.fastq.gz \
	LEADING:20 TRAILING:20 MINLEN:30 \
	ILLUMINACLIP:[adapter_file.fa:2:40:10]
	done < d_ordii_sra_list.txt

## D. stephensi
while read -a line
	do
	trimmomatic PE -phred33 -threads 5 \
	./sra_downloads/${line[0]}_1.fastq.gz \
	./sra_downloads/${line[0]}_2.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_unpaired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	./cleaned_reads/trimmed_unpaired_${line[0]}_2.fastq.gz \
	LEADING:20 TRAILING:20 MINLEN:30 \
	ILLUMINACLIP:[adapter_file.fa:2:40:10]
	done < d_stephensi_sra_list.txt


## --------------------------------
## Run FASTQC on cleaned/trimmed read files
## D. ordii
while read -a line
	do
	fastqc -t 5 -o fastqc_output/ \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_unpaired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	./cleaned_reads/trimmed_unpaired_${line[0]}_2.fastq.gz \
	done < d_ordii_sra_list.txt
	
## D. stephensi
while read -a line
	do
	fastqc -t 5 -o fastqc_output/ \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_unpaired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	./cleaned_reads/trimmed_unpaired_${line[0]}_2.fastq.gz \
	done < d_stephensi_sra_list.txt


