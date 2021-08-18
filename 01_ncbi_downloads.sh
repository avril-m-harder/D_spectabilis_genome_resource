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
module load sra


## --------------------------------
## NCBI downloads
# mkdir sra_downloads/
cd sra_downloads/

## D. spectabilis HiFi reads
cat ../d_spectabilis_sra_list.txt | while read line
	do
	fastq-dump --split-files --origfmt --gzip ${line}
	done

## D. ordii Illumina reads
cat ../d_ordii_sra_list.txt | while read line
	do
	fastq-dump --split-files --origfmt --gzip ${line}
	done

## D. stephensi Illumina reads
cat ../d_stephensi_sra_list.txt | while read line
	do
	fastq-dump --split-files --origfmt --gzip ${line}
	done
	
## C. canadensis genome assembly (GCA_001984765.1)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/984/765/GCA_001984765.1_C.can_genome_v1.0/GCA_001984765.1_C.can_genome_v1.0_genomic.fna.gz

gunzip -c GCA_001984765.1_C.can_genome_v1.0_genomic.fna.gz > \
GCA_001984765.1_C.can_genome_v1.0_genomic.fa

## D. spectabilis genome assembly (GCA_019054845.1)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/054/845/GCA_019054845.1_ASM1905484v1/GCA_019054845.1_ASM1905484v1_genomic.fna.gz

gunzip -c GCA_019054845.1_ASM1905484v1_genomic.fna.gz > \
GCA_019054845.1_ASM1905484v1_genomic.fa
