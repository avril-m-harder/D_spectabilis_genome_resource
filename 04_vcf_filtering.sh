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
PROJ=assem_stats

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
module load bcftools/1.10.2
module load samtools/1.11
module load vcftools/0.1.14


## --------------------------------
## Merge VCF files for same references for filtering
cd variant_files/
gzip *.vcf

bcftools merge unfilt_d_spectabilis_Dgenome.vcf.gz \
unfilt_d_ordii_Dgenome.vcf.gz \
unfilt_d_stephensi_Dgenome.vcf.gz \
--threads 5 \
-o unfilt_3_spp_Dgenome.vcf.gz

bcftools merge unfilt_d_spectabilis_Cgenome.vcf.gz \
unfilt_d_ordii_Cgenome.vcf.gz \
unfilt_d_stephensi_Cgenome.vcf.gz \
--threads 5 \
-o unfilt_3_spp_Cgenome.vcf.gz

## --------------------------------
## Filter merged VCF file -- not using this for now, vcftools instead
# bcftools filter \
# --SnpGap 3 \
# --threads 5 \
# --exclude FORMAT/DP[0-2] < 10 \
# -o filtered_3_spp_Dgenome.vcf.gz \
# unfilt_3_spp_Dgenome.vcf.gz

## --max-missing 1 ==> no missing data allowed
## D genome
vcftools --vcf unfilt_3_spp_Dgenome.vcf.gz \
--out filtered_3_spp_Dgenome \
--recode --recode-INFO-all \
--remove-indels \
--max-missing 1 \
--minQ 30 \
--minDP 10 \
--site-mean-depth

## C genome
vcftools --vcf unfilt_3_spp_Cgenome.vcf.gz \
--out filtered_3_spp_Cgenome \
--recode --recode-INFO-all \
--remove-indels \
--max-missing 1 \
--minQ 30 \
--minDP 10 \
--site-mean-depth


## --------------------------------
## Split back in to sample VCFs and convert to FASTQ files for PSMC input






# vcfutils.pl vcf2fq -Q 30 [VCF] > cns.fq


## --------------------------------
## Copy results back to project output directory (in home)

