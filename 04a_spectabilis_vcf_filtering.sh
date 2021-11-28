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
## Filter sample VCF files 
module load vcftools/0.1.14
module load bcftools/1.13
# module load samtools

# vcftools --gzvcf unfilt_d_spectabilis_Dgenome.vcf.gz \
# --out filtered_d_spectabilis_Dgenome \
# --recode --recode-INFO-all \
# --remove-indels \
# --minQ 30 \
# --minDP 10

bgzip -c filtered_d_spectabilis_Dgenome.recode.vcf > \
filtered_d_spectabilis_Dgenome.recode.vcf.gz

bcftools index filtered_d_spectabilis_Dgenome.recode.vcf.gz


## --------------------------------
## Check VCF stats

vcftools --gzvcf filtered_d_spectabilis_Dgenome.recode.vcf.gz \
--depth --out filtered_d_spectabilis_Dgenome_depth

## --------------------------------
mail -s 'PSMC spectabilis VCF filtering finished' avrilharder@gmail.com <<< 'bcftools finished'












