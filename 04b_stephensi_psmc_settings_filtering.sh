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
module load bcftools/1.13


## --------------------------------
## Try filtering on depth again
cd ./variant_files/

bcftools view -i 'INFO/DP>=11 & INFO/DP<=68' \
--threads 20 \
--output-type z \
--output d_stephensi_psmc_depth_filtered.vcf.gz \
filtcontigs_d_stephensi.recode.vcf



mail -s 'PSMC stephensi VCF filtering finished' avrilharder@gmail.com <<< 'bcftools finished'












