#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - SMALL queue      |
#   |    - 5 CPU + 4 Gb     |
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
#    -- script to be run in /home/projectdir/scripts
#    -- project directory (of same name as script) in /home/
#    -- /input/ and /output/ subdirs within project dir

## cd into working scratch directory
#cd /scratch/aubaxh002_psmc/align_files


## --------------------------------
## Load module
#module load samtools


## --------------------------------
## Run flagstat to get info on alignment files
#samtools flagstat -@ 5 \
#sorted_d_spectabilis_Dgenome.bam -O tsv > \
#sorted_d_spectabilis_Dgenome_STATS.tsv

#samtools flagstat -@ 5 \
#sorted_d_spectabilis_Cgenome.bam -O tsv > \
#sorted_d_spectabilis_Cgenome_STATS.tsv


## --------------------------------
## Load module
module load bcftools


## --------------------------------
## Run bcftools stats on variant files
cd /scratch/aubaxh002_psmc/variant_files

bcftools stats --threads 5 \
raw_d_spectabilis_Dgenome.vcf > raw_d_spectabilis_Dgenome.txt

bcftools stats --threads 5 \
raw_d_spectabilis_Cgenome.vcf >	raw_d_spectabilis_Cgenome.txt
