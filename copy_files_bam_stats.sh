#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - SMALL queue      |
#   |    - 4 CPU + def Gb   |
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

##  Set username
USER=aubaxh002

## Set project name
PROJ=psmc

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/


## --------------------------------
## Load module
module load samtools/1.11

## --------------------------------
## Run flagstat to check final sorted BAM files for all species/genome combos
cd align_files/

# samtools flagstat -@ 3 --output-fmt tsv sorted_d_stephensi_Cgenome.bam > \
# /home/aubaxh002/psmc/input/sorted_d_stephensi_Cgenome_STATS.txt
# 
# samtools flagstat -@ 3 --output-fmt tsv sorted_d_stephensi_Dgenome.bam > \
# /home/aubaxh002/psmc/input/sorted_d_stephensi_Dgenome_STATS.txt
# 
# samtools flagstat -@ 3 --output-fmt tsv sorted_d_ordii_Cgenome.bam > \
# /home/aubaxh002/psmc/input/sorted_d_ordii_Cgenome_STATS.txt
# 
# samtools flagstat -@ 3 --output-fmt tsv sorted_d_ordii_Dgenome.bam > \
# /home/aubaxh002/psmc/input/sorted_d_ordii_Dgenome_STATS.txt
# 
# samtools flagstat -@ 3 --output-fmt tsv sorted_d_spectabilis_Cgenome.bam > \
# /home/aubaxh002/psmc/input/sorted_d_spectabilis_Cgenome_STATS.txt
# 
# samtools flagstat -@ 3 --output-fmt tsv sorted_d_spectabilis_Dgenome.bam > \
# /home/aubaxh002/psmc/input/sorted_d_spectabilis_Dgenome_STATS.txt


## --------------------------------
## Copy final sorted BAM files to home directory
cp sorted_d*.bam /home/aubaxh002/psmc/input/
