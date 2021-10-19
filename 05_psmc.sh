#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - LARGE queue      |
#   |    - 1 CPU + 60 Gb    |
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
module load psmc/2016-1-21
cd ./final_psmc_input


## --------------------------------
## Convert .fastq consensus files to .psmcfa, run PSMC, plot results
## See Robinson et al. 2021 for their parameter selection process

# for i in d_ordii \
# 		 d_spectabilis \
# 		 d_stephensi

for i in d_spectabilis \
 		 d_stephensi
	do
	## 'r' shouldn't have much of an impact because PSMC estimates it? I think?
	psmc -N25 -t10 -r5 -p "4+25*2+4+6" -o ${i}.psmc ${i}.psmcfa
	
	psmc_plot.pl \
	-g 1 \
	-T ${i} \
	${i} ${i}.psmc
	done



## --------------------------------
## Copy results back to project output directory (in home)

mail -s 'PSMC finished' avrilharder@gmail.com <<< 'PSMC finished'