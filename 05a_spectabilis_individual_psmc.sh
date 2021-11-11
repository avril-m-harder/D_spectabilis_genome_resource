#!/bin/bash
#
#   +------------------------+
#   |  USE:                  |
#   |    - SMALL queue       |
#   |    - 1 CPU + def Gb    |
#   +------------------------+
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
# cp /home/$USER/$PROJ/input/* /scratch/${USER}_${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/


## --------------------------------
## Load modules 
module load bcftools/1.13
module load psmc/2016-1-21


## --------------------------------
## For each species, run PSMC on filtered VCF file (not filtered jointly 
## based on across-sample missingness)
##
## Convert .fastq consensus files to .psmcfa, run PSMC, plot results
## See Robinson et al. 2021 for their parameter selection process

cd ./final_psmc_input/

## convert VCF to FASTQ; -d set to 1/3 avg. depth, -D set to 2X avg. depth
vcfutils.pl vcf2fq -d 8 -D 50 \
../variant_files/filtered_d_spectabilis_Dgenome.recode.vcf | \
gzip > d_spectabilis.fq.gz

## convert FASTQ to PSMCFA
fq2psmcfa d_spectabilis.fq.gz > d_spectabilis.psmcfa
	
## run PSMC using selected settings
psmc -N25 -t10 -r5 -p "8+25*2+2+4" -o d_spectabilis_N25_t10_r5.psmc d_spectabilis.psmcfa
## plot the results
psmc_plot.pl -g 0.5 -u 2.2e-09 d_spectabilis_N25_t10_r5_0.5 d_spectabilis_N25_t10_r5.psmc
psmc_plot.pl -g 1 -u 2.2e-09 d_spectabilis_N25_t10_r5_1 d_spectabilis_N25_t10_r5.psmc

mail -a d_spectabilis_N25_t10_r5_0.5.eps -s 'initial PSMC spectabilis finished' \
avrilharder@gmail.com <<< 'PSMC finished'

## --------------------------------
## Bootstrapping - start with 50 rounds, see how long that takes
## splitfa to split chromosomes into shorter segments (default == 500,000)
splitfa d_spectabilis.psmcfa 100000 > split_d_spectabilis.psmcfa

## run 50 rounds of bootstrapping (set by -b option)
for round in {1..50}
	do
	psmc -N25 -t10 -r5 -p "8+25*2+2+4" -b -o d_spectabilis_${round}.psmc split_d_spectabilis.psmcfa
	done


## --------------------------------
## Copy results back to project output directory (in home)

mail -s 'PSMC boots spectabilis finished' avrilharder@gmail.com <<< 'PSMC finished'