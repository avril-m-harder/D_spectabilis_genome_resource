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
module purge
module load bcftools/1.13

cd ./variant_files/

## --------------------------------
## Merge to filter based on missingness (none allowed)
bcftools merge \
biallelic_filtered_d_ordii_Dgenome.recode.vcf.gz \
biallelic_filtered_d_stephensi_Dgenome.recode.vcf.gz \
biallelic_filtered_d_spectabilis_Dgenome.recode.vcf.gz \
--info-rules DP:join \
--no-index \
--output-type z \
--output unfiltered_3_spp.vcf.gz

module purge
module load vcftools/0.1.14

vcftools \
--gzvcf unfiltered_3_spp.vcf.gz \
--max-missing 1 \
--out filtered_3_spp \
--recode --recode-INFO-all
	

## --------------------------------
## Get list of loci retained after missingness filtering
##### CHECK FORMATTING OF VCF, GET TAB-SEP CHROM/POS LIST OF SITES TO KEEP
##### FOR FILTERING WITH vcftools --positions

## filename = loci_no_missing_data.txt


## --------------------------------
## Filter sample VCFs based on sites filtered for missingness and convert to FASTQ files 
## for PSMC input

# for i in d_ordii \
# 		 d_spectabilis \
# 		 d_stephensi
# 	do
# 	vcftools \
# 	--gzvcf biallelic_filtered_${i}_Dgenome.recode.vcf.gz \
# 	--positions loci_no_missing_data.txt \
# 	--out missing_filtered_${i} \
# 	--recode --recode-INFO-all
# 	
# 	mv missing_filtered_${i}.recode.vcf.gz ../final_psmc_input/
# 	done
# 
# # run stats on each file to make sure they look good
# module purge
# module load bcftools/1.13
# 
# cd ../final_psmc_input
# 
# module purge
# module load bcftools/1.13
# 
# for i in d_ordii \
# 		 d_stephensi \
# 		 d_spectabilis
# 	do
# 	bcftools stats missing_filtered_${i}.recode.vcf.gz > missing_filtered_${i}_STATS.txt
# 	done
# 
# gunzip -c missing_filtered_d_ordii.recode.vcf.gz > final_d_ordii.vcf
# gunzip -c missing_filtered_d_stephensi.recode.vcf.gz > final_d_stephensi.vcf
# gunzip -c missing_filtered_d_spectabilis.recode.vcf.gz > final_d_spectabilis.vcf
# 
# vcfutils.pl vcf2fq final_d_ordii.vcf | gzip > \
# d_ordii.fq.gz
# 
# vcfutils.pl vcf2fq final_d_stephensi.vcf | gzip > \ 
# d_stephensi.fq.gz 
# 
# vcfutils.pl vcf2fq final_d_spectabilis.vcf | gzip > \
# d_spectabilis.fq.gz


## --------------------------------
## Copy results back to project output directory (in home)
# cp *.fq.gz /home/aubaxh002/psmc/output/

mail -s 'PSMC VCF filtering finished' avrilharder@gmail.com <<< 'bcftools finished'












