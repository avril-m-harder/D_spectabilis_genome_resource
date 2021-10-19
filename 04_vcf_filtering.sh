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
# module load bcftools/1.13


## --------------------------------
## Index, get stats on sample VCF files
cd variant_files/

# bgzip unfilt_d_spectabilis_Dgenome.vcf
# bgzip unfilt_d_stephensi_Dgenome.vcf
# bgzip unfilt_d_ordii_Dgenome.vcf

# bcftools index unfilt_d_spectabilis_Dgenome.vcf.gz
# bcftools index unfilt_d_stephensi_Dgenome.vcf.gz
# bcftools index unfilt_d_ordii_Dgenome.vcf.gz

# bcftools stats unfilt_d_ordii_Dgenome.vcf.gz > unfilt_d_ordii_Dgenome.vcf.gz_STATS.txt
# bcftools stats unfilt_d_spectabilis_Dgenome.vcf.gz > unfilt_d_spectabilis_Dgenome.vcf.gz_STATS.txt
# bcftools stats unfilt_d_stephensi_Dgenome.vcf.gz > unfilt_d_stephensi_Dgenome.vcf.gz_STATS.txt

## this results in loss of INFO details. should do something else.
# bcftools merge unfilt_d_spectabilis_Dgenome.vcf.gz \
# unfilt_d_ordii_Dgenome.vcf.gz \
# unfilt_d_stephensi_Dgenome.vcf.gz \
# --threads 5 \
# -O z \
# -o unfilt_3_spp_Dgenome.vcf.gz


## --------------------------------
## Filter sample VCF files (merging before filtering results in loss of INFO details)
# module load vcftools/0.1.14
# module load samtools
# 
# for i in d_ordii \
# 		 d_stephensi \
# 		 d_spectabilis
# 	do
# 	vcftools --gzvcf unfilt_${i}_Dgenome.vcf.gz \
# 	--out filtered_${i}_Dgenome \
# 	--recode --recode-INFO-all \
# 	--remove-indels \
# 	--minQ 30 \
# 	--minDP 15
# 	done
# 
# module purge
# module load bcftools/1.13
# 
# for i in d_ordii \
# 		 d_stephensi \
# 		 d_spectabilis
# 	do
# 	bgzip filtered_${i}_Dgenome.recode.vcf
# 	bcftools stats filtered_${i}_Dgenome.recode.vcf.gz > filtered_${i}_STATS.txt
# 	done
# 	
# for i in d_ordii \
# 		 d_stephensi \
# 		 d_spectabilis
# 	do
# 	bcftools view -M 2 \
# 	--output-type z \
# 	--output biallelic_filtered_${i}_Dgenome.recode.vcf.gz \
# 	filtered_${i}_Dgenome.recode.vcf.gz
# 	
# 	bcftools stats biallelic_filtered_${i}_Dgenome.recode.vcf.gz > \
# 	biallelic_filtered_${i}_STATS.txt
# 	done

## --------------------------------
## Merge to filter based on missingness (none allowed)
# module purge
# module load bcftools/1.13
# 
# for i in d_ordii \
# 		 d_stephensi \
# 		 d_spectabilis
# 	do
# 	bcftools index biallelic_filtered_${i}_Dgenome.recode.vcf.gz
# 	done
# 
# bcftools merge \
# 	biallelic_filtered_d_ordii_Dgenome.recode.vcf.gz \
# 	biallelic_filtered_d_stephensi_Dgenome.recode.vcf.gz \
# 	biallelic_filtered_d_spectabilis_Dgenome.recode.vcf.gz \
# 	--info-rules DP:join \
# 	--no-index \
# 	--output-type z \
# 	--output unfiltered_3_spp.vcf.gz
# 
# module purge
# module load vcftools/0.1.14
# 	vcftools \
# 	--gzvcf unfiltered_3_spp.vcf.gz \
# 	--max-missing 1 \
# 	--out filtered_3_spp \
# 	--recode --recode-INFO-all 


## --------------------------------
## Split back in to sample VCFs and convert to FASTQ files for PSMC input
# module purge
# module load bcftools/1.13

# bgzip filtered_3_spp.recode.vcf
# 
# bcftools stats filtered_3_spp.recode.vcf.gz > \
# filtered_3_spp.recode_STATS.txt
# 
# bcftools +split -Oz -o ../final_psmc_input filtered_3_spp.recode.vcf.gz

## run stats on each file to make sure they look good
module purge
module load bcftools/1.13

cd ../final_psmc_input
# 
# for i in d_ordii \
# 		 d_stephensi \
# 		 d_spectabilis
# 	do
# 	bcftools stats final_${i}_3spp_split.vcf.gz > final_3spp_split_${i}_STATS.txt
# 	done

# gunzip -c final_d_ordii_3spp_split.vcf.gz > final_d_ordii_3spp_split.vcf
# gunzip -c final_d_stephensi_3spp_split.vcf.gz > final_d_stephensi_3spp_split.vcf
# gunzip -c final_d_spectabilis_3spp_split.vcf.gz > final_d_spectabilis_3spp_split.vcf

# vcfutils.pl vcf2fq -d 10 -Q 30 final_d_ordii_3spp_split.vcf | gzip > \
# d_ordii.fq.gz

vcfutils.pl vcf2fq -d 10 -Q 30 final_d_stephensi_3spp_split.vcf | gzip > \ 
d_stephensi.fq.gz 

# vcfutils.pl vcf2fq -d 10 -Q 30 final_d_spectabilis_3spp_split.vcf | gzip > \
# d_spectabilis.fq.gz


## --------------------------------
## Copy results back to project output directory (in home)
cp *.fq.gz /home/aubaxh002/psmc/output/

mail -s 'PSMC VCF filtering finished' avrilharder@gmail.com <<< 'bcftools finished'












