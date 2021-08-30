#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - LARGE queue      |
#   |    - 20 CPU + 50 Gb   |
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

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/

## --------------------------------
## Load modules
module load samtools/1.11
module load bcftools/1.10.2
MINIMAP=/scratch/aubaxh002_minimap2/minimap2/minimap2


## --------------------------------
## Index reference genomes for Minimap2
# $MINIMAP -t 10 -x map-hifi \
# -d ./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.mmi \
# ./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa
# 
# $MINIMAP -t 10 -x map-hifi \
# -d ./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.mmi \
# ./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa


## --------------------------------
## Align D. spectabilis reads to D. spectabilis reference...
# mkdir align_files/
# 
# cat d_spectabilis_sra_list.txt | while read line
# 	do
# 	$MINIMAP -ax map-hifi -t 10 \
# 	./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.mmi \
# 	./sra_downloads/${line}_1.fastq.gz > ./align_files/${line}_Dgenome.sam
# 	
# 	samtools view -@ 9 -S -b ./align_files/${line}_Dgenome.sam > \
# 	./align_files/${line}_Dgenome.bam
# 	done

## ... and to C. canadensis reference
# cat d_spectabilis_sra_list.txt | while read line
# 	do
# 	$MINIMAP -ax map-hifi -t 10 \
# 	./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.mmi \
# 	./sra_downloads/${line}_1.fastq.gz > ./align_files/${line}_Cgenome.sam
# 	
# 	samtools view -@ 9 -S -b ./align_files/${line}_Cgenome.sam > \
# 	./align_files/${line}_Cgenome.bam
# 	done


## --------------------------------
## Merge and sort BAM files
# cd /scratch/aubaxh002_psmc/align_files/

## D genome alignment
# samtools merge -@ 9 d_spectabilis_Dgenome.bam \
# SRR14662548_Dgenome.bam \
# SRR14662549_Dgenome.bam \
# SRR14662550_Dgenome.bam \
# SRR14662551_Dgenome.bam \
# SRR14662552_Dgenome.bam
# 
# samtools sort -@ 9 -o sorted_d_spectabilis_Dgenome.bam d_spectabilis_Dgenome.bam

## C genome alignment
# samtools merge -@ 9 d_spectabilis_Cgenome.bam \
# SRR14662548_Cgenome.bam \
# SRR14662549_Cgenome.bam \
# SRR14662550_Cgenome.bam \
# SRR14662551_Cgenome.bam \
# SRR14662552_Cgenome.bam
# 
# samtools sort -@ 9 -o sorted_d_spectabilis_Cgenome.bam d_spectabilis_Cgenome.bam

## copy to home directory
# cp sorted_d_spectabilis_Dgenome.bam /home/aubaxh002/psmc/output/
# cp sorted_d_spectabilis_Cgenome.bam /home/aubaxh002/psmc/output/


## --------------------------------
## Index reference genomes for samtools/bcftools
cd /scratch/aubaxh002_psmc/variant_files

# samtools faidx ../sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa
# samtools faidx ../sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa


## --------------------------------
## Generate VCF for each alignment
## D genome alignment
bcftools mpileup --threads 20 \
-f ../sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa \
--output raw_d_spectabilis_Dgenome.vcf \
--output-type v \
../align_files/sorted_d_spectabilis_Dgenome.bam

bcftools call --threads 20 -c \
--output-type v \
--output unfilt_d_spectabilis_Dgenome.vcf \
raw_d_spectabilis_Dgenome.vcf


## C genome alignment
bcftools mpileup --threads 20 \
-f ../sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa \
--output raw_d_spectabilis_Cgenome.vcf \
--output-type v \
../align_files/sorted_d_spectabilis_Cgenome.bam

bcftools call --threads 20 -c \
--output-type v \
--output unfilt_d_spectabilis_Cgenome.vcf \
raw_d_spectabilis_Cgenome.vcf
