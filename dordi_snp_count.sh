#!/bin/bash
#
cd /scratch/aubaxh002_psmc/variant_files/

## count total number of sites to make sure totals match up at the end
gzip -cd d_ordii_psmc_depth_filtered.vcf.gz | grep -v "^#" | wc -l


## count reference homozygous sites
gzip -cd d_ordii_psmc_depth_filtered.vcf.gz | grep "0/0" > dordi_hom_ref_sites.txt
wc -l dordi_hom_ref_sites.txt


## count alternate allele homozygous sites
gzip -cd d_ordii_psmc_depth_filtered.vcf.gz | grep "1/1" > dordi_alt_ref_sites.txt
wc -l dordi_alt_ref_sites.txt


## count heterozygous sites
gzip -cd d_ordii_psmc_depth_filtered.vcf.gz | grep "0/1" > dordi_het_sites.txt
wc -l dordi_het_sites.txt