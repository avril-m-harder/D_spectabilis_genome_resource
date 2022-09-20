# D_spectabilis_genome_resource

### Scripts and files associated with genome resource article for *Dipodomys spectabilis*, currently available at [Genome Biology and Evolution](https://academic.oup.com/gbe/article/14/1/evac005/6506520): "High-quality reference genome for an arid-adapted mammal, the banner-tailed kangaroo rat (<i>Dipodomys spectabilis</i>)"

#### Scripts
* 01_ncbi_downloads.sh: Downloads HiFi reads for *D. spectabilis*, Illumina reads for *D. stephensi* and *D. ordii*, and reference assembly for *D. spectabilis*.

* 02a_align_hifi.sh: Aligns HiFi reads to the reference genome using Minimap 2 and generates VCF file.

* 02b_qc_illumina.sh: Runs FastQC on raw Illumina reads, Trimmomatic to clean and trim, then FastQC on cleaned reads. 

* 03b_align_illumina.sh: Aligns cleaned Illumina reads to reference genome and generates VCF files.

* 04a_\*_vcf_filtering.sh: Quality filters VCF files to keep SNPs on specified contigs.

* 04b_\*_psmc_settings_filtering.sh: Filters VCF files using depth threshholds set for PSMC analyses.

* 05a_\*_individual_psmc.sh: Runs PSMC analysis + 50 bootstrap replicates.

* \*_snp_count.sh: Uses VCF files produced by \*_psmc_settings_filtering.sh to count heterozygous and homozygous genotypes.

* minimap2_install.sh: Downloads and installs Minimap 2 from Github.


#### Files
* \*_sra_list.txt: Accession numbers for SRA read downloads.

* filtered_contigs.list: List of contigs >= 100 kb, excluding contigs likely originating from sex chromosomes.
