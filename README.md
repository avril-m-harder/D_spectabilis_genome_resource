# D_spectabilis_genome_resource

### Scripts and files associated with genome resource article for *Dipodomys spectabilis*, currently available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.11.23.469702v1.full) and submitted for peer-review.

#### Scripts
* 01_ncbi_downloads.sh: Downloads HiFi reads for *D. spectabilis*, Illumina reads for *D. stephensi* and *D. ordii*, and reference assembly for *D. spectabilis*.

* 02a_align_hifi.sh: Aligns HiFi reads to the reference genome using Minimap 2 and generates VCF file.

* 02b_qc_illumina.sh: Runs FastQC on raw Illumina reads, Trimmomatic to clean and trim, then FastQC on cleaned reads. 

* 03b_align_illumina.sh: Aligns cleaned Illumina reads to reference genome and generates VCF files.

* 04a_\*_vcf_filtering.sh: 

* 04b_\*_psmc_settings_filtering.sh:

* 05a_\*_individual_psmc.sh:

* \*_snp_count.sh:

* minimap2_install.sh:


#### Files
* \*_sra_list.txt: 

* filtered_contigs.list: 
