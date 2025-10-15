# cfDNA-donor-fraction

---

### Motivation  
This repository contains the source code and resources for **cfDNA-donor-fraction**, a reproducible bioinformatics pipeline for estimating the **donor-derived cell-free DNA (cfDNA)** fraction in recipient blood samples.  
The workflow is designed for **post-transplant monitoring** and implements a standardized genomic approach using a panel of 40 informative SNPs to quantify the proportion of donor cfDNA.  

The pipeline performs complete processing of sequencing data, including quality control, alignment, deduplication, recalibration, joint genotyping of donor and recipient, and computation of donor fraction based on allele-specific read counts.  
It integrates open-source tools from the **GATK**, **samtools**, and **Picard** ecosystems to ensure analytical reproducibility, traceability, and clinical-grade robustness.

---

### Code availability  
The pipeline is implemented as a fully automated **Bash** workflow (`cfDNA-donor-fraction.sh`) with embedded Python utilities for parsing variant tables and computing donor fractions.  
Key modules include:
- **Quality control:** *fastp*, *FastQC*, and *MultiQC* with integrated reports (RAW vs. FILTERED).  
- **Alignment and recalibration:** *BWA-MEM*, *samtools*, *Picard*, and *GATK BQSR*.  
- **Joint genotyping:** *GATK HaplotypeCaller* and *GenotypeGVCFs* for donor–recipient signatures.  
- **Donor fraction estimation:** *GATK Mutect2* and custom Python scripts to compute allele fractions across informative loci, including heterozygous anchor sites and dosage correction.  
- **Comprehensive QC metrics:** GC bias, base distribution, and coverage metrics integrated into a single MultiQC report.  

The pipeline is modular, portable, and compatible with HPC environments.  
It has been validated on controlled cfDNA mixtures to assess sensitivity and accuracy across donor fractions.

**License:** GNU General Public License v3.0 (GPLv3).  
**Version:** v9 (2025).

---

### Data availability  
This repository provides only example scripts and configuration files.  
Sequencing data from transplant samples are subject to data protection regulations and cannot be publicly released.  
Users can run the workflow with their own FASTQ files and reference resources (`hg38.fa`, `dbSNP`, and a custom 40-SNP VCF panel provided under `supporting_files/`).  

The output includes:
- Donor/recipient genotyping tables (`RD.signature.tsv`)  
- Donor–receptor variant maps (`donor_map.tsv`)  
- Informative anchor sites (`anchor_sites.tsv`)  
- Per-sample donor fraction tables (`*.per_site.tsv`)  
- Global donor percentage summary (`final_summary.tsv`)  
- Consolidated QC report (`multiqc_report.html`)

---

### Credits  
This pipeline was developed by **[Lucas-Ruiz, Fernando]**, Biomedical Research Institute of Murcia, Spain.  
Contact: **[your.email@institution.org]**  

The project builds upon open-source components from the **Broad Institute GATK**, **samtools**, **Picard**, and **MultiQC** frameworks.  
All analyses and automation logic were designed and implemented by the author to support cfDNA quantification in transplant and liquid biopsy studies.
