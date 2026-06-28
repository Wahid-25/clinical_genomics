Here's a rewritten version of your guide with the same technical content but different wording and improved formatting. You can use this in your project documentation, report, or GitHub README.

---

# Clinical Genomics Variant Calling Pipeline – Technical Documentation

## Pipeline Overview

This project implements an automated clinical genomics workflow using **Nextflow** to analyze DNA sequencing data. The pipeline processes raw sequencing reads by aligning them to a reference genome, generating an indexed alignment, and identifying genomic variants using a deep learning-based variant caller.

The workflow is fully containerized with **Apptainer/Singularity**, ensuring portability, reproducibility, and consistent execution across computing environments.

---

# Workflow Stages

## Step 1: Read Alignment and BAM Sorting (`ALIGN_AND_SORT`)

### Purpose

The first stage aligns raw sequencing reads to the reference genome and immediately sorts the resulting alignments into a BAM file.

### Software Used

* **minimap2** – Performs sequence alignment.
* **samtools sort** – Sorts the alignment records.

### Inputs

* `reads.fastq` – Raw sequencing reads.
* `ref_genome.fa` – Reference genome.
* `ref_genome.fa.fai` – FASTA index for rapid reference access.

### Output

* `aligned.bam` – Sorted alignment file.

### Command

```bash
minimap2 -a ref_genome.fa reads.fastq | samtools sort -o aligned.bam -
```

### Resource Allocation

* CPUs: **8**
* Memory: **24 GB**
* Container: **minimap2.sif**

---

## Step 2: BAM Index Generation (`INDEX_BAM`)

### Purpose

Creates an index for the BAM file, enabling rapid access to genomic regions during downstream analysis.

### Software Used

* **samtools index**

### Input

* `aligned.bam`

### Output

* `aligned.bam.bai`

### Command

```bash
samtools index aligned.bam
```

### Resource Allocation

* CPUs: **4**
* Memory: **8 GB**
* Container: **minimap2.sif**

---

## Step 3: Variant Calling (`CALL_VARIANTS`)

### Purpose

Detects genomic variants from the aligned sequencing reads using a deep learning approach specifically optimized for long-read sequencing.

### Software Used

* **Clair3**

### Inputs

* `aligned.bam`
* `aligned.bam.bai`
* `ref_genome.fa`
* `ref_genome.fa.fai`

### Sequencing Platform

* Oxford Nanopore Technology (**ONT**)

### Output

* `output.vcf`

### Command

```bash
run_clair3.sh \
--bam_fn=aligned.bam \
--ref_fn=ref_genome.fa \
--output=. \
--threads=8 \
--platform=ont
```

### Resource Allocation

* CPUs: **8**
* Memory: **24 GB**
* Container: **clair3.sif**

---

# Pipeline Workflow

```
reads.fastq
       +
ref_genome.fa
       +
ref_genome.fa.fai
          │
          ▼
ALIGN_AND_SORT
          │
          ▼
aligned.bam
          │
          ▼
INDEX_BAM
          │
          ▼
aligned.bam.bai
          │
          ▼
CALL_VARIANTS
          │
          ▼
output.vcf
```

---

# Project Directory Structure

## Input Data

**D:\clinical_genomics\data\**

```
reads.fastq
ref_genome.fa
ref_genome.fa.fai
```

These files provide the sequencing reads, the reference genome, and its corresponding index.

---

## Container Files

**D:\clinical_genomics\containers\**

```
minimap2.sif
minimap2.def
clair3.sif
clair3.def
```

The `.sif` files contain executable software environments, while the `.def` files define how each container was built.

---

## Output Directory

**D:\clinical_genomics\results\**

The pipeline automatically generates this directory during execution.

Expected outputs include:

```
aligned.bam
aligned.bam.bai
output.vcf
```

---

# Nextflow Configuration (`nextflow.config`)

The configuration file defines pipeline parameters and execution settings.

## Pipeline Parameters

| Parameter        | Description                        |
| ---------------- | ---------------------------------- |
| `reads`          | Path to the input FASTQ file       |
| `ref_genome`     | Reference genome FASTA file        |
| `ref_genome_fai` | FASTA index file                   |
| `outdir`         | Directory for storing output files |

### Execution Settings

```groovy
errorStrategy = 'retry'
maxRetries = 1
cleanup = false
apptainer.enabled = true
apptainer.autoMounts = true
```

### Configuration Details

* Automatically retries failed processes once.
* Retains intermediate working directories for debugging.
* Executes every process inside an Apptainer/Singularity container.
* Automatically mounts required host directories into containers.

---

# Resource Allocation

| Process        | Container    | CPUs | Memory |
| -------------- | ------------ | ---- | ------ |
| ALIGN_AND_SORT | minimap2.sif | 8    | 24 GB  |
| INDEX_BAM      | minimap2.sif | 4    | 8 GB   |
| CALL_VARIANTS  | clair3.sif   | 8    | 24 GB  |

---

# Executing the Workflow

Run the pipeline using:

```bash
nextflow run main.nf
```

### Execution Sequence

1. Raw FASTQ reads are aligned against the reference genome.
2. The resulting BAM file is sorted.
3. A BAM index is created.
4. Clair3 analyzes the indexed BAM file to detect genomic variants.
5. Final outputs are saved in the `results/` directory.

Each stage executes independently within its assigned container while utilizing the CPU and memory resources defined in the configuration.

---

# Output Files

After successful completion, the following files are available in the `results/` directory.

| File              | Description                                                         |
| ----------------- | ------------------------------------------------------------------- |
| `aligned.bam`     | Sorted sequence alignments                                          |
| `aligned.bam.bai` | BAM index used for rapid genomic access                             |
| `output.vcf`      | Variant Call Format (VCF) file containing detected genomic variants |

---

# Workflow Definition (`main.nf`)

The main Nextflow script coordinates three sequential processes:

* **ALIGN_AND_SORT** aligns sequencing reads to the reference genome and produces a sorted BAM file.
* **INDEX_BAM** generates an index for the BAM file to support efficient data retrieval.
* **CALL_VARIANTS** runs Clair3 on the aligned reads to identify genetic variants and produces a VCF output.

Nextflow channels automatically transfer outputs from one process to the next, enabling efficient workflow execution while maintaining modularity.

---

# Pipeline Highlights

* **Automated Workflow:** Executes all stages from alignment to variant detection without manual intervention.
* **Containerized Execution:** Each process runs inside its dedicated Apptainer/Singularity container.
* **Reproducible Analysis:** Fixed software environments and configuration ensure consistent results across systems.
* **Fault Tolerance:** Failed tasks are automatically retried once before termination.
* **Efficient Resource Management:** CPU cores and memory are allocated individually for each workflow stage.
* **Flexible Configuration:** Input files and output locations can be modified easily through `nextflow.config` without changing the workflow code.

This version conveys the same information as your original document while using different wording, clearer organization, and a more professional documentation style suitable for a technical report, thesis, or GitHub repository.
