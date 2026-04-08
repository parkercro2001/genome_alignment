# Geiler-Samerotte Lab — Genome Alignment and Variant Calling Pipeline

This pipeline takes paired-end Illumina FASTQ files, aligns them to a reference genome, and calls SNPs and indels using the GATK4 best practices workflow. It is designed to run on the **ASU Sol HPC cluster**.

---

## Repository Contents

| File | Description |
|---|---|
| `alignment_pipeline.sh` | Main pipeline script — aligns reads and calls variants |
| `submit_job.sh` | SLURM job submission script |
| `sample_list.txt` | One sample name per line — drives the job array |

---

## Before You Begin

### 1. ASU Sol account and scratch space
You must have an active ASU Research Computing account and access to the Lynch lab allocation (`grp_mlynch11`). Contact Kara if you need to be added.

All analyses should be run in your scratch directory, which is located at:
```
/scratch/YOURASURITE/
```
Scratch storage is not backed up — copy important output files to your home directory or lab storage when done.

### 2. Copy pipeline files to scratch
Copy the three files in this repository to your project directory on scratch:
```bash
scp alignment_pipeline.sh submit_job.sh sample_list.txt \
    YOURASURITE@sol.asu.edu:/scratch/YOURASURITE/YOUR_PROJECT_DIR/
```

### 3. Upload your reference genome
Place your reference genome FASTA file in your project directory on scratch. Reference genomes are available from:
- **S. cerevisiae (sacCer3):** UCSC Genome Browser or SGD
- **Human (GRCh38):** Ensembl or NCBI
- **Mouse (GRCm39):** Ensembl or NCBI

### 4. Upload your trimmed FASTQ files
Create a `trimmed/` subdirectory inside your project directory and place your trimmed, paired-end FASTQ files there:
```
YOUR_PROJECT_DIR/
├── trimmed/
│   ├── SAMPLENAME_R1_001.fastq.gz
│   ├── SAMPLENAME_R2_001.fastq.gz
│   └── ...
├── alignment_pipeline.sh
├── submit_job.sh
├── sample_list.txt
└── your_reference_genome.fa
```

---

## Input File Requirements

### FASTQ naming convention
Your FASTQ files **must** follow this exact naming format:

```
SAMPLENAME_R1_001.fastq.gz
SAMPLENAME_R2_001.fastq.gz
```

- `SAMPLENAME` is the unique identifier for each sample (e.g., `SG_UDP0258V2_S1`)
- `_R1_001` and `_R2_001` indicate read 1 and read 2 of the paired-end reads
- Files must be gzip-compressed (`.fastq.gz`)

The sample name used in the filename must match exactly what is listed in `sample_list.txt`.

### Sample list format
`sample_list.txt` should contain one sample name per line, with no header:
```
SG_UDP0258V2_S1
SG_UDP0259V2_S2
SG_UDP0260V2_S3
```

Each line corresponds to one SLURM array task. Line 1 → task 1, line 2 → task 2, etc.

---

## Configuration

Before submitting, edit the configuration section at the top of **both** scripts.

### In `alignment_pipeline.sh`
Set the three variables in the CONFIGURATION block:
```bash
REF_GENOME="saccer3.fa"      # Filename of your reference genome (must be in BASE_DIR)
BASE_DIR="/scratch/YOURUSERNAME/YOUR_PROJECT_DIR"
TRIMMED_DIR="${BASE_DIR}/trimmed"
SAMPLE_LIST="${BASE_DIR}/sample_list.txt"
```

Also update the `--mail-user` line with your ASU email.

### In `submit_job.sh`
Update:
- `--array=1` → `--array=1-N` where N is the number of samples in your sample list
- `--mail-user` with your ASU email
- The `BASE_DIR` path in the `bash` command at the bottom

---

## How to Run

1. SSH into Sol:
   ```bash
   ssh YOURASURITE@sol.asu.edu
   ```

2. Navigate to your project directory:
   ```bash
   cd /scratch/YOURASURITE/YOUR_PROJECT_DIR
   ```

3. Submit the job array:
   ```bash
   sbatch submit_job.sh
   ```

4. Monitor job status:
   ```bash
   squeue -u YOURASURITE
   ```

5. Check for errors if a job fails:
   ```bash
   cat slurm.JOBID.err
   ```

Each sample runs as a separate job in the array, so all samples are processed in parallel. Output files for each sample are written to a subdirectory named after the sample.

---

## Output Files

For each sample, a subdirectory is created containing:

| File | Description |
|---|---|
| `SAMPLE_aligned_reads.sam` | Raw aligned reads |
| `SAMPLE_sorted_dedup_reads.bam` | Sorted BAM with duplicates marked |
| `SAMPLE_dedup_metrics.txt` | Duplication rate statistics |
| `SAMPLE_alignment_metrics.txt` | Alignment summary (% mapped, etc.) |
| `SAMPLE_insert_metrics.txt` | Insert size distribution statistics |
| `SAMPLE_insert_size_histogram.pdf` | Insert size histogram plot |
| `SAMPLE_depth_out.txt` | Per-base coverage depth |
| `SAMPLE_recal_reads.bam` | BQSR-corrected BAM (primary analysis input) |
| `SAMPLE_recalibration_plots.pdf` | Before/after BQSR quality plots |
| `SAMPLE_raw_variants_recal.vcf` | Raw variants called from recalibrated BAM |
| `SAMPLE_filtered_snps_final.vcf` | Final filtered SNPs |
| `SAMPLE_filtered_indels_final.vcf` | Final filtered indels |

The files you will use most often for downstream analysis are:
- `SAMPLE_filtered_snps_final.vcf` — final SNP calls
- `SAMPLE_filtered_indels_final.vcf` — final indel calls
- `SAMPLE_recal_reads.bam` — recalibrated alignments

---

## Pipeline Steps Explained

### Step 1 — Index the reference genome
BWA and GATK require index files to efficiently look up positions in the reference genome. This only needs to run once per reference genome. If you are running multiple samples against the same reference, the index will be rebuilt each time — this is harmless but slow. If you run many samples, consider running the indexing step separately once first.

### Step 2 — Align reads to the reference (BWA-MEM)
Maps paired-end reads to the reference genome. The output is a SAM file with each read's position in the genome. The read group (`@RG`) tag is added here and is required by GATK for downstream processing.

### Step 3 — Mark duplicates and sort (GATK MarkDuplicatesSpark)
PCR amplification during library preparation can create copies of the same DNA fragment. These duplicates inflate variant allele frequencies and must be flagged. Reads are also coordinate-sorted, which is required for BAM indexing and variant calling.

### Step 4 — Collect QC metrics
Generates quality control statistics about the alignment, including the percentage of reads that mapped, average coverage depth, and insert size distribution. Review these before proceeding — low alignment rates or unexpected insert sizes can indicate a problem with the sample or reference genome.

### Step 5 — Initial variant calling (GATK HaplotypeCaller)
Identifies candidate SNPs and indels across the genome. This first pass is not the final variant set — it is used to build a training set for base quality score recalibration in the next steps.

### Step 6 — Hard filter variants (for BQSR training set)
Applies quality score thresholds to the raw variants to produce a high-confidence subset. These filters are intentionally permissive at this stage — they just need to remove the worst artifacts so BQSR has clean training data.

### Step 7 — Base Quality Score Recalibration (BQSR)
The quality scores assigned by the sequencer to each base call are often systematically biased. BQSR builds a model of these biases using the known-variant sites from Step 6, then corrects the quality scores in the BAM file. `AnalyzeCovariates` generates a PDF showing the before/after comparison.

### Step 8 — Final variant calling (GATK HaplotypeCaller)
Runs HaplotypeCaller again on the recalibrated BAM. This is the primary output of the pipeline.

### Step 9 — Hard filter final variants
Applies stricter quality thresholds to the post-BQSR variants. These are your final SNP and indel calls for downstream analysis.

---

## Reference

Pipeline based on the GATK4 variant calling best practices:
https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
