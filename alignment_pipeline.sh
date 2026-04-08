#!/bin/bash

#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --mail-user=YOUREMAIL@asu.edu # Mail-to address

# ============================================================
# CONFIGURATION - Edit these variables before running
# ============================================================

# Path to your reference genome FASTA file.
# The file must be located in BASE_DIR.
# Examples by organism:
#   Yeast (S. cerevisiae S288C): saccer3.fa
#   Human (GRCh38):              GRCh38.fa
#   Mouse (GRCm39):              GRCm39.fa
REF_GENOME="saccer3.fa"

# Your working directory on scratch (where this script, your reference genome,
# your sample list, and your trimmed/ folder all live)
BASE_DIR="/scratch/YOURUSERNAME/YOUR_PROJECT_DIR"

# Directory containing your trimmed paired-end FASTQ files
TRIMMED_DIR="${BASE_DIR}/trimmed"

# Path to your sample list file
SAMPLE_LIST="${BASE_DIR}/sample_list.txt"

# ============================================================
# DO NOT EDIT BELOW THIS LINE
# ============================================================

module purge
module load mamba/latest
source activate base
module load samtools-1.17-xb
module load gatk-4.4.0.0-7u
module load bwa-0.7.17-2t
module load bowtie2-2.5.1-y6
module load openjdk-17.0.3_7-4s

# Read the sample name for this array task
bn=$(sed -n "${SLURM_ARRAY_TASK_ID}"p "$SAMPLE_LIST")
echo "Processing sample: $bn"

# Create output directory for this sample
out_dir="${BASE_DIR}/${bn}/"
mkdir -p "$out_dir"

# ============================================================
# STEP 1: Index the reference genome
# Only needs to run once per reference genome. GATK and BWA
# require these index files to efficiently look up sequences.
# ============================================================
bwa index "${BASE_DIR}/${REF_GENOME}"
gatk CreateSequenceDictionary -R "${BASE_DIR}/${REF_GENOME}"
samtools faidx "${BASE_DIR}/${REF_GENOME}"

# ============================================================
# STEP 2: Align reads to the reference genome
# Uses BWA-MEM to map paired-end reads to the reference.
# The @RG (read group) tag is required by GATK downstream.
# ============================================================
bwa mem \
    -K 100000000 \
    -Y \
    -R "@RG\tID:${bn}\tLB:${bn}\tPL:ILLUMINA\tPM:HISEQ\tSM:${bn}" \
    "${BASE_DIR}/${REF_GENOME}" \
    "${TRIMMED_DIR}/${bn}_R1_001.fastq.gz" \
    "${TRIMMED_DIR}/${bn}_R2_001.fastq.gz" \
    > "${out_dir}${bn}_aligned_reads.sam"

# ============================================================
# STEP 3: Mark duplicates and sort
# Identifies read pairs that are PCR duplicates (identical
# start/end positions). Duplicates are flagged rather than
# removed, and the BAM is coordinate-sorted.
# ============================================================
gatk MarkDuplicatesSpark \
    -I "${out_dir}${bn}_aligned_reads.sam" \
    -M "${out_dir}${bn}_dedup_metrics.txt" \
    -O "${out_dir}${bn}_sorted_dedup_reads.bam"

# ============================================================
# STEP 4: Collect alignment and insert size metrics
# Quality control stats: what fraction of reads aligned,
# and the distribution of insert sizes between read pairs.
# ============================================================
gatk CollectAlignmentSummaryMetrics \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -I "${out_dir}${bn}_sorted_dedup_reads.bam" \
    -O "${out_dir}${bn}_alignment_metrics.txt"

gatk CollectInsertSizeMetrics \
    -I "${out_dir}${bn}_sorted_dedup_reads.bam" \
    -O "${out_dir}${bn}_insert_metrics.txt" \
    -H "${out_dir}${bn}_insert_size_histogram.pdf"

samtools depth -a "${out_dir}${bn}_sorted_dedup_reads.bam" \
    > "${out_dir}${bn}_depth_out.txt"

# ============================================================
# STEP 5: Initial variant calling (pre-BQSR)
# HaplotypeCaller identifies SNPs and indels. This first pass
# is used to build a set of known variants for base quality
# score recalibration (BQSR) in the next steps.
# ============================================================
gatk HaplotypeCaller \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -I "${out_dir}${bn}_sorted_dedup_reads.bam" \
    -O "${out_dir}${bn}_raw_variants.vcf"

gatk SelectVariants \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -V "${out_dir}${bn}_raw_variants.vcf" \
    --select-type-to-include SNP \
    -O "${out_dir}${bn}_raw_snps.vcf"

gatk SelectVariants \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -V "${out_dir}${bn}_raw_variants.vcf" \
    --select-type-to-include INDEL \
    -O "${out_dir}${bn}_raw_indels.vcf"

# ============================================================
# STEP 6: Hard filter variants (for BQSR training set)
# Applies quality thresholds to the raw variants to produce
# a high-confidence set used as known sites for BQSR.
# ============================================================
gatk VariantFiltration \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -V "${out_dir}${bn}_raw_snps.vcf" \
    -O "${out_dir}${bn}_filtered_snps.vcf" \
    -filter-name "QD_filter"           -filter "QD < 2.0" \
    -filter-name "FS_filter"           -filter "FS > 60.0" \
    -filter-name "MQ_filter"           -filter "MQ < 40.0" \
    -filter-name "SOR_filter"          -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter"    -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

gatk VariantFiltration \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -V "${out_dir}${bn}_raw_indels.vcf" \
    -O "${out_dir}${bn}_filtered_indels.vcf" \
    -filter-name "QD_filter"  -filter "QD < 2.0" \
    -filter-name "FS_filter"  -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0"

gatk SelectVariants \
    --exclude-filtered \
    -V "${out_dir}${bn}_filtered_snps.vcf" \
    -O "${out_dir}${bn}_bqsr_snps.vcf"

gatk SelectVariants \
    --exclude-filtered \
    -V "${out_dir}${bn}_filtered_indels.vcf" \
    -O "${out_dir}${bn}_bqsr_indels.vcf"

# ============================================================
# STEP 7: Base Quality Score Recalibration (BQSR)
# Corrects systematic errors in the base quality scores
# assigned by the sequencer. Uses the filtered variants from
# Step 6 as known sites to distinguish true variants from
# sequencing errors during model training.
# ============================================================
gatk BaseRecalibrator \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -I "${out_dir}${bn}_sorted_dedup_reads.bam" \
    --known-sites "${out_dir}${bn}_bqsr_snps.vcf" \
    --known-sites "${out_dir}${bn}_bqsr_indels.vcf" \
    -O "${out_dir}${bn}_recal_data.table"

gatk ApplyBQSR \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -I "${out_dir}${bn}_sorted_dedup_reads.bam" \
    -bqsr "${out_dir}${bn}_recal_data.table" \
    -O "${out_dir}${bn}_recal_reads.bam"

gatk BaseRecalibrator \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -I "${out_dir}${bn}_recal_reads.bam" \
    --known-sites "${out_dir}${bn}_bqsr_snps.vcf" \
    --known-sites "${out_dir}${bn}_bqsr_indels.vcf" \
    -O "${out_dir}${bn}_post_recal_data.table"

gatk AnalyzeCovariates \
    -before "${out_dir}${bn}_recal_data.table" \
    -after  "${out_dir}${bn}_post_recal_data.table" \
    -plots  "${out_dir}${bn}_recalibration_plots.pdf"

# ============================================================
# STEP 8: Final variant calling (post-BQSR)
# Runs HaplotypeCaller again on the recalibrated BAM.
# This is the primary output used for downstream analysis.
# ============================================================
gatk HaplotypeCaller \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -I "${out_dir}${bn}_recal_reads.bam" \
    -O "${out_dir}${bn}_raw_variants_recal.vcf"

gatk SelectVariants \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -V "${out_dir}${bn}_raw_variants_recal.vcf" \
    --select-type-to-include SNP \
    -O "${out_dir}${bn}_raw_snps_recal.vcf"

gatk SelectVariants \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -V "${out_dir}${bn}_raw_variants_recal.vcf" \
    --select-type-to-include INDEL \
    -O "${out_dir}${bn}_raw_indels_recal.vcf"

# ============================================================
# STEP 9: Hard filter final variants
# Applies stricter thresholds to the post-BQSR variants.
# These are the final filtered variants for analysis.
# ============================================================
gatk VariantFiltration \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -V "${out_dir}${bn}_raw_snps_recal.vcf" \
    -O "${out_dir}${bn}_filtered_snps_final.vcf" \
    -filter-name "QD_filter"             -filter "QD < 8.0" \
    -filter-name "FS_filter"             -filter "FS > 55.0" \
    -filter-name "MQ_filter"             -filter "MQ < 45.0" \
    -filter-name "SOR_filter"            -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter"      -filter "MQRankSum < -5.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -3.0"

gatk VariantFiltration \
    -R "${BASE_DIR}/${REF_GENOME}" \
    -V "${out_dir}${bn}_raw_indels_recal.vcf" \
    -O "${out_dir}${bn}_filtered_indels_final.vcf" \
    -filter-name "QD_filter"  -filter "QD < 2.0" \
    -filter-name "FS_filter"  -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0"

echo "Done processing sample: $bn"
