#!/bin/bash

# ============================================================
# SLURM Job Submission Script
# Submit this script with: sbatch submit_job.sh
# ============================================================

#SBATCH --job-name=alignment_pipeline
#SBATCH --array=1                    # Set to 1-N where N = number of lines in sample_list.txt
                                     # Example for 8 samples: --array=1-8
#SBATCH --ntasks=1
#SBATCH -t 0-12:00:00               # Walltime limit (d-hh:mm:ss). Increase for large genomes.
#SBATCH -N 1                        # Number of nodes
#SBATCH -c 8                        # Number of cores
#SBATCH --mem=32G                   # Memory per node
#SBATCH -q grp_mlynch11
#SBATCH --partition=general
#SBATCH -o slurm.%j.out             # STDOUT log file
#SBATCH -e slurm.%j.err             # STDERR log file
#SBATCH --mail-type=ALL             # Email on job start, end, and failure
#SBATCH --mail-user=YOUREMAIL@asu.edu

# Path to the alignment pipeline script on scratch
bash /scratch/YOURUSERNAME/YOUR_PROJECT_DIR/alignment_pipeline.sh
