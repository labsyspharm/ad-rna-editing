#!/bin/bash
#SBATCH -J rnaseq
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)
#SBATCH -o %x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %x_%j.err                 # File to which STDERR will be written, including job ID (%j)

set -eu
module load java/jdk-11.0.11

set -x

nextflow run nf-core/rnaseq -r 3.18.0 -profile singularity \
  -params-file nf-params.json -c ~/o2.config \
  -with-trace -with-report "report_${SLURM_JOB_ID}.html" \
  -with-timeline "timeline_${SLURM_JOB_ID}.html" \
  "$@"
