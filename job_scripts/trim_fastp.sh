#!/bin/bash
#SBATCH -J fastp
#SBATCH -c 4                             # Request one core
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)
#SBATCH -o trimmed/logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e trimmed/logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fastp

set -eux

PREFIX_PATH="$1"
PREFIX_FILE="${PREFIX_PATH##*/}"

fastp -i "${PREFIX_PATH}_1.fastq.gz" -I "${PREFIX_PATH}_2.fastq.gz" \
  -o "trimmed/${PREFIX_FILE}_1.fastq.gz" -O "trimmed/${PREFIX_FILE}_1_unpaired.fastq.gz" \
  --unpaired1 "trimmed/${PREFIX_FILE}_2.fastq.gz" --unpaired2 "trimmed/${PREFIX_FILE}_2_unpaired.fastq.gz" \
  --cut_right \
  --detect_adapter_for_pe \
  -w 4 \
  --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --json "trimmed/${PREFIX_FILE}_report.json" --html "trimmed/${PREFIX_FILE}_report.html" \
  --overrepresentation_analysis \
  --report_title "$PREFIX_FILE"
