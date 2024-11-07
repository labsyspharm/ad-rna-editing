#!/bin/bash
#SBATCH -J salmon
#SBATCH -c 5                               # Request one core
#SBATCH -t 0-0:40                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

# USAGE EXAMPLE
# sbatch salmon.sh /n/scratch3/users/c/ch305/rna-editing/rna-seq/fastq/Homo_sapiens.GRCh38.gentrome_including_variants_index 962_TCX_2

# use python provided by miniconda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon

# sets options how shell script is executed:
# e = exit immediately on error
# u = treat using undeclared variables as error
# x = print commands to log before running them
set -eux

PREFIX_PATH="$2"
PREFIX_FILE="${PREFIX_PATH##*/}"


# $1 contains path to salmon index
# $2 contains sample name
# trimming outputs two fastq files. one contains paired reads,
# the other reads whose pair got lost during trimming.
# I'm using `cat` to concatenate them back together for alignment with salmon

salmon quant -i "$1" -l A \
  --seqBias --gcBias -p 5 \
  -1 "${PREFIX_PATH}_R1_001.fastq.gz" \
  -2 "${PREFIX_PATH}_R2_001.fastq.gz" \
  -o "quants/${PREFIX_FILE}"

# fastq/azenta/C9Total-1 fastq/azenta/C9dsRIP-1 fastq/azenta/CntrlTotal-1 fastq/azenta/CntrldsRIP-1
# fastq/bauer/P516_MAvet18947_A07v1_C9Total1_S1 fastq/bauer/P516_MAvet18947_B07v1_CntrlTotal1_S2 fastq/bauer/P516_MAvet18947_C07v1_C9dsRIP1_S3 fastq/bauer/P516_MAvet18947_D07v1_CntrldsRIP1_S4 fastq/bauer/P516_MAvet18947_E07v1_C9Total2_S5 fastq/bauer/P516_MAvet18947_F07v1_CntrlTotal2_S6 fastq/bauer/P516_MAvet18947_G07v1_C9dsRIP2_S7 fastq/bauer/P516_MAvet18947_H07v1_CntrldsRIP2_S8
