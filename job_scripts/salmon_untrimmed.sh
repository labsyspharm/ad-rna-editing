#!/bin/bash
#SBATCH -J salmon
#SBATCH -c 5                               # Request one core
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

#
# comm -1 -3 <(ls -1 quants | grep -e .done -e .lock | sed s/.done//g | sed s/.lock//g | sort | uniq) <(ls -1 fastq | grep _1.fastq.gz | sed s/_1.fastq.gz//g | sort | uniq) |
#  parallel sbatch salmon_untrimmed.sh {}

# USAGE EXAMPLE
# sbatch salmon_untrimmed.sh /n/scratch3/users/c/ch305/rna-editing/rna-seq/fastq/Homo_sapiens.GRCh38.gentrome_including_variants_index 962_TCX_2

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

# make a temporary lock file so I can see that this sample is
# currently being worked on
echo "$SLURM_JOB_ID" > "quants/${PREFIX_FILE}.lock"

salmon quant -i "$1" -l A \
  --seqBias --gcBias --posBias -p 5 \
  -1 "${PREFIX_PATH}_1.fastq.gz" \
  -2 "${PREFIX_PATH}_2.fastq.gz" \
  -o "quants/${PREFIX_FILE}"

# Marking sample as done on the file system
rm "quants/${PREFIX_FILE}.lock"
echo "$SLURM_JOB_ID" > "quants/${PREFIX_FILE}.done"
