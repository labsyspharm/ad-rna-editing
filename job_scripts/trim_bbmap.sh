#!/bin/bash
#SBATCH -J trim
#SBATCH -c 4                               # Request one core
#SBATCH -t 0-0:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=1G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

# comm -1 -3 <(ls -1 trimmed | grep _1.fastq.gz$ | sed s/_1.fastq.gz//g | sort | uniq) <(ls -1 fastq | grep _1.fastq.gz | sed s/_1.fastq.gz//g | sort | uniq) |
#  parallel sbatch trim_bbmap.sh {}

module load java

# USAGE
# sbatch trim.sh

set -eux

BBMAP_PATH=/home/ch305/software/bbmap

PREFIX_PATH="$1"
PREFIX_FILE="${PREFIX_PATH##*/}"

echo "$SLURM_JOB_ID" > "trimmed/${PREFIX_FILE}.lock"

${BBMAP_PATH}/bbduk.sh -Xmx1g in1="${PREFIX_PATH}_1.fastq.gz" in2="${PREFIX_PATH}_2.fastq.gz" \
  out1="trimmed/${PREFIX_FILE}_1.fastq.gz" out2="trimmed/${PREFIX_FILE}_2.fastq.gz" \
  ref="${BBMAP_PATH}/resources/adapters.fa" \
  ktrim=r k=23 mink=13 hdist=1 stats="trimmed/${PREFIX_FILE}_trim_stats.txt" tbo tpe \
  threads=4

rm "trimmed/${PREFIX_FILE}.lock"
echo "$SLURM_JOB_ID" > "trimmed/${PREFIX_FILE}.done"
