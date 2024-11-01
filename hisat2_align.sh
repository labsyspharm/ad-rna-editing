#!/bin/bash
#SBATCH -J hisat2
#SBATCH -c 10                               # Request one core
#SBATCH -t 0-1:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=30G                         # Memory total in MiB (for all cores)
#SBATCH -o hisat2/logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hisat2/logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)


hisat-3n-build -p 10 --snp genome.snp --haplotype genome.haplotype \
 --ss Homo_sapiens.GRCh38.107.ss --exon Homo_sapiens.GRCh38.107.exon \
 --seed 42 --base-change A,G --repeat-index \
 Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38_index_3n

set -eu
module load gcc/9.2.0 hisat2/2.2.1 samtools/1.15.1

set -x

FASTQ_PREFIX_PATH="$2"
FASTQ_PREFIX="${FASTQ_PREFIX_PATH##*/}"

LOCKFILE="hisat2/${FASTQ_PREFIX}.lock"
touch "$LOCKFILE"

UNPAIRED="$(mktemp).fastq.gz"
cat "${FASTQ_PREFIX_PATH}_1_unpaired.fastq.gz" "${FASTQ_PREFIX_PATH}_2_unpaired.fastq.gz" > "$UNPAIRED"

hisat2 -x "$1" \
    -p 10 \
    -1 "${FASTQ_PREFIX_PATH}_1.fastq.gz" -2 "${FASTQ_PREFIX_PATH}_2.fastq.gz" \
    -U "$UNPAIRED" \
    -S "hisat2/${FASTQ_PREFIX}.sam"

samtools sort -m 2G -@ 10 -o "hisat2/${FASTQ_PREFIX}.bam" "hisat2/${FASTQ_PREFIX}.sam"
rm "hisat2/${FASTQ_PREFIX}.sam"

samtools index -@ 10 "hisat2/${FASTQ_PREFIX}.bam"

rm "$LOCKFILE"
touch "hisat2/${FASTQ_PREFIX}.done"
