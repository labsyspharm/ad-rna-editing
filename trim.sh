#!/bin/bash
#SBATCH -J trim
#SBATCH -c 4                               # Request one core
#SBATCH -t 0-1:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)
#SBATCH -o %x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %x_%j.err                 # File to which STDERR will be written, including job ID (%j)
module load java

set -eux

TRIMMOMATIC_PATH=/home/ch305/software/trimmomatic/0.39

PREFIX_PATH="$1"
PREFIX_FILE="${PREFIX_PATH##*/}"

touch "trimmed/${PREFIX_FILE}.lock"

java -Xmx20g -jar ${TRIMMOMATIC_PATH}/trimmomatic-0.39.jar PE \
  -threads 4 -phred33 \
  -summary "trimmed/${PREFIX_FILE}_summary.txt" \
  "${PREFIX_PATH}_1.fastq.gz" "${PREFIX_PATH}_2.fastq.gz" \
  "trimmed/${PREFIX_FILE}_1.fastq.gz" "trimmed/${PREFIX_FILE}_1_unpaired.fastq.gz" \
  "trimmed/${PREFIX_FILE}_2.fastq.gz" "trimmed/${PREFIX_FILE}_2_unpaired.fastq.gz" \
  ILLUMINACLIP:${TRIMMOMATIC_PATH}/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

rm "trimmed/${PREFIX_FILE}.lock"
touch "trimmed/${PREFIX_FILE}.done"
