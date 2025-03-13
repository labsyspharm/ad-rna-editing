#!/bin/bash
#SBATCH -J sprint
#SBATCH -c 12                               # Request one core
#SBATCH -t 2-0:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

#  find ~/scratch/dsrna/SRP185789/output/fastq -name '*_1.fastq.gz' | sed -E 's/\.[^.]+(\.[^.]+)*$//' | sed 's/_1$//' | head -n 1| while read path; do sbatch sprint.sh "$path" ~/scratch/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.fa; done

# python 2.7
# conda create -n sprint python=2.7
# conda activate sprint
# python setup.py install
# conda install -c bioconda samtools==1.2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate sprint

filename="${1##*/}"

sprint main -1 "${1}_1.fastq" -2 "${1}_2.fastq" \
  -p 12 \
  "$2" "$filename" "$3" samtools

echo "done"
