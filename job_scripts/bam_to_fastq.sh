#!/bin/bash
#SBATCH -J samtools_bam_to_fastq
#SBATCH -c 10                               # Request one core
#SBATCH -t 0-0:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=80G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/samtools_bam_to_fastq_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/samtools_bam_to_fastq_%j.err                 # File to which STDERR will be written, including job ID (%j)
set -eux

module load samtools/1.15.1
module load pigz/2.3.4

# comm -1 -3 <(ls -1 fastq | grep _1.fastq.gz$ | sed s/_1.fastq.gz//g | sort | uniq) <(ls -1 raw | grep .bam$ | sed s/.bam//g | sort | uniq) |
#  parallel sbatch bam_to_fastq.sh {}

# for i in SRR8571937  SRR8571938 SRR8571939 SRR8571940 SRR8571941 SRR8571942 SRR8571944 SRR8571945 SRR8571947 SRR8571948 SRR8571949 SRR8571950 SRR8571951 SRR8571952
for i in $@
do
if [ ! -f "raw/${i}_sorted.bam" ]; then
  samtools sort -n -@ 10 -m 5G -o "raw/${i}_sorted.bam" "raw/${i}.bam"
fi
if [ ! -f "fastq/${i}_1.fastq.gz" ]; then
  samtools fastq -1 >(pigz -p 4 -c > "fastq/${i}_1.fastq.gz") \
    -2 >(pigz -p 4 -c > "fastq/${i}_2.fastq.gz") \
    -0 "fastq/${i}_weird.fastq" \
    -@ 4 "raw/${i}_sorted.bam"
fi
done

