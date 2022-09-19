set -eux
for i in SRR8571937 SRR8571938 SRR8571939 SRR8571940 SRR8571941 SRR8571942 SRR8571944 SRR8571945 SRR8571947 SRR8571948 SRR8571949 SRR8571950 SRR8571951 SRR8571952
do
salmon quant -i Homo_sapiens.GRCh38.gentrome_including_variants_index -l A \
  --seqBias --gcBias --posBias --recoverOrphans -p 20 \
  -1 <(cat trimmed/${i}_1.fastq.gz trimmed/${i}_1_unpaired.fastq.gz) -2 <(cat trimmed/${i}_2.fastq.gz trimmed/${i}_2_unpaired.fastq.gz) -o quants/${i}
done
