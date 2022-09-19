set -eux

UNPAIRED=$(mktemp).fastq.gz
for i in SRR8571937  SRR8571938 SRR8571939 SRR8571940 SRR8571941 SRR8571942 SRR8571944 SRR8571945 SRR8571947 SRR8571948 SRR8571949 SRR8571950 SRR8571951 SRR8571952
do
cat trimmed/${i}_1_unpaired.fastq.gz trimmed/${i}_2_unpaired.fastq.gz > ${UNPAIRED}
hisat2 -x hisat2/Homo_sapiens.GRCh38_index \
    -p 12 \
    -1 trimmed/${i}_1.fastq.gz -2 trimmed/${i}_2.fastq.gz \
    -U $UNPAIRED \
    -S hisat2/${i}.sam
done

