set -eux
TRIMMOMATIC_PATH=/home/ch305/software/trimmomatic/0.39

parallel -j 6 --bar \
    java -Xmx40g -jar ${TRIMMOMATIC_PATH}/trimmomatic-0.39.jar PE -threads 20 -phred33 \
        -summary trimmed/{}_summary.txt \
        {}_1.fastq.gz {}_2.fastq.gz \
        trimmed/{}_1.fastq.gz trimmed/{}_1_unpaired.fastq.gz \
        trimmed/{}_2.fastq.gz trimmed/{}_2_unpaired.fastq.gz \
        ILLUMINACLIP:${TRIMMOMATIC_PATH}/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        ::: SRR8571951 SRR8571952
#        ::: SRR8571937  SRR8571938 SRR8571939 SRR8571940 SRR8571941 SRR8571942 SRR8571944 SRR8571945 SRR8571947 SRR8571948 SRR8571949 SRR8571950 SRR8571951 SRR8571952
