# https://github.com/hamidghaedi/REdiscoverTE
# https://github.com/ucsffrancislab/REdiscoverTE

wget http://research-pub.gene.com/REdiscoverTEpaper/software/REdiscoverTE_1.0.1.tar.gz

zcat rollup_annotation/REdiscoverTE_whole_transcriptome_hg38.fa.xz > genome.fa

conda activate salmon
salmon index \
	-t genome.fa \
	--threads 15 \
	-i REdiscoverTE

find . -maxdepth 1 -type d -not -name . -exec echo tar -cvzf {}.tar.gz -C {} ./ \;

echo sample$'\t'quant_sf_path > quants.tsv
ls -1 quants/*/quant.sf | awk -F/ '{split($0, a, "/");print a[2]"\t"$0}' >> quants.tsv

# https://github.com/ucsffrancislab/REdiscoverTE/blob/master/rollup.R
Rscript REdiscoverTE/rollup.R --metadata=quants.tsv --datadir=REdiscoverTE/rollup_annotation --threads=3 --assembly=hg38 --outdir=rollups
