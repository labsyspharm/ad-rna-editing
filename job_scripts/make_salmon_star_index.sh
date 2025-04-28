# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

wget ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
wget ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

# Remove transcript versions from cDNA files

gzip -cd Homo_sapiens.GRCh38.cdna.all.fa.gz | \
  sed -E 's/(>ENST[0-9]+)\.[0-9]+/\1/' | \
  pigz -c -p 12 > Homo_sapiens.GRCh38.cdna.all.no_version.fa.gz

gzip -cd Homo_sapiens.GRCh38.ncrna.fa.gz | \
  sed -E 's/(>ENST[0-9]+)\.[0-9]+/\1/' | \
  pigz -c -p 12 > Homo_sapiens.GRCh38.ncrna.no_version.fa.gz

# Concatenate ncRNA and cDNA
cat Homo_sapiens.GRCh38.cdna.all.no_version.fa.gz \
  Homo_sapiens.GRCh38.ncrna.no_version.fa.gz > \
  Homo_sapiens.GRCh38.transcriptome.no_version.fa.gz

# Create list of all transcripts in fasta
gzip -cd Homo_sapiens.GRCh38.transcriptome.no_version.fa.gz | \
  sed -nE 's/>(ENST[0-9]+).*/\1/p' | sort | uniq > \
  transcripts.txt

# Make GTF file only containing those transcripts
gzip -cd Homo_sapiens.GRCh38.113.gtf.gz | \
  grep -Ff transcripts.txt |
  pigz -c -p 12 > Homo_sapiens.GRCh38.113.filtered.gtf.gz

# Salmon authors recommend including the whole genome as decoy sequences, which
# are taken into account for finding alignments but are not actually quantified

gzip -cd Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | \
  grep "^>" | \
  cut -d " " -f 1 | \
  sed -e 's/>//g' > decoys.txt

# Concatenate cDNA and whole genome
cat Homo_sapiens.GRCh38.transcriptome.no_version.fa.gz \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > \
  Homo_sapiens.GRCh38.gentrome.fa.gz

salmon index -t Homo_sapiens.GRCh38.gentrome.fa.gz \
  -d decoys.txt -p 12 \
  -i Homo_sapiens.GRCh38.gentrome_salmon

gunzip -cd Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa
gunzip -cd Homo_sapiens.GRCh38.113.gtf.gz > Homo_sapiens.GRCh38.113.gtf
mkdir Homo_sapiens.GRCh38.113_star
STAR --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir Homo_sapiens.GRCh38.113_star \
  --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf
