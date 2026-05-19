#!/bin/bash

set -euo pipefail

# ---- 1. Args ----
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <sample_name> <R1.fastq.gz> <R2.fastq.gz>" >&2
  exit 1
fi
SAMPLE="$1"
R1="$2"
R2="$3"

for f in "$R1" "$R2"; do
  [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

# ---- 2. Paths (EDIT IF NEEDED) ----
STAR_INDEX=Homo_sapiens.GRCh38.113_star
GTF=Homo_sapiens.GRCh38.113.gtf.gz

# Containers (read-only, on apptainer allowed-paths)
SHARED=/n/app/containers/shared/nf-core/rnaseq/3.23.0
STAR_IMG=$(ls "${SHARED}"/*star-2.7*.img    2>/dev/null | head -1)
SAMTOOLS_IMG=$(ls "${SHARED}"/*samtools-1*.img 2>/dev/null | head -1)
[[ -f "$STAR_IMG"     ]] || { echo "ERROR: STAR image not found under $SHARED"     >&2; exit 1; }
[[ -f "$SAMTOOLS_IMG" ]] || { echo "ERROR: samtools image not found under $SHARED" >&2; exit 1; }

# Output dir per-sample (so multiple runs don't collide)
OUTDIR=runs/star_only/${SAMPLE}
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Resources — match your SLURM allocation
THREADS="${SLURM_CPUS_ON_NODE:-12}"

# ---- 3. Decompress GTF if needed (STAR wants plain text) ----
if [[ "$GTF" == *.gz ]]; then
  GTF_UNCOMPRESSED="${OUTDIR}/$(basename "${GTF%.gz}")"
  if [[ ! -f "$GTF_UNCOMPRESSED" ]]; then
    echo "Decompressing GTF..."
    gunzip -c "$GTF" > "$GTF_UNCOMPRESSED"
  fi
  GTF="$GTF_UNCOMPRESSED"
fi

# ---- 4. Sanity prints ----
echo "=== STAR alignment: ${SAMPLE} ==="
echo "Started:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "Host:       $(hostname)"
echo "Threads:    ${THREADS}"
echo "Output dir: ${OUTDIR}"
echo "R1:         ${R1}"
echo "R2:         ${R2}"
echo "Index:      ${STAR_INDEX}"
echo "GTF:        ${GTF}"
echo "STAR img:   ${STAR_IMG}"
echo "samtools:   ${SAMTOOLS_IMG}"
echo

# ---- 5. STAR alignment → sorted BAM directly ----
echo "=== Running STAR ==="
apptainer exec --cleanenv \
  --bind /n/scratch:/n/scratch \
  --bind /n/app:/n/app \
  "${STAR_IMG}" \
  STAR \
    --runThreadN          "${THREADS}" \
    --genomeDir           "${STAR_INDEX}" \
    --readFilesIn         "${R1}" "${R2}" \
    --readFilesCommand    zcat \
    --sjdbGTFfile         "${GTF}" \
    --outSAMtype          BAM SortedByCoordinate \
    --outSAMunmapped      Within \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin      20 \
    --alignIntronMax      1000000 \
    --alignMatesGapMax    1000000 \
    --alignSJoverhangMin  8 \
    --alignSJDBoverhangMin 1 \
    --twopassMode         Basic \
    --outFileNamePrefix   "${SAMPLE}." \
    --limitBAMsortRAM     30000000000

# STAR's output BAM
BAM="${OUTDIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
[[ -f "$BAM" ]] || { echo "ERROR: STAR output BAM not found" >&2; exit 1; }

# ---- 6. Index the BAM ----
echo
echo "=== Indexing BAM ==="
apptainer exec --cleanenv \
  --bind /n/scratch:/n/scratch \
  "${SAMTOOLS_IMG}" \
  samtools index -@ "${THREADS}" "${BAM}"

# ---- 7. Quick stats so you know it worked ----
echo
echo "=== flagstat ==="
apptainer exec --cleanenv \
  --bind /n/scratch:/n/scratch \
  "${SAMTOOLS_IMG}" \
  samtools flagstat -@ "${THREADS}" "${BAM}" | tee "${SAMPLE}.flagstat.txt"

echo
echo "=== Outputs ==="
ls -lh "${OUTDIR}"

echo
echo "=== Finished: $(date '+%Y-%m-%d %H:%M:%S') ==="
