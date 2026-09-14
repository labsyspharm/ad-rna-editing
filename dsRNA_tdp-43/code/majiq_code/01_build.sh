#!/usr/bin/env bash
# 01_build.sh -- generate the config and run `majiq build` (de novo ON).
set -euo pipefail
cd "$(dirname "$0")"; source params.sh
: "${MAJIQ_LICENSE_FILE:?set MAJIQ_LICENSE_FILE (see 00_install_majiq.sh)}"
[ -d "$BAMDIR" ] || { echo "BAMDIR not found: $BAMDIR"; exit 1; }

# --- optional: fetch GENCODE v47 GFF3 (chr-prefixed) if REF_GFF3 is missing ---
if [ ! -s "$REF_GFF3" ]; then
  echo "REF_GFF3 missing; fetching GENCODE v47 GFF3 ..."
  mkdir -p "$(dirname "$REF_GFF3")"
  url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz"
  curl -fSL "$url" -o "${REF_GFF3}.gz" && gunzip -f "${REF_GFF3}.gz"
fi

# --- contig-naming sanity check: BAM vs GFF3 must agree (chr1 vs 1) ---
bam1=$(ls "$BAMDIR"/*.bam | head -1)
bam_has_chr=$(samtools idxstats "$bam1" 2>/dev/null | grep -qE '^chr' && echo yes || echo no)
gff_has_chr=$(grep -m1 -vE '^#' "$REF_GFF3" | grep -qE '^chr' && echo yes || echo no)
echo "contig style -> BAM chr-prefixed: $bam_has_chr | GFF3 chr-prefixed: $gff_has_chr"
[ "$bam_has_chr" = "$gff_has_chr" ] || {
  echo "!! MISMATCH: BAM and GFF3 use different contig naming."
  echo "   Use an Ensembl GFF3 (no chr) if your BAMs are Ensembl-aligned, or a"
  echo "   chr-prefixed GFF3 (GENCODE) if chr-prefixed. Aborting."; exit 1; }

# --- generate the builder config from the 14-sample table ---
python make_config.py \
    --sample-map "$SAMPLE_MAP" --bamdir "$BAMDIR" \
    --genome "$GENOME" --strandness "$STRANDNESS" \
    --grp1-name "$GRP1_NAME" --grp2-name "$GRP2_NAME" \
    --out "${BUILD_DIR}/config.ini"

# --- build. de novo junction/exon detection is ON by default (what we want). ---
mkdir -p "$BUILD_DIR"
majiq build "$REF_GFF3" \
    -c "${BUILD_DIR}/config.ini" \
    -j "$THREADS" \
    -o "$BUILD_DIR"
    # keep defaults: de novo ON, IR ON. To match a junction-only analysis add
    # --disable-ir . To require more evidence for novel junctions tune
    # --min-denovo / --minreads / --minpos .

echo "build done -> $BUILD_DIR (per-sample .majiq + splicegraph.sql)"
ls -la "$BUILD_DIR"
