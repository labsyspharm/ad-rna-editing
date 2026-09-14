#!/usr/bin/env bash
# 03_voila_tsv.sh -- export the LSV table (all LSVs, incl. de novo junctions).
# --show-all keeps non-significant LSVs so 04 can apply its own crypTE filters.
set -euo pipefail
cd "$(dirname "$0")"; source params.sh

SG="${BUILD_DIR}/splicegraph.sql"
VOILA=$(ls "${DPSI_DIR}"/*.deltapsi.voila 2>/dev/null | head -1 || true)
[ -s "$SG" ]    || { echo "no splicegraph.sql in $BUILD_DIR (run 01)"; exit 1; }
[ -s "$VOILA" ] || { echo "no .deltapsi.voila in $DPSI_DIR (run 02)"; exit 1; }

mkdir -p "$VOILA_DIR"
OUT="${VOILA_DIR}/${GRP1_NAME}_${GRP2_NAME}.deltapsi.tsv"
voila tsv "$SG" "$VOILA" -f "$OUT" --show-all -j "$THREADS"

echo "voila tsv -> $OUT"
echo "columns:"; grep -vE '^#' "$OUT" | head -1 | tr '\t' '\n' | nl
