#!/usr/bin/env bash
# 03_voila_tsv.sh -- export the LSV table (all LSVs, incl. de novo junctions).
# --show-all keeps non-significant LSVs so 04 can apply its own crypTE filters.
set -euo pipefail
cd "$(dirname "$0")"; source params.sh; source synapse_io.sh

SG="${BUILD_DIR}/splicegraph.sql"
VOILA=$(ls "${DPSI_DIR}"/*.deltapsi.voila 2>/dev/null | head -1 || true)
[ -s "$SG" ]    || { echo "no splicegraph.sql in $BUILD_DIR (run 01)"; exit 1; }
[ -s "$VOILA" ] || { echo "no .deltapsi.voila in $DPSI_DIR (run 02)"; exit 1; }

mkdir -p "$VOILA_DIR"
OUT="${VOILA_DIR}/${GRP1_NAME}_${GRP2_NAME}.deltapsi.tsv"
voila tsv "$SG" "$VOILA" -f "$OUT" --show-all -j "$THREADS"

echo "voila tsv -> $OUT"
echo "columns:"; grep -vE '^#' "$OUT" | head -1 | tr '\t' '\n' | nl

# --- upload the LSV table to Synapse voila/ folder (no-op unless USE_SYNAPSE=1) ---
# The Synapse entity is named TDPneg_TDPpos_LSV_deltapsi.tsv (syn78075015); this step
# writes TDPneg_TDPpos.deltapsi.tsv, so copy to the Synapse name to create a new VERSION.
if syn_enabled; then
  SYN_TSV="${VOILA_DIR}/${GRP1_NAME}_${GRP2_NAME}_LSV_deltapsi.tsv"
  cp -f "$OUT" "$SYN_TSV"
  syn_store "$SYN_TSV" "$SYN_VOILA_FOLDER"
fi
