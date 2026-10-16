#!/usr/bin/env bash
# run_all.sh -- build -> deltapsi -> voila tsv -> crypTE-Exon calls.
# Prereqs: 00_install_majiq.sh done, conda env active, params.sh edited,
#          MAJIQ_LICENSE_FILE + BAMs + REF_GFF3 in place. See README.md.
set -euo pipefail
cd "$(dirname "$0")"; source params.sh; source synapse_io.sh

bash 01_build.sh       # fetches liu_sample_groups.csv; stores the .sj bundle
bash 02_deltapsi.sh
bash 03_voila_tsv.sh   # stores the LSV deltapsi tsv

# --- fetch the RepeatMasker table from Synapse for the crypTE step ---
# (04_crypte_from_lsv.py is pure-compute / local-path by design, so run_all manages
#  its Synapse I/O. no-op unless USE_SYNAPSE=1.)
RMSK_USE="$RMSK_CSV"
if syn_enabled; then
  syn_get "$SYN_RMSK" "$(dirname "$RMSK_CSV")"
  RMSK_USE="$(dirname "$RMSK_CSV")/repeatmasker_raw.csv.gz"   # pandas reads .gz transparently
fi

python 04_crypte_from_lsv.py \
    --voila-tsv "${VOILA_DIR}/${GRP1_NAME}_${GRP2_NAME}.deltapsi.tsv" \
    --rmsk-csv  "$RMSK_USE" \
    --out-dir   "$CRYPTE_DIR" \
    --min-dpsi "$MIN_DPSI" --min-prob "$MIN_PROB" --slop "$TE_ENDPOINT_SLOP"

# --- upload crypTE-Exon outputs to Synapse outputs/ folder (no-op unless USE_SYNAPSE=1) ---
if syn_enabled; then
  syn_store "${CRYPTE_DIR}/crypte_exon_junctions.tsv.gz" "$SYN_OUT_FOLDER" "$SYN_VOILA_TSV"
  syn_store "${CRYPTE_DIR}/crypte_exon_candidates.tsv"   "$SYN_OUT_FOLDER" "$SYN_VOILA_TSV"
  syn_store "${CRYPTE_DIR}/crypte_summary.txt"           "$SYN_OUT_FOLDER" "$SYN_VOILA_TSV"
fi

echo "ALL DONE. crypTE-Exon outputs in: $CRYPTE_DIR"
