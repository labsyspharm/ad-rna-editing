#!/usr/bin/env bash
# run_all.sh -- build -> deltapsi -> voila tsv -> crypTE-Exon calls.
# Prereqs: 00_install_majiq.sh done, conda env active, params.sh edited,
#          MAJIQ_LICENSE_FILE + BAMs + REF_GFF3 in place. See README.md.
set -euo pipefail
cd "$(dirname "$0")"; source params.sh

bash 01_build.sh
bash 02_deltapsi.sh
bash 03_voila_tsv.sh

python 04_crypte_from_lsv.py \
    --voila-tsv "${VOILA_DIR}/${GRP1_NAME}_${GRP2_NAME}.deltapsi.tsv" \
    --rmsk-csv  "$RMSK_CSV" \
    --out-dir   "$CRYPTE_DIR" \
    --min-dpsi "$MIN_DPSI" --min-prob "$MIN_PROB" --slop "$TE_ENDPOINT_SLOP"

echo "ALL DONE. crypTE-Exon outputs in: $CRYPTE_DIR"
