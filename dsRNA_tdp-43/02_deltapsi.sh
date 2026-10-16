#!/usr/bin/env bash
# 02_deltapsi.sh -- quantify dPSI between TDPneg (grp1) and TDPpos (grp2).
# dPSI = PSI(grp1) - PSI(grp2); positive = higher inclusion in TDPneg.
# Note: deltapsi treats the two arms as independent groups of 7; it does NOT
# model the patient pairing (leverage that downstream if needed).
set -euo pipefail
cd "$(dirname "$0")"; source params.sh
: "${MAJIQ_LICENSE_FILE:?set MAJIQ_LICENSE_FILE}"

# collect per-group .majiq files (built in step 01) from the sample table
mapfile -t TN < <(awk -F, 'NR>1 && $3=="TN"{print "'"$BUILD_DIR"'/"$1".majiq"}' "$SAMPLE_MAP")
mapfile -t TP < <(awk -F, 'NR>1 && $3=="TP"{print "'"$BUILD_DIR"'/"$1".majiq"}' "$SAMPLE_MAP")
echo "grp1 ${GRP1_NAME}: ${#TN[@]} files | grp2 ${GRP2_NAME}: ${#TP[@]} files"
for f in "${TN[@]}" "${TP[@]}"; do [ -s "$f" ] || { echo "missing $f (run 01_build.sh)"; exit 1; }; done

mkdir -p "$DPSI_DIR"
majiq deltapsi \
    -grp1 "${TN[@]}" \
    -grp2 "${TP[@]}" \
    -n "$GRP1_NAME" "$GRP2_NAME" \
    -o "$DPSI_DIR" \
    -j "$THREADS"

echo "deltapsi done -> $DPSI_DIR"
ls -la "$DPSI_DIR"
