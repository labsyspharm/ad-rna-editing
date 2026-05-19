#!/usr/bin/env bash
set -uo pipefail

# --- self-contained paths (no reliance on parent shell) ---
WORK_DIR=~/data
SRC_DIR="$WORK_DIR/input_files"
OUT_DIR="$WORK_DIR/fastq_from_bam"
TMP_DIR="$OUT_DIR/tmp"
SAMPLE_LIST="$WORK_DIR/pcc_bam_samples.txt"

cd "$WORK_DIR"
mkdir -p "$OUT_DIR" "$TMP_DIR" logs

# --- sanity checks ---
command -v samtools >/dev/null || { echo "ERROR: samtools not in PATH" >&2; exit 1; }
command -v pigz     >/dev/null || { echo "ERROR: pigz not in PATH"     >&2; exit 1; }
[ -d "$SRC_DIR" ]    || { echo "ERROR: $SRC_DIR missing"     >&2; exit 1; }
[ -f "$SAMPLE_LIST" ] || { echo "ERROR: $SAMPLE_LIST missing" >&2; exit 1; }

while IFS= read -r i; do
	  i="${i//$'\r'/}"                 # strip stray CR
	    [ -z "$i" ] || [ "$i" = "NA" ] && continue

	      echo "=== $(date '+%H:%M:%S') $i ==="

	        # locate BAM across the three naming variants in your folder
		  bam=""
		    for cand in \
			          "$SRC_DIR/${i}-PCC.final.bam" \
				        "$SRC_DIR/${i}-PCC.bam"       \
					      "$SRC_DIR/${i}.final.bam"     \
					            "$SRC_DIR/${i}.bam"
		      do
			          if [ -f "$cand" ]; then bam="$cand"; break; fi
				    done
				      if [ -z "$bam" ]; then
					          echo "  SKIP: no BAM for $i" >&2
						      continue
						        fi
							  echo "  BAM: $bam"

							    # resume-safe
							      if [ -s "$OUT_DIR/${i}_1.fastq.gz" ] && [ -s "$OUT_DIR/${i}_2.fastq.gz" ]; then
								          echo "  already done"; continue
									    fi

									      sorted="$TMP_DIR/${i}_sorted.bam"
									        samtools sort -n -@ 8 -m 4G -o "$sorted" "$bam"

										  samtools fastq \
											      -1 >(pigz -p 3 -c > "$OUT_DIR/${i}_1.fastq.gz") \
											          -2 >(pigz -p 3 -c > "$OUT_DIR/${i}_2.fastq.gz") \
												      -0 >(pigz -p 3 -c > "$OUT_DIR/${i}_unpaired.fastq.gz") \
												          -s >(pigz -p 3 -c > "$OUT_DIR/${i}_singleton.fastq.gz") \
													      -@ 4 "$sorted"

										    rm -f "$sorted"
										      echo "  done $i"
									      done < "$SAMPLE_LIST" 2>&1 | tee -a logs/bam_to_fastq_interactive.log
