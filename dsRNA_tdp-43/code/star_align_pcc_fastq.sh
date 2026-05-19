#!/bin/bash
set -uo pipefail

# ---- 1. Paths ----
INPUT_DIR=fastq
BAM_DIR=bam
SAMPLE_CSV=pcc_with_fastq.csv
STAR_INDEX=Homo_sapiens.GRCh38.113_star
GTF=${BAM_DIR}/Homo_sapiens.GRCh38.113.gtf
STAR_IMG=89d-data.img
THREADS=16

# ---- 2. Sanity checks ----
[[ -d "$INPUT_DIR"  ]] || { echo "ERROR: INPUT_DIR missing: $INPUT_DIR"   >&2; exit 1; }
[[ -d "$BAM_DIR"    ]] || { echo "ERROR: BAM_DIR missing: $BAM_DIR"       >&2; exit 1; }
[[ -d "$STAR_INDEX" ]] || { echo "ERROR: STAR_INDEX missing: $STAR_INDEX" >&2; exit 1; }
[[ -f "$GTF"        ]] || { echo "ERROR: GTF missing: $GTF"               >&2; exit 1; }
[[ -f "$STAR_IMG"   ]] || { echo "ERROR: STAR image missing: $STAR_IMG"   >&2; exit 1; }
[[ -f "$SAMPLE_CSV" ]] || { echo "ERROR: SAMPLE_CSV missing: $SAMPLE_CSV" >&2; exit 1; }

cd "$BAM_DIR"

# ---- 3. Build sample list from CSV (specimenID stripped of -PCC) ----
SAMPLE_LIST="$(mktemp)"
tail -n +2 "$SAMPLE_CSV" \
	  | awk -F, '$1 != "NA" && $1 != "" { sub("-PCC$","",$1); print $1 }' \
	    > "$SAMPLE_LIST"
	    n_samples=$(wc -l < "$SAMPLE_LIST")

	    # ---- 4. Helpers ----
	    log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
	    warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $*" >&2; }
	    bam_exists() {
		      local sample="$1"
		        compgen -G "${BAM_DIR}/${sample}*Aligned.sortedByCoord.out.bam" > /dev/null
		}

	# ---- 5. Header ----
	log "================================================================"
	log "START — STAR batch alignment (pcc_with_fastq.csv)"
	log "Host             : $(hostname)"
	log "SLURM_JOB_ID     : ${SLURM_JOB_ID:-N/A}"
	log "Threads          : ${THREADS}"
	log "Samples queued   : ${n_samples}"
	log "INPUT_DIR        : ${INPUT_DIR}"
	log "BAM_DIR          : ${BAM_DIR}"
	log "SAMPLE_CSV       : ${SAMPLE_CSV}"
	log "STAR_INDEX       : ${STAR_INDEX}"
	log "GTF              : ${GTF}"
	log "STAR_IMG         : $(basename "${STAR_IMG}")"
	log "================================================================"

	# ---- 6. Main loop ----
	declare -i n_total=0 n_skip_bam=0 n_skip_r1=0 n_skip_r2=0 n_run_ok=0 n_run_fail=0

	while IFS= read -r SAMPLE; do
		  SAMPLE="${SAMPLE//$'\r'/}"
		    [[ -z "$SAMPLE" || "$SAMPLE" == "NA" ]] && continue
		      n_total+=1

		        # locate R1: matches "{SAMPLE}_S<digits>_R1_001.fastq.gz"
			  R1=$(compgen -G "${INPUT_DIR}/${SAMPLE}_S*_R1_001.fastq.gz" | head -1)
			    if [[ -z "$R1" ]]; then
				        warn "[SKIP-R1 ] ${SAMPLE} — no R1 matched ${SAMPLE}_S*_R1_001.fastq.gz"
					    n_skip_r1+=1
					        continue
						  fi

						    R2="${R1/_R1_001/_R2_001}"
						      if [[ ! -f "$R2" ]]; then
							          warn "[SKIP-R2 ] ${SAMPLE} — R2 not found: $(basename "$R2")"
								      n_skip_r2+=1
								          continue
									    fi

									      if bam_exists "$SAMPLE"; then
										          existing=$(ls "${BAM_DIR}/${SAMPLE}"*Aligned.sortedByCoord.out.bam 2>/dev/null | head -1)
											      log "[SKIP-BAM] ${SAMPLE} — exists: $(basename "$existing")"
											          n_skip_bam+=1
												      continue
												        fi

													  # ---- run STAR ----
													    log "[RUN     ] ${SAMPLE} — starting"
													      log "           R1=$(basename "$R1")"
													        log "           R2=$(basename "$R2")"

														  rm -rf "${SAMPLE}._STARtmp" "${SAMPLE}._STARgenome" "${SAMPLE}._STARpass1" 2>/dev/null || true
														    start_ts=$(date +%s)

														      apptainer exec --cleanenv \
															          --bind /n/scratch:/n/scratch \
																      --bind /n/app:/n/app \
																          --bind /tmp:/tmp \
																	      "${STAR_IMG}" \
																	          STAR \
																		        --runThreadN            "${THREADS}" \
																			      --genomeDir             "${STAR_INDEX}" \
																			            --readFilesIn           "${R1}" "${R2}" \
																				          --readFilesCommand      zcat \
																					        --sjdbGTFfile           "${GTF}" \
																						      --outSAMtype            BAM SortedByCoordinate \
																						            --outSAMunmapped        Within \
																							          --outFilterMultimapNmax 20 \
																								        --outFilterMismatchNmax 999 \
																									      --outFilterMismatchNoverReadLmax 0.04 \
																									            --alignIntronMin        20 \
																										          --alignIntronMax        1000000 \
																											        --alignMatesGapMax      1000000 \
																												      --alignSJoverhangMin    8 \
																												            --alignSJDBoverhangMin  1 \
																													          --twopassMode           Basic \
																														        --outTmpDir             "/tmp/star_${SAMPLE}_$$" \
																															      --outFileNamePrefix     "${SAMPLE}." \
																															            --limitBAMsortRAM       30000000000
														        rc=$?
															  end_ts=$(date +%s)
															    elapsed=$((end_ts - start_ts))

															      if [[ $rc -eq 0 ]] && bam_exists "$SAMPLE"; then
																          bam_size=$(du -h "${BAM_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam" 2>/dev/null | cut -f1)
																	      log "[OK      ] ${SAMPLE} — done in ${elapsed}s, BAM size=${bam_size}"
																	          rm -rf "${SAMPLE}._STARtmp" "${SAMPLE}._STARgenome" "${SAMPLE}._STARpass1" 2>/dev/null || true
																		      n_run_ok+=1
																		        else
																				    warn "[FAIL    ] ${SAMPLE} — STAR exit=$rc, elapsed=${elapsed}s"
																				        n_run_fail+=1
																					  fi
																				  done < "$SAMPLE_LIST"

																				  rm -f "$SAMPLE_LIST"

																				  # ---- 7. Summary ----
																				  log "================================================================"
																				  log "END — STAR batch summary"
																				  log "Total in list           : ${n_total}"
																				  log "Skipped (BAM exists)    : ${n_skip_bam}"
																				  log "Skipped (no R1 match)   : ${n_skip_r1}"
																				  log "Skipped (no R2)         : ${n_skip_r2}"
																				  log "Ran successfully        : ${n_run_ok}"
																				  log "Failed                  : ${n_run_fail}"
																				  log "================================================================"
