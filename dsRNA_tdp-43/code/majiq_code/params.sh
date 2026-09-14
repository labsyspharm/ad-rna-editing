# ============================================================================
# params.sh  --  edit these, then run the numbered scripts in order.
# Source of truth for all paths/params. Every script does `source params.sh`.
# ============================================================================

# ---- INPUTS YOU MUST PROVIDE -----------------------------------------------
# 1) Directory holding the 14 Liu sorted BAMs + their .bam.bai indexes.
#    BAM basenames must match the Run IDs in liu_sample_groups.csv, i.e.
#    SRR8571937.bam, SRR8571937.bam.bai, ...  (see README for renaming/indexing)
export BAMDIR="/path/to/liu/bams"

# 2) MAJIQ academic license file. Obtain by accepting the license at
#    https://majiq.biociphers.org/  (download -> academic). MAJIQ will not run
#    without it. This is the ONE thing that cannot be scripted for you.
export MAJIQ_LICENSE_FILE="/path/to/majiq_license_academic_official.lic"

# 3) GFF3 annotation. MAJIQ needs GFF3 (not GTF). Its contig names MUST match
#    the BAM (chr1 vs 1). 01_build.sh can fetch GENCODE v47 (chr-prefixed) for
#    you; set to that path after fetching, or point at your own.
export REF_GFF3="${PWD}/refs/gencode.v47.annotation.gff3"

# ---- REFERENCE FOR THE crypTE STEP -----------------------------------------
# RepeatMasker table (the repeatmasker_raw.csv already in this project works).
export RMSK_CSV="${PWD}/refs/repeatmasker_raw.csv"

# ---- PARAMETERS ------------------------------------------------------------
export GENOME="hg38"            # label for VOILA UCSC links only
# Liu used the NuGEN Ovation RNA-Seq System V2 -> NON-directional library.
# So strandness is None (unstranded). CONFIRM empirically (see README §strand);
# a wrong value silently zeroes junction coverage.
export STRANDNESS="None"        # one of: None | forward | reverse
export THREADS=16

export SAMPLE_MAP="${PWD}/liu_sample_groups.csv"
export OUTROOT="${PWD}/results"
export BUILD_DIR="${OUTROOT}/build"
export DPSI_DIR="${OUTROOT}/deltapsi"
export VOILA_DIR="${OUTROOT}/voila"
export CRYPTE_DIR="${OUTROOT}/crypte"

# deltapsi contrast: grp1 - grp2, so positive dPSI = higher inclusion in TDPneg
# (the cryptic-splicing direction expected on TDP-43 loss).
export GRP1_NAME="TDPneg"
export GRP2_NAME="TDPpos"

# crypTE / significance thresholds (04_crypte_from_lsv.py)
export MIN_DPSI=0.10            # |E[dPSI]| threshold for "changing"
export MIN_PROB=0.90            # probability_changing threshold
export TE_ENDPOINT_SLOP=2       # bp tolerance when testing junction end in TE
