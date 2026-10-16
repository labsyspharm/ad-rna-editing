#!/usr/bin/env python3
"""
Strand-aware AS-pattern classification of cryptic junctions  (Liu 14-sample cohort)
===================================================================================
Produces:
  * cryptic_junctions_AS_pattern_14samples.csv.gz   (one AS-pattern label per junction)

Assigns every cryptic junction a MAJIQ-Modulizer-style alternative-splicing category
(Cassette / Alt 5 / Alt 3 / AFE / ALE / Multi Exon Spanning / Complex) from its splice-site
geometry against the splicegraph's annotated + de-novo exon boundaries.

STRAND-AWARE (important): donor vs acceptor is decided by genomic geometry AND strand.
On the '+' strand the genomic-left splice site is the donor (Alt 5) and the right is the
acceptor (Alt 3); on the '-' strand these are swapped. The geometry tests are identical
for both strands; only the Alt 5 / Alt 3 label mapping flips on '-'.
(An earlier version decided Alt 5/Alt 3 by coordinate only and mislabeled minus-strand
junctions; this is the corrected classifier.)

INPUTS
  cryptic_junctions_TE_annotated_14samples.csv.gz   cryptic junctions (needs columns:
        gene_id, chrom, junc_start, junc_end, crypte_class_position, strand)   [artifact c1ef386e]
  <MAJIQ splicegraph>  results/build_off14/splicegraph.sql  (SQLite; the `exon` table with
        gene_id, start, end, annotated). This is a MAJIQ `build` output; regenerate with:
            majiq build <GFF3> -c config.ini -j 16 --disable-ir --min-experiments 1 \
                 --incremental -o results/build_off14
        reusing the per-sample .sj files (no re-quantification needed for the splicegraph).

USAGE:  python 01_crypte_AS_pattern_classify.py  [path/to/splicegraph.sql]
"""
import sys, sqlite3, pandas as pd
import synapse_io, syn_ids   # Synapse I/O (fetch inputs / store outputs)

C_PATH  = synapse_io.fetch(syn_ids.INPUTS["cryptic_junctions_TE_annotated_14samples.csv.gz"])  # Synapse input
SG_PATH = sys.argv[1] if len(sys.argv) > 1 else "results/build_off14/splicegraph.sql"  # LOCAL MAJIQ output (not in the plots Synapse folder)
OUT     = "cryptic_junctions_AS_pattern_14samples.csv.gz"

C = pd.read_csv(C_PATH, low_memory=False); C["chrom"] = C.chrom.astype(str)

# ---- splicegraph exon boundaries, per gene ----
con = sqlite3.connect(f"file:{SG_PATH}?mode=ro", uri=True)
allex = pd.read_sql("SELECT gene_id gid,start,end,annotated FROM exon", con); con.close()
annx = allex[(allex.annotated == 1) & (allex.start > 0) & (allex.end > 0)]
novx = allex[(allex.annotated == 0) & (allex.start > 0) & (allex.end > 0) & (allex.end > allex.start)]
ann_ends   = {g: set(v) for g, v in annx.groupby("gid").end}      # annotated exon ends   (genomic-left splice target)
ann_starts = {g: set(v) for g, v in annx.groupby("gid").start}    # annotated exon starts (genomic-right splice target)
nov_starts = {g: set(v) for g, v in novx.groupby("gid").start}    # de-novo exon starts
nov_ends   = {g: set(v) for g, v in novx.groupby("gid").end}      # de-novo exon ends
EMPTY = set()

def near(x, S):                       # tolerant membership (+/-2 bp)
    return any((x + d) in S for d in (0, -1, 1, -2, 2))

def classify_sa(gid, s, e, cpos, strand):
    if cpos == "TSS":                 # TE <=50 bp from a transcript 5' start
        return "AFE"
    if cpos == "ApA":                 # TE <=50 bp from a transcript 3' end
        return "ALE"
    ae = ann_ends.get(gid, EMPTY); ast = ann_starts.get(gid, EMPTY)
    ns = nov_starts.get(gid, EMPTY); ne = nov_ends.get(gid, EMPTY)
    if near(e, ns) or near(s, ne):    # junction into/out of a de-novo exon -> cassette
        return "Cassette"
    nd = not near(s, ae)              # genomic-LEFT site is novel (not an annotated exon end)
    na = not near(e, ast)            # genomic-RIGHT site is novel (not an annotated exon start)
    minus = strand in ("-", "C")
    if nd and not na:                 # left site novel: donor on '+', acceptor on '-'
        return "Alt 3" if minus else "Alt 5"
    if na and not nd:                 # right site novel: acceptor on '+', donor on '-'
        return "Alt 5" if minus else "Alt 3"
    if (not nd) and (not na):         # both sites annotated -> exon-skip across annotated exons
        return "Multi Exon Spanning"
    return "Complex"                  # both sites novel, not flanking a de-novo exon

C["as_pattern"] = [classify_sa(g, int(s), int(e), cp, st)
                   for g, s, e, cp, st in zip(C.gene_id, C.junc_start, C.junc_end, C.crypte_class_position, C.strand)]
C[["gene_id", "chrom", "junc_start", "junc_end", "crypte_class_position", "as_pattern"]].to_csv(OUT, index=False, compression="gzip")
print("saved", OUT, "| counts:", dict(C.as_pattern.value_counts()))

# upload the AS-pattern intermediate to Synapse (intermediates/)
synapse_io.store(OUT, syn_ids.FOLDERS["intermediates"],
                 used=[syn_ids.INPUTS["cryptic_junctions_TE_annotated_14samples.csv.gz"]],
                 executed="01_crypte_AS_pattern_classify.py")
