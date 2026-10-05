#!/usr/bin/env python3

import sys, sqlite3, pandas as pd

C_PATH  = "cryptic_junctions_TE_annotated_14samples.csv.gz"
SG_PATH = sys.argv[1] if len(sys.argv) > 1 else "results/build_off14/splicegraph.sql"
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
