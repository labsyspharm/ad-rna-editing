#!/usr/bin/env python3

import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

INPUTS = {
    "de":         "de_te_locus_14.csv.gz",
    "ann":        "liu_te_annotation.csv.gz",
    "C":          "cryptic_junctions_TE_annotated_14samples.csv.gz",
    "IR":         "retained_introns_TE_annotated_14samples.csv.gz",
    "as_pattern": "cryptic_junctions_AS_pattern_14samples.csv.gz",
}
OUT_FIG = "Fig_excess_TE_misplicing_refined_14samples.png"
OUT_TAB = "excess_TE_misplicing_refined_events_14samples.csv"
KEEP = {"Cassette", "Alt 5", "Alt 3", "AFE", "ALE"}   # TE-forming categories kept
PROB, DPSI = 0.90, 0.10                               # MAJIQ significance thresholds

def apply_figure_style():
    mpl.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8,
        "axes.titlesize": 8.5, "legend.fontsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})

# ---------------------------------------------------------------- load TE-locus expression DE
de = pd.read_csv(INPUTS["de"]).rename(columns={"id": "GeneID"})
ann = pd.read_csv(INPUTS["ann"])
de = de.merge(ann[["GeneID", "Chr", "Start", "End"]], on="GeneID", how="left")
de["chrom"] = de.Chr.astype(str).str.replace(r"^chr", "", regex=True)
de["nlfdr"] = -np.log10(np.clip(de.FDR, 1e-30, 1))

# ---------------------------------------------------------------- cryptic junctions: sig TE events, both directions
C = pd.read_csv(INPUTS["C"], low_memory=False); C["chrom"] = C.chrom.astype(str)
ap = pd.read_csv(INPUTS["as_pattern"], dtype={"chrom": str})   # chrom MUST be str or the join silently drops rows
C = C.merge(ap[["gene_id", "chrom", "junc_start", "junc_end", "as_pattern"]],
            on=["gene_id", "chrom", "junc_start", "junc_end"], how="left")
pcol = [c for c in C.columns if "prob" in c.lower()][0]
cte = C[(C.splice_site_in_TE == True) & (C[pcol] >= PROB) & (C.delta_TDPloss.abs() >= DPSI)].copy()
cte["sdir"] = np.where(cte.delta_TDPloss > 0, "up", "down")
loc = cte.te_locus.astype(str).str.extract(r"^([^:]+):(\d+)-(\d+)")       # TE locus coords for each event
cte["tl_chrom"] = loc[0].str.replace(r"^chr", "", regex=True)
cte["tl_start"] = pd.to_numeric(loc[1], errors="coerce")
cte["tl_end"]   = pd.to_numeric(loc[2], errors="coerce")
cte = cte.dropna(subset=["tl_start", "tl_end"]).astype({"tl_start": int, "tl_end": int})
cte_keep = cte[cte.as_pattern.isin(KEEP)].copy()                          # drop Complex + Multi Exon Spanning

# ---------------------------------------------------------------- retained introns: sig TE events, both directions
IR = pd.read_csv(INPUTS["IR"], low_memory=False); IR["chrom"] = IR.chrom.astype(str)
irp = [c for c in IR.columns if "prob" in c.lower()][0]
irs = IR[(IR.te_fraction > 0) & (IR[irp] >= PROB) & (IR.delta_TDPloss.abs() >= DPSI)].copy()
irs["sdir"] = np.where(irs.delta_TDPloss > 0, "up", "down")
irs["chrom"] = irs.chrom.str.replace(r"^chr", "", regex=True)

# ---------------------------------------------------------------- overlap each event with a TE-expression locus
# interval index over DE loci, per chromosome
DEI = {}
for c, g in de.dropna(subset=["Start", "End"]).groupby("chrom"):
    g = g.sort_values("Start")
    DEI[str(c)] = (g.Start.to_numpy(np.int64), g.End.to_numpy(np.int64), g.index.to_numpy(),
                   int((g.End - g.Start).max()) if len(g) else 0)

def de_hits(chrom, s, e):
    """Return [(de_row_index, overlap_bp), ...] for DE loci overlapping [s,e)."""
    d = DEI.get(str(chrom))
    if d is None:
        return []
    ss, ee, idx, ml = d
    lo = np.searchsorted(ss, s - ml, "left"); hi = np.searchsorted(ss, e, "right")
    return [(int(idx[i]), int(min(e, ee[i]) - max(s, ss[i]))) for i in range(lo, hi) if ee[i] > s and ss[i] < e]

# panel-C data (crypTE events only): dPSI vs expression of the max-overlap TE locus
pc = []
for r in cte_keep.itertuples(index=False):
    hits = de_hits(r.tl_chrom, int(r.tl_start), int(r.tl_end))
    if hits:
        i = max(hits, key=lambda h: h[1])[0]
        pc.append((de.loc[i, "logFC"], float(r.delta_TDPloss), bool(de.loc[i, "FDR"] < 0.05)))

# event-centric overlay points: one representative (max-overlap) locus per event, crypTE + Alt-Intron
pts = []
for r in cte_keep.itertuples(index=False):
    hits = de_hits(r.tl_chrom, int(r.tl_start), int(r.tl_end))
    if hits:
        i = max(hits, key=lambda h: h[1])[0]
        pts.append((de.loc[i, "logFC"], de.loc[i, "nlfdr"], r.sdir, "crypTE", bool(de.loc[i, "FDR"] < 0.05), de.loc[i, "logFC"] > 0))
for r in irs.itertuples(index=False):
    hits = de_hits(r.chrom, int(r.intron_start), int(r.intron_end))
    if hits:
        i = max(hits, key=lambda h: h[1])[0]
        pts.append((de.loc[i, "logFC"], de.loc[i, "nlfdr"], r.sdir, "IR", bool(de.loc[i, "FDR"] < 0.05), de.loc[i, "logFC"] > 0))
PTS = pd.DataFrame(pts, columns=["logFC", "nlfdr", "sdir", "type", "esig", "exprup"])
PTS.to_csv(OUT_TAB, index=False)

nup = int((PTS.sdir == "up").sum()); ndn = int((PTS.sdir == "down").sum())
up_c = [int((cte_keep.sdir == "up").sum()),   int((irs.sdir == "up").sum())]
dn_c = [int((cte_keep.sdir == "down").sum()), int((irs.sdir == "down").sum())]
excess = set(de.index[(de.FDR < 0.05) & (de.logFC > 0)]); dep = set(de.index[(de.FDR < 0.05) & (de.logFC < 0)])
up = de.loc[list(excess)]; dn = de.loc[list(dep)]
P = pd.DataFrame(pc, columns=["te_logFC", "dPSI", "esig"])
xlo = np.floor(P.te_logFC.min() - 0.3); xhi = np.ceil(P.te_logFC.max() + 0.3)
print(f"events overlaid: {len(PTS)} (up {nup}, down {ndn}) | expr-significant {int(PTS.esig.sum())} ({100*PTS.esig.mean():.0f}%)")
print(f"keep crypTE: up {up_c[0]} down {dn_c[0]} | Alt-Intron: up {up_c[1]} down {dn_c[1]}")

# ---------------------------------------------------------------- 3-panel figure
apply_figure_style()
NS, RED, BLU, UPc, DNc = "#d9d9d9", "#c0392b", "#2c6fbb", "#d94801", "#238b8b"
fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.3))

ax = axes[0]
gsub = de.sample(min(40000, len(de)), random_state=0)
ax.scatter(gsub.logFC, gsub.nlfdr, s=2, c=NS, alpha=0.4, lw=0, rasterized=True)
ax.scatter(dn.logFC, dn.nlfdr, s=5, c=BLU, alpha=0.4, lw=0, label=f"expr down ({len(dep)})")
ax.scatter(up.logFC, up.nlfdr, s=5, c=RED, alpha=0.4, lw=0, label=f"expr up = excess TE ({len(excess)})")
pu = PTS[PTS.sdir == "up"]; pd_ = PTS[PTS.sdir == "down"]
ax.scatter(pu.logFC, pu.nlfdr, s=34, facecolors="none", edgecolors=UPc, lw=1.2, label=f"mis-spliced UP on loss ({nup} events)")
ax.scatter(pd_.logFC, pd_.nlfdr, s=34, facecolors="none", edgecolors=DNc, lw=1.2, label=f"mis-spliced DOWN on loss ({ndn} events)")
ax.axhline(-np.log10(0.05), ls="--", lw=0.7, c="#888"); ax.axvline(0, lw=0.5, c="#bbb"); ax.set_xlim(-8, 8)
ax.set_xlabel("TE-locus expression log$_2$FC (TDP-43 loss)"); ax.set_ylabel("$-$log$_{10}$FDR")
ax.set_title("A  TE-locus expression + mis-splicing events", fontsize=8.4)
ax.legend(loc="upper left", fontsize=5.7, markerscale=1.2, handletextpad=0.3, labelspacing=0.3)

ax = axes[1]; x = np.arange(2); w = 0.38
ax.bar(x - w/2, up_c, w, color=UPc, label="up on TDP-43 loss (gained)")
ax.bar(x + w/2, dn_c, w, color=DNc, label="down on loss (lost)")
for i in range(2):
    ax.text(i - w/2, up_c[i] + 2, up_c[i], ha="center", fontsize=6.6, fontweight="bold", color=UPc)
    ax.text(i + w/2, dn_c[i] + 2, dn_c[i], ha="center", fontsize=6.6, fontweight="bold", color=DNc)
ax.set_xticks(x); ax.set_xticklabels(["crypTE\n(Cassette/Alt/AFE/ALE)", "Alternative Intron"], fontsize=6.8)
ax.set_ylabel("significant TE-forming events"); ax.set_ylim(0, max(up_c + dn_c) * 1.25)
ax.set_title("B  Refined mis-splicing set (Complex + MES removed)", fontsize=8.4); ax.legend(loc="upper right", fontsize=6.2)

ax = axes[2]; ns_ = P[~P.esig]; sg_ = P[P.esig]
ax.scatter(ns_.te_logFC, ns_.dPSI, s=20, c=NS, edgecolors="#999", lw=0.4, label="expr n.s.")
ax.scatter(sg_.te_logFC, sg_.dPSI, s=28, c=RED, lw=0, label="expr sig (FDR<0.05)")
ax.axhline(0, lw=0.6, c="#bbb"); ax.axvline(0, lw=0.5, c="#bbb")
ax.set_xlabel("TE-locus expression log$_2$FC"); ax.set_ylabel("crypTE splicing ΔΨ (TDP-43 loss)")
ax.set_xlim(xlo, xhi); ax.set_ylim(-1, 1)
ax.text(0.97, 0.97, "↑ gained on loss", transform=ax.transAxes, ha="right", va="top", fontsize=5.8, color=UPc)
ax.text(0.97, 0.03, "↓ lost on loss", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=DNc)
ax.set_title("C  Splicing ≠ expression (both directions)", fontsize=8.4); ax.legend(loc="lower left", fontsize=6.2)
ax.text(0.03, 0.5, f"n={len(P)} keep-cat\ncrypTE events\n{100*P.esig.mean():.0f}% expr-sig", transform=ax.transAxes, fontsize=5.8, va="center")

fig.suptitle("Mis-splicing on the excess-TE landscape — TE-forming events only (Complex + Multi Exon Spanning removed), both splicing directions", fontsize=9.2, y=1.02)
fig.tight_layout()
fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
plt.close(fig)
print("saved:", OUT_FIG, "and", OUT_TAB)
