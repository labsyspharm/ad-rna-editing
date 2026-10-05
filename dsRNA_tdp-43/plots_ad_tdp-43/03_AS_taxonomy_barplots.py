#!/usr/bin/env python3

import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch

FA_PATH = "full_AS_pattern_classification_14samples.csv"
OUT     = "Fig_AS_taxonomy_modulizer_names.png"

def apply_figure_style():
    mpl.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8,
        "axes.titlesize": 8.5, "legend.fontsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 7.5,
        "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
apply_figure_style()

fa = pd.read_csv(FA_PATH)

# significance rate = significant / QUANTIFIED (equal denominator; panel C x-label).
# Must use *_quant, NOT *_all: unquantified events (no PSI estimate) can never be significant,
# so dividing by *_all deflates the non-TE bars (TE_quant == TE_all, so TE bars are unaffected).
rt = fa[fa.pattern != "TOTAL"].copy()
rt["TE_rate"]    = 100.0 * rt["TE_sig"]    / rt["TE_quant"].replace(0, np.nan)
rt["nonTE_rate"] = 100.0 * rt["nonTE_sig"] / rt["nonTE_quant"].replace(0, np.nan)
rt = rt.set_index("pattern")[["TE_rate", "nonTE_rate"]]

# our internal category -> exact MAJIQ Modulizer name
MODU = {"Cassette / skipped exon": "Cassette", "Alt 5' splice site": "Alt 5", "Alt 3' splice site": "Alt 3",
        "Alt first exon / cryptic TSS": "AFE", "Alt last exon / cryptic ApA": "ALE",
        "Novel exon-skip (annotated sites)": "Multi Exon Spanning", "Complex (both sites novel)": "Complex",
        "Intron retention (RI)": "Alternative Intron"}
# potential to form a TE-derived feature (drives y-tick label colour)
FORM = {"Cassette": "forms", "AFE": "forms", "ALE": "forms", "Alt 5": "partial", "Alt 3": "partial",
        "Multi Exon Spanning": "no", "Complex": "no", "Alternative Intron": "intron"}
FCOL = {"forms": "#2e7d32", "partial": "#d17a00", "no": "#8a8a8a", "intron": "#2c6fbb"}

d = fa[fa.pattern != "TOTAL"].copy(); d["modu"] = d.pattern.map(MODU)
d = d.sort_values("total_all"); y = np.arange(len(d))
TEc, NTEc = "#e07b39", "#9aa7b4"
fig, axes = plt.subplots(1, 3, figsize=(15.6, 5.2))

# A -- all events
ax = axes[0]
ax.barh(y, d.TE_all, color=TEc, label="TE"); ax.barh(y, d.nonTE_all, left=d.TE_all, color=NTEc, label="non-TE")
for i, (t, te) in enumerate(zip(d.total_all, d.TE_all)):
    ax.text(t * 1.06, i, f"{t:,} ({100*te/t:.0f}% TE)", va="center", fontsize=6.2)
ax.set_yticks(y); ax.set_yticklabels(d.modu); ax.set_xscale("log")
ax.set_xlim(100, ax.get_xlim()[1] * 4); ax.set_xlabel("all cryptic events (log)")
ax.set_title("A  All events by Modulizer type", fontsize=8.6); ax.legend(loc="lower right", fontsize=6.5)

# B -- significant events
ax = axes[1]
ax.barh(y, d.TE_sig, color=TEc); ax.barh(y, d.nonTE_sig, left=d.TE_sig, color=NTEc)
for i, (t, te) in enumerate(zip(d.total_sig, d.TE_sig)):
    if t > 0:
        ax.text(t + 2, i, f"{t}" + (f" ({te} TE)" if te > 0 else ""), va="center", fontsize=6.2)
ax.set_yticks(y); ax.set_yticklabels([]); ax.set_xlabel("significant events (↑ on loss)")
ax.set_xlim(0, d.total_sig.max() * 1.4); ax.set_title("B  Significant events (count)", fontsize=8.6)

# C -- significance rate (equal, quantified denominator)
ax = axes[2]; h = 0.36
teR  = [rt.loc[p, "TE_rate"]    for p in d.pattern]
nteR = [rt.loc[p, "nonTE_rate"] for p in d.pattern]
ax.barh(y + h/2, teR, height=h, color=TEc, label="TE")
ax.barh(y - h/2, [0 if pd.isna(v) else v for v in nteR], height=h, color=NTEc, label="non-TE")
for i, (t, n) in enumerate(zip(teR, nteR)):
    if not pd.isna(t):
        ax.text(t + 0.05, i + h/2, f"{t:.2f}%", va="center", fontsize=5.6)
    ax.text((0.05 if pd.isna(n) else n + 0.05), i - h/2, ("(TE only)" if pd.isna(n) else f"{n:.2f}%"),
            va="center", fontsize=5.4, color=("#999" if pd.isna(n) else "#333"), style=("italic" if pd.isna(n) else "normal"))
ax.set_yticks(y); ax.set_yticklabels([]); ax.set_xlabel("% of QUANTIFIED events significant")
ax.set_xlim(0, 3.3); ax.set_title("C  Significance rate (equal denominator)", fontsize=8.6); ax.legend(loc="lower right", fontsize=6.5)

for lab in axes[0].get_yticklabels():                    # colour type labels by TE-forming potential
    lab.set_color(FCOL[FORM[lab.get_text()]]); lab.set_fontweight("bold")
key = [Patch(fc=FCOL["forms"],   label="forms a TE-derived exon (Cassette/AFE/ALE)"),
       Patch(fc=FCOL["partial"], label="partial — TE extends an exon edge (Alt 5/3)"),
       Patch(fc=FCOL["no"],      label="no new TE-exon (Multi Exon Spanning, Complex)"),
       Patch(fc=FCOL["intron"],  label="TE retained inside an intron (Alt Intron)")]
fig.legend(handles=key, loc="lower center", ncol=4, fontsize=6.6, frameon=False, bbox_to_anchor=(0.5, -0.04),
           title="y-axis label colour = potential to FORM a TE-derived feature", title_fontsize=7)
fig.suptitle("Cryptic-splicing taxonomy in MAJIQ Modulizer event-type names — Liu 14-sample cohort (min-experiments=1)", fontsize=9.6, y=1.0)
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(OUT, dpi=200, bbox_inches="tight"); plt.close(fig)
print("saved", OUT)
