"""
02_crossmodal_4panels.py
-------------
2x2 cross-modal figures relating the molecular cryptic-exon call and the autopsy
TDP-43 stage.
Molecular TDP-43 positive is >=1 STMN2/UNC13A cryptic exon (== permissive class_low).
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


import data_prep as dp
import plot_utils as pu

OUTDIR = os.environ.get("ROSMAP_OUTDIR", "figures")
os.makedirs(OUTDIR, exist_ok=True)

INK = pu.INK; SC = pu.SC; C_MOL = pu.C_MOL; C_AUT = pu.C_AUT





def _labels(a, bars, fmt="{:.0f}", dy=1):
    for b in bars:
        h = b.get_height()
        a.annotate(fmt.format(h), (b.get_x() + b.get_width()/2, h), ha="center",
                   va="bottom", fontsize=9.5, xytext=(0, dy), textcoords="offset points", color=INK)






def _panel_count_dist(a, counts, xlabels, cut_x, neg_txt, pos_txt, pos_color,
                      xlabel, ylabel, title, ymax=300):
    
    x = np.arange(len(counts))
    bars = a.bar(x, counts, color=SC[:len(counts)], edgecolor="white", width=0.7)
    _labels(a, bars)
    a.axvline(cut_x, ls="--", lw=1.2, color="#888")
    a.text(cut_x/2, max(counts)*0.9, neg_txt, ha="center", fontsize=9, color="#555")
    a.text(cut_x + (len(counts)-1-cut_x)/2, max(counts)*0.9, pos_txt, ha="center",
           fontsize=9, color=pos_color, fontweight="bold")
    a.set_xticks(x); a.set_xticklabels(xlabels); a.set_xlabel(xlabel); a.set_ylabel(ylabel)
    a.set_title(title, loc="left"); a.set_ylim(0, ymax); a.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)


def _panel_concordance(a, ct, row_labels, col_labels, xlabel, ylabel, title):
    M = ct.values.astype(float); tot = M.sum()
    a.imshow(M, cmap="Blues", aspect="auto")
    for i in range(2):
        for j in range(2):
            a.text(j, i, f"{int(M[i, j])}\n{M[i, j]/tot*100:.0f}%", ha="center", va="center",
                   fontsize=15, fontweight="bold", color="white" if M[i, j] > M.max()*0.55 else INK)
    a.set_xticks([0, 1]); a.set_xticklabels(col_labels); a.set_yticks([0, 1]); a.set_yticklabels(row_labels)
    a.set_xlabel(xlabel); a.set_ylabel(ylabel)
    agree = (M[0, 0] + M[1, 1]) / tot * 100 if row_labels[0].startswith(("beyond", "molecular positive")) else (M[0, 1] + M[1, 0]) / tot * 100
    a.set_title(title.format(n=int(tot), agree=agree), loc="left")
    for s in a.spines.values():
        s.set_visible(False)


def _fig_grid():
    fig, ax = plt.subplots(2, 2, figsize=(14.5, 10.4))
    fig.subplots_adjust(hspace=0.42, wspace=0.24, top=0.9, bottom=0.12, left=0.075, right=0.975)
    return fig, ax


def molecular_4panel(df, panelD="dx", stem="ROSMAP_TDP43_figure_PCC_cryptic1plus",
                     title="Molecular TDP-43 (≥1 cryptic exon = positive) — ROSMAP PCC RNA-seq samples",
                     region_note="PCC samples only (n=546)."):
    pu.setup()
    df = df.copy(); mp = df.dropna(subset=["tdp_st4"])
    fig, ax = _fig_grid()

    # cryptic-exon cnt distribution
    counts = [int((df.n_expressed == k).sum()) for k in [0, 1, 2, 3]]
    _panel_count_dist(ax[0, 0], counts, ["0", "1", "2", "3"], 0.5,
                      f"negative\n(n={counts[0]})",
                      f"positive: ≥1 cryptic exon (n={sum(counts[1:])})", C_MOL,
                      "Cryptic exons detected (n_expressed)", "specimens (n)",
                      "A   Cryptic-exon count distribution")

    ### molecular ->> positive vs autopsy stage
    a = ax[0, 1]; x = np.arange(4)
    pos = [mp[mp.tdp_st4 == s]["molpos"].mean()*100 for s in [0, 1, 2, 3]]
    ns = [int((mp.tdp_st4 == s).sum()) for s in [0, 1, 2, 3]]
    _labels(a, a.bar(x, pos, color=C_MOL, edgecolor="white", width=0.6), "{:.0f}%")
    a.set_xticks(x); a.set_xticklabels([f"0\nnone\n(n={ns[0]})", f"1\namyg\n(n={ns[1]})",
                                        f"2\n+limbic\n(n={ns[2]})", f"3\n+neocort.\n(n={ns[3]})"])
    a.set_xlabel("Autopsy TDP-43 stage (tdp_st4)"); a.set_ylabel("% molecular TDP-43 positive")
    a.set_title("B   Molecular TDP-43+ vs autopsy stage", loc="left"); a.set_ylim(0, 105)
    a.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    ### concordance 
    ct = pd.crosstab(mp["molpos"], (mp.tdp_st4 >= 2)).reindex(index=[1.0, 0.0], columns=[False, True])
    _panel_concordance(ax[1, 0], ct, ["molecular\npositive", "molecular\nnegative"],
                       ["≤ amygdala\n(stage 0–1)", "beyond amygdala\n(stage 2–3)"],
                       "Autopsy TDP-43", "Molecular (≥1 cryptic exon)",
                       "C   Cross-modal concordance  (n={n}, agreement {agree:.0f}%)")

    ### panel
    a = ax[1, 1]
    if panelD == "dx":
        order = ["NCI", "MCI", "AD", "Other"]; colors = [pu.DX_COLORS[g] for g in order]
        xlabel = "Clinical diagnosis (cogdx)"; dtitle = "D   Molecular TDP-43+ by diagnosis"
    else:
        order = ["Braak 0–II", "Braak III–IV", "Braak V–VI"]
        colors = [pu.BRAAK3_COLORS[g] for g in order]
        xlabel = "Braak stage group (tau)"; dtitle = "D   Molecular TDP-43+ by Braak stage"
    gcol = "dx" if panelD == "dx" else "braak_grp3"
    vals = [df[df[gcol] == g]["molpos"].mean()*100 for g in order]
    ns = [int((df[gcol] == g).sum()) for g in order]
    x = np.arange(len(order))
    _labels(a, a.bar(x, vals, color=colors, edgecolor="white", width=0.62), "{:.0f}%")
    a.set_xticks(x); a.set_xticklabels([f"{g.replace('Braak ', '')}\n(n={n})" for g, n in zip(order, ns)])
    a.set_ylabel("% molecular TDP-43 positive"); a.set_xlabel(xlabel)
    a.set_title(dtitle, loc="left"); a.set_ylim(0, 105); a.grid(axis="y", alpha=.25)
    
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.965, color=INK)
    fig.text(0.075, 0.028, region_note + " Molecular TDP-43 positive = ≥1 STMN2/UNC13A cryptic exon; "
             "panels A–C relate it to autopsy tdp_st4.", fontsize=8.4, color="#555")
    pu.save(fig, stem, OUTDIR)


    
def autopsy_4panel_byBraak(df):
    pu.setup()
    mp = df.dropna(subset=["tdp_st4"]).copy()
    fig, ax = _fig_grid()

    counts = [int((mp.tdp_st4 == s).sum()) for s in [0, 1, 2, 3]]
    _panel_count_dist(ax[0, 0], counts, ["0\nnone", "1\namyg", "2\n+limbic", "3\n+neocort."], 1.5,
                      f"≤ amygdala\n(n={sum(counts[:2])})",
                      f"beyond amygdala\n(n={sum(counts[2:])})", C_AUT,
                      "Autopsy TDP-43 stage (tdp_st4)", "specimens (n)",
                      "A   Autopsy TDP-43 stage distribution")

    # beyond amygdala vs molecular cryptic-exon count
    a = ax[0, 1]; x = np.arange(4)
    by = [(mp[mp.n_expressed == k].tdp_st4 >= 2).mean()*100 for k in [0, 1, 2, 3]]
    ns = [int((mp.n_expressed == k).sum()) for k in [0, 1, 2, 3]]
    _labels(a, a.bar(x, by, color=C_AUT, edgecolor="white", width=0.6), "{:.0f}%")
    a.set_xticks(x); a.set_xticklabels([f"{k}\n(n={n})" for k, n in zip([0, 1, 2, 3], ns)])
    a.set_xlabel("Molecular cryptic-exon count (n_expressed)")
    a.set_ylabel("% TDP-43 beyond amygdala (autopsy)")

    
    a.set_title("B   Autopsy TDP-43 vs molecular cryptic exons", loc="left"); a.set_ylim(0, 60)
    a.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    ct = pd.crosstab((mp.tdp_st4 >= 2), mp["molpos"]).reindex(index=[True, False], columns=[1.0, 0.0])
    _panel_concordance(ax[1, 0], ct, ["beyond\namygdala", "≤ amygdala"],
                       ["molecular\npositive", "molecular\nnegative"],
                       "Molecular (≥1 cryptic exon)", "Autopsy TDP-43",
                       "C   Cross-modal concordance  (n={n}, agreement {agree:.0f}%)")

    ###   beyond amygdala by Braak
    a = ax[1, 1]; order = ["Braak 0–II", "Braak III–IV", "Braak V–VI"]; x = np.arange(3)
    vals = [(mp[mp.braak_grp3 == g].tdp_st4 >= 2).mean()*100 for g in order]
    means = [mp[mp.braak_grp3 == g].tdp_st4.mean() for g in order]
    ns = [int((mp.braak_grp3 == g).sum()) for g in order]
    bars = a.bar(x, vals, color=[pu.BRAAK3_COLORS[g] for g in order], edgecolor="white", width=0.62)
    for b, m in zip(bars, means):
        h = b.get_height()
        a.annotate(f"{h:.0f}%", (b.get_x()+b.get_width()/2, h), ha="center", va="bottom",
                   fontsize=10.5, fontweight="bold", xytext=(0, 2), textcoords="offset points", color=INK)
        a.annotate(f"mean {m:.2f}", (b.get_x()+b.get_width()/2, 1.2), ha="center", va="bottom",
                   fontsize=8.3, color="#555")
    a.set_xticks(x); a.set_xticklabels([f"{g.replace('Braak ', '')}\n(n={n})" for g, n in zip(order, ns)])
    a.set_ylabel("% TDP-43 beyond amygdala (autopsy)"); a.set_xlabel("Braak stage group (tau)")
    a.set_title("D   Autopsy TDP-43 by Braak stage", loc="left"); a.set_ylim(0, 60); a.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    fig.suptitle("Autopsy TDP-43 pathology — ROSMAP PCC RNA-seq samples", fontsize=15,
                 fontweight="bold", y=0.965, color=INK)
    fig.text(0.075, 0.028, "PCC samples with an autopsy TDP-43 stage. Panels A–C relate autopsy TDP-43 to the "
             "molecular call; panel D distributes it by Braak (tau) stage — a clear gradient (48% at Braak V–VI).",
             fontsize=8.4, color="#555")
    pu.save(fig, "ROSMAP_TDP43_figure_PCC_autopsy_byBraak", OUTDIR)


def main():
    pcc = dp.pcc_master()
    linked = dp._derive_common(dp.build_master())
    linked = linked[linked.tdp43_match == "matched"].copy()  # panel D uses linked-with-cogdx below
    all_linked = dp._derive_common(dp.build_master())

    molecular_4panel(all_linked[all_linked.projid.notna()], panelD="dx",
                     stem="ROSMAP_TDP43_figure_cryptic1plus",
                     title="Molecular TDP-43 (≥1 cryptic exon = positive) across modalities — ROSMAP RNA-seq specimens",
                     region_note="All linked specimens (n≈1,432; 602 with both measures).")
    
    molecular_4panel(pcc, panelD="dx")
    
    molecular_4panel(pcc, panelD="braak", stem="ROSMAP_TDP43_figure_PCC_byBraak",
                     title="Molecular TDP-43 (≥1 cryptic exon = positive) — ROSMAP PCC RNA-seq samples")
    autopsy_4panel_byBraak(pcc)


if __name__ == "__main__":
    main()
