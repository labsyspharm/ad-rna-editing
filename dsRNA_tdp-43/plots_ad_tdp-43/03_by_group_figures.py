"""
03_by_group_figures.py
----------------------
Composition + rate figures for the PCC cohort, grouped by clinical diagnosis
and by Braak stage, for both the autopsy TDP-43 stage and the molecular
cryptic-exon measure for here plotting
"""
import os
import numpy as np
import pandas as pd


from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Patch



import data_prep as dp
import plot_utils as pu


OUTDIR = os.environ.get("ROSMAP_OUTDIR", "figures")
os.makedirs(OUTDIR, exist_ok=True)

TDP_LEVELS = [0, 1, 2, 3]
TDP_LABELS = ["0  none", "1  amygdala", "2  +limbic", "3  +neocort."]
CRY_LABELS = ["0  none", "1  exon", "2  exons", "3  exons"]



DX_ORDER = ["NCI", "MCI", "AD", "Other"]
BR3 = ["Braak 0–II", "Braak III–IV", "Braak V–VI"]


def by_diagnosis(pcc):
    pu.two_panel_by_group(
        pcc, "dx", DX_ORDER, "tdp_st4", TDP_LEVELS, TDP_LABELS,
        "TDP-43 stage (tdp_st4)", 2, "% TDP-43 beyond amygdala (stage 2–3)",
        "Autopsy TDP-43 pathology by clinical diagnosis - ROSMAP PCC RNA-seq samples",
        "PCC samples with an autopsy TDP-43 stage; participant-level.",
        "D   Autopsy TDP-43 stage composition", "D′   Advanced TDP-43 by diagnosis",
        "ROSMAP_TDP43_autopsy_by_diagnosis_PCC", OUTDIR,
        group_colors=pu.DX_COLORS, xlabel="Clinical diagnosis (cogdx)")

    pu.two_panel_by_group(
        pcc, "dx", DX_ORDER, "n_expressed", TDP_LEVELS, CRY_LABELS,
        "Cryptic exons (n_expressed)", 1, "% molecular TDP-43 positive (≥1 exon)",
        "Molecular TDP-43 (cryptic exons) by clinical diagnosis — ROSMAP PCC RNA-seq samples",
        "PCC samples with a clinical diagnosis; cryptic-exon count composition + % positive.",
        "M   Cryptic-exon count composition", "M′   Molecular TDP-43+ by diagnosis",
        "ROSMAP_TDP43_molecular_by_diagnosis_PCC", OUTDIR,
        group_colors=pu.DX_COLORS, right_ymax=100, xlabel="Clinical diagnosis (cogdx)")


def by_braak(pcc):
    pu.two_panel_by_group(
        pcc, "braak_grp3", BR3, "tdp_st4", TDP_LEVELS, TDP_LABELS,
        "TDP-43 stage (tdp_st4)", 2, "% TDP-43 beyond amygdala (stage 2–3)",
        "Autopsy TDP-43 pathology by Braak stage — ROSMAP PCC RNA-seq samples",
        "PCC samples with an autopsy TDP-43 stage. Braak: 0–II early, III–IV limbic, V–VI isocortical.",
        "D   Autopsy TDP-43 stage composition", "D′   Advanced TDP-43 by Braak",
        "ROSMAP_TDP43_autopsy_by_braak_PCC", OUTDIR,
        group_colors=pu.BRAAK3_COLORS, xlabel="Braak stage group (tau)")

    pu.two_panel_by_group(
        pcc, "braak_grp3", BR3, "n_expressed", TDP_LEVELS, CRY_LABELS,
        "Cryptic exons (n_expressed)", 1, "% molecular TDP-43 positive (≥1 exon)",
        "Molecular TDP-43 (cryptic exons) by Braak stage — ROSMAP PCC RNA-seq samples",
        "All PCC samples. Note the near-flat molecular signal across tau burden.",
        "M   Cryptic-exon count composition", "M′   Molecular TDP-43+ by Braak",
        "ROSMAP_TDP43_molecular_by_braak_PCC", OUTDIR,
        group_colors=pu.BRAAK3_COLORS, right_ymax=100, xlabel="Braak stage group (tau)")


def braak_binary_vs_autopsy(pcc):
    # do something
    pu.setup()
    GRP = ["Control (Braak 0–III)", "Case (Braak IV–VI)"]
    mp = pcc.dropna(subset=["tdp_st4"])
    ctrl = mp[mp.braak_bin == GRP[0]]["tdp_st4"]; case = mp[mp.braak_bin == GRP[1]]["tdp_st4"]
    _, pw = stats.mannwhitneyu(ctrl, case, alternative="two-sided")
    chi_b, pc_b, _, _ = stats.chi2_contingency(pd.crosstab(mp.braak_bin, (mp.tdp_st4 >= 2)))
    print(f"Braak-binary vs autopsy: mean {ctrl.mean():.2f} vs {case.mean():.2f} | "
          f"Wilcoxon p={pw:.2e} | Chi2(beyond) p={pc_b:.2e}")

    BC = {"Control (Braak 0–III)": "#8Fb0c4", "Case (Braak IV–VI)": "#C44E52"}
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 6.1))
    fig.subplots_adjust(wspace=0.34, top=0.82, bottom=0.17, left=0.08, right=0.975)

    a = ax[0]
    tab = pd.crosstab(mp.braak_bin, mp.tdp_st4).reindex(GRP)
    prop = tab.div(tab.sum(1), 0) * 100; nper = tab.sum(1).astype(int)
    x = np.arange(2); bottom = np.zeros(2)


    for s in range(4):
        vals = prop[float(s)].values if float(s) in prop else np.zeros(2)
        a.bar(x, vals, bottom=bottom, color=pu.SC[s], edgecolor="white", width=0.58)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 6:
                a.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center", fontsize=9.5,
                       color="white" if s >= 3 else pu.INK, fontweight="bold" if s >= 2 else "normal")
        bottom += vals
    a.set_xticks(x); a.set_xticklabels([f"{g}\n(n={nper[g]})" for g in GRP], fontsize=9.5)
    a.set_ylabel("% of participants"); a.set_ylim(0, 100)
    a.set_title("A   Autopsy TDP-43 stage composition", loc="left")
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)
    a.legend(handles=[Patch(fc=pu.SC[i], ec="white", label=l) for i, l in
                      enumerate(["TDP 0 none", "TDP 1 amygdala", "TDP 2 +limbic", "TDP 3 +neocort."])],
             frameon=False, fontsize=8.3, loc="lower center", bbox_to_anchor=(0.5, -0.32), ncol=2)

    a = ax[1]; sem = lambda s: s.std() / np.sqrt(len(s))
    
    for i, g in enumerate(GRP):
        s = mp[mp.braak_bin == g]["tdp_st4"]
        a.bar(i, s.mean(), 0.5, yerr=sem(s), capsize=4, color=BC[g], edgecolor="white",
              error_kw=dict(lw=1.2, ecolor="#555"), zorder=2)
        a.scatter(np.random.normal(i, 0.06, len(s)), s + np.random.normal(0, 0.05, len(s)),
                  s=6, color="#33333330", zorder=3)
        a.text(i, -0.30, f"n={len(s)}\n{(s >= 2).mean()*100:.0f}% beyond amyg",
               ha="center", fontsize=8, color="#666")
    star = lambda p: "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"
    a.plot([0, 0, 1, 1], [2.02, 2.10, 2.10, 2.02], lw=1.1, color="#444")
    a.text(0.5, 2.13, f"Wilcoxon p={pw:.1e} {star(pw)}  •  χ²(beyond) p={pc_b:.1e}",
           ha="center", va="bottom", fontsize=8.6)
    a.set_xticks([0, 1]); a.set_xticklabels(["Control\n(Braak 0–III)", "Case\n(Braak IV–VI)"])
    a.set_ylabel("Autopsy TDP-43 stage (tdp_st4)"); a.set_ylim(0, 2.7)
    a.set_title("B   TDP-43 stage: control vs case", loc="left"); a.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    fig.suptitle("Braak control vs case (0–III vs IV–VI) vs autopsy TDP-43 — ROSMAP PCC (n=546)",
                 fontsize=13.5, fontweight="bold", y=0.955, color=pu.INK)
    fig.text(0.08, 0.015, "High-tau cases carry significantly more TDP-43 pathology than controls "
             "— the opposite of the molecular cryptic-exon calls, which did not differ by Braak.",
             fontsize=8.4, color="#555")
    pu.save(fig, "ROSMAP_TDP43_braakBinary_vs_autopsy_PCC", OUTDIR)


def main():
    pcc = dp.pcc_master()
    by_diagnosis(pcc)
    by_braak(pcc)
    braak_binary_vs_autopsy(pcc)
    
    ### full cohort here
    cs = dp.radc_full()
    pu.two_panel_by_group(
        cs, "dx", DX_ORDER, "tdp_st4", TDP_LEVELS, TDP_LABELS,
        "TDP-43 stage (tdp_st4)", 2, "% TDP-43 beyond amygdala (stage 2–3)",
        "Autopsy TDP-43 pathology by clinical diagnosis — full RADC cohort",
        "All RADC participants with tdp_st4 + cogdx (participant-level).",
        "D   Autopsy TDP-43 stage composition", "D′   Advanced TDP-43 by diagnosis",
        "ROSMAP_TDP43_autopsy_by_diagnosis", OUTDIR,
        group_colors=pu.DX_COLORS, xlabel="Clinical diagnosis (cogdx)")


if __name__ == "__main__":
    main()
