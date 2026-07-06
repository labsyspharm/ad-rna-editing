"""
05_stats_forest_mediation.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.mediation import Mediation
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import data_prep as dp
import plot_utils as pu

OUTDIR = os.environ.get("ROSMAP_OUTDIR", "figures")
os.makedirs(OUTDIR, exist_ok=True)
INK = pu.INK

## do only for PCC here
def association_battery(pcc):
    d = pcc.copy()
    # do something
    d["beyond"] = (d.tdp_st4 >= 2).astype("float"); d.loc[d.tdp_st4.isna(), "beyond"] = np.nan
    d["classhigh"] = (d.n_expressed >= 2).astype(int)
    R = []

    def sp(x, y, label, blk):
        s = d[[x, y]].dropna(); rho, p = stats.spearmanr(s[x], s[y])
        R.append((blk, label, len(s), "Spearman", f"ρ={rho:+.3f}", p))

    def mw(gcol, vcol, g1, g2, label, blk):
        a = d[d[gcol] == g1][vcol].dropna(); b = d[d[gcol] == g2][vcol].dropna()
        U, p = stats.mannwhitneyu(a, b, alternative="two-sided"); rb = 1 - 2*U/(len(a)*len(b))
        R.append((blk, label, len(a)+len(b), "Mann-Whitney", f"rank-biserial={-rb:+.3f}", p))

    def tab2(bx, by, label, blk):
        s = d[[bx, by]].dropna(); ct = pd.crosstab(s[bx], s[by])
        if ct.shape == (2, 2):
            orr, p = stats.fisher_exact(ct.values); R.append((blk, label, len(s), "Fisher 2×2", f"OR={orr:.2f}", p))
        else:
            chi, p, _, _ = stats.chi2_contingency(ct); R.append((blk, label, len(s), "chi²", f"χ²={chi:.1f}", p))

    d["tdp_group"] = d["tdp_st4"].map(lambda v: "Stage 0" if v == 0 else ("Stage 3" if v == 3 else np.nan))
    # for pathology
    sp("tdp_st4", "braaksc", "TDP-43 stage vs Braak", "1 Pathology")
    sp("tdp_st4", "ceradsc", "TDP-43 stage vs CERAD", "1 Pathology")
    sp("tdp_st4", "cogdx", "TDP-43 stage vs diagnosis", "1 Pathology")
    sp("tdp_st4", "age_death", "TDP-43 stage vs age at death", "1 Pathology")
    tab2("beyond", "apoe4b", "TDP-43 beyond-amygdala vs APOE e4", "1 Pathology")
    tab2("beyond", "msex", "TDP-43 beyond-amygdala vs sex(male)", "1 Pathology")
    
    # molecular vs biology here now
    sp("n_expressed", "tdp_st4", "Cryptic count vs TDP-43 stage", "2 Molecular↔biology")
    mw("tdp_group", "n_expressed", "Stage 0", "Stage 3", "Cryptic: Stage 0 vs Stage 3", "2 Molecular↔biology")
    tab2("classhigh", "beyond", "Cryptic ≥2 vs TDP-43 beyond-amygdala", "2 Molecular↔biology")
    sp("n_expressed", "cogdx", "Cryptic count vs diagnosis", "2 Molecular↔biology")
    tab2("STMN2short_b", "beyond", "STMN2 cryptic vs TDP-43 beyond-amygdala", "2 Molecular↔biology")
    
    # technical driver
    sp("n_expressed", "RIN", "Cryptic count vs RIN", "3 Technical driver")
    sp("STMN2short_b", "RIN", "STMN2 cryptic vs RIN", "3 Technical driver")
    sp("n_expressed", "pmi", "Cryptic count vs PMI", "3 Technical driver")
    sp("tdp_st4", "RIN", "TDP-43 stage vs RIN [control]", "3 Technical driver")
    # marker co-occurrence for now
    tab2("STMN2short_b", "UNC13A_CE1_b", "STMN2 vs UNC13A-CE1", "4 Marker co-occur")
    tab2("UNC13A_CE1_b", "UNC13A_CE2_b", "UNC13A-CE1 vs UNC13A-CE2", "4 Marker co-occur")

    res = pd.DataFrame(R, columns=["block", "comparison", "n", "test", "effect_size", "p"])
    res["q_FDR"] = multipletests(res["p"], method="fdr_bh")[1]
    res = res.sort_values("p")
    res.to_csv(f"{OUTDIR}/ROSMAP_TDP43_PCC_stats_results.csv", index=False)
    print(res.to_string(index=False))
    return res


## forest ploit for PCC - in 2 panels here
def _fit(df, y, xs):
    s = df[[y]+xs].apply(pd.to_numeric, errors="coerce").dropna().astype(float)
    m = sm.Logit(s[y], sm.add_constant(s[xs])).fit(disp=0)
    return m, len(s), int(s[y].sum())


def pcc_forest(pcc):
    pu.setup()
    d = pcc.copy()
    d["beyond"] = (d.tdp_st4 >= 2).astype(float); d.loc[d.tdp_st4.isna(), "beyond"] = np.nan
    d["male"] = d["msex"]; d["cerad_plaque"] = 4 - d["ceradsc"]
    d["AD"] = (d.cogdx.isin([4, 5])).astype(float); d.loc[d.cogdx.isna(), "AD"] = np.nan

    mA, nA, eA = _fit(d, "beyond", ["braaksc", "cerad_plaque", "apoe4b", "male"])
    A = [("braaksc", "Braak stage (tangles)"), ("cerad_plaque", "CERAD (plaque load)"),
         ("apoe4b", "APOE ε4 carrier"), ("male", "Male sex")]
    dz = (d[["AD", "tdp_st4", "braaksc", "cerad_plaque"]]
          .apply(pd.to_numeric, errors="coerce").dropna().astype(float))
    for c in ["tdp_st4", "braaksc", "cerad_plaque"]:
        dz[c] = (dz[c] - dz[c].mean()) / dz[c].std()
    mB = sm.Logit(dz["AD"], sm.add_constant(dz[["tdp_st4", "braaksc", "cerad_plaque"]])).fit(disp=0)
    B = [("tdp_st4", "TDP-43 stage"), ("braaksc", "Braak (tangles)"), ("cerad_plaque", "CERAD (plaques)")]

    fig, axs = plt.subplots(1, 2, figsize=(13.6, 5.4))
    fig.subplots_adjust(wspace=0.55, top=0.80, bottom=0.15, left=0.20, right=0.965)
    pu.forest(axs[0], mA, A, "A  Predictors of advanced TDP-43",
              f"beyond amygdala · adjusted · n={nA}, events={eA}", 5, [1, 2, 4])
    pu.forest(axs[1], mB, B, "B  Independent contributors to AD dementia",
              f"odds per 1 SD · n={len(dz)}", 3, [1, 2, 3])
    fig.suptitle("PCC — multivariable models of TDP-43 pathology and dementia",
                 fontsize=14, fontweight="bold", y=0.955, color=INK)
    pu.save(fig, "ROSMAP_TDP43_PCC_forest", OUTDIR)



# full cohort forests with mediation here for now....

def _mediate(cs, exp, med, cov, n_rep=800):
    cols = ["AD", exp, med] + cov
    dd = cs[cols].dropna().copy()
    om = smf.logit(f"AD ~ {exp} + {med}" + (" + " + " + ".join(cov) if cov else ""), data=dd)
    mm = smf.ols(f"{med} ~ {exp}" + (" + " + " + ".join(cov) if cov else ""), data=dd)
    s = Mediation(om, mm, exp, med).fit(n_rep=n_rep).summary()
    pm = s.loc["Prop. mediated (average)", "Estimate"] * 100
    lo = s.loc["Prop. mediated (average)", "Lower CI bound"] * 100
    hi = s.loc["Prop. mediated (average)", "Upper CI bound"] * 100
    return pm, lo, hi


def fullcohort_summary(cs):
    pu.setup()
    mA, nA, _ = _fit(cs, "beyond", ["age10", "braaksc", "cerad_plaque", "apoe4b", "male"])
    A = [("age10", "Age at death (per decade)"), ("apoe4b", "APOE ε4 carrier"),
         ("cerad_plaque", "CERAD (plaques)"), ("braaksc", "Braak (tangles)"), ("male", "Male sex")]
    mB, nB, _ = _fit(cs, "AD", ["tdp_st4", "HS", "braaksc", "cerad_plaque", "apoe4b", "male", "age10"])
    B = [("HS", "Hippocampal sclerosis"), ("age10", "Age (per decade)"),
         ("braaksc", "Braak (tangles)"), ("tdp_st4", "TDP-43 stage"),
         ("cerad_plaque", "CERAD (plaques)"), ("apoe4b", "APOE ε4")]

    print("Running mediation bootstraps (this takes ~1 min)...")
    
    MED = {"Age at death": {"TDP-43": _mediate(cs, "age10", "tdp_st4", ["male"]),
                            "Tau (Braak)": _mediate(cs, "age10", "braaksc", ["male"])},
           "APOE ε4": {"TDP-43": _mediate(cs, "apoe4b", "tdp_st4", ["age10", "male"]),
                       "Tau (Braak)": _mediate(cs, "apoe4b", "braaksc", ["age10", "male"])}}

    fig = plt.figure(figsize=(15.5, 5.8))
    
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.05, 0.95], wspace=0.62)
    fig.subplots_adjust(top=0.79, bottom=0.16, left=0.135, right=0.975)
    pu.forest(fig.add_subplot(gs[0]), mA, A, "A  Predictors of advanced TDP-43",
              f"beyond amygdala · n={nA}", 3.2, [1, 2, 3])
    pu.forest(fig.add_subplot(gs[1]), mB, B, "B  Contributors to AD dementia",
              f"adjusted incl. hippocampal sclerosis · n={nB}", 6, [1, 2, 3, 4, 6])
    ax2 = fig.add_subplot(gs[2]); exps = list(MED); w = 0.36; x = np.arange(len(exps))
    COL = {"TDP-43": "#C44E52", "Tau (Braak)": "#4C72B0"}

    # do something here for above
    for i, med in enumerate(["TDP-43", "Tau (Braak)"]):
        vals = [MED[e][med][0] for e in exps]
        los = [MED[e][med][0]-MED[e][med][1] for e in exps]; his = [MED[e][med][2]-MED[e][med][0] for e in exps]
        ax2.bar(x+(i-0.5)*w, vals, w, yerr=[los, his], capsize=4, color=COL[med], label=med,
                edgecolor="white", error_kw=dict(lw=1.2, ecolor="#555"))
        for xi, v in zip(x+(i-0.5)*w, vals):
            ax2.text(xi, v+2.5, f"{v:.0f}%", ha="center", fontsize=9, fontweight="bold")
    ax2.set_xticks(x); ax2.set_xticklabels(exps); ax2.set_ylim(0, 60)
    ax2.set_ylabel("% of effect on dementia mediated")
    ax2.set_title("C  Mediation of dementia risk", loc="left", fontsize=11.5, fontweight="bold")
    ax2.legend(frameon=False, fontsize=8.8, loc="upper right"); ax2.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        ax2.spines[s].set_visible(False)
    fig.suptitle("Full RADC cohort (n≈2,000) — TDP-43 pathology, hippocampal sclerosis, and dementia risk",
                 fontsize=14.5, fontweight="bold", y=0.95, color=INK)
    pu.save(fig, "ROSMAP_TDP43_fullcohort_summary", OUTDIR)



## Braak-stratified

def paper_approach(pcc):
    pu.setup()
    d = pcc.copy()
    d["tdp_pos"] = d["tdp_st4"].apply(lambda v: np.nan if pd.isna(v)
                                      else ("TDP-43 positive" if v >= 1 else "TDP-43 negative"))
    BORD = ["0–III", "IV", "V", "VI"]; sem = lambda s: s.std()/np.sqrt(len(s))
    fig, ax = plt.subplots(1, 2, figsize=(13, 6.0))
    fig.subplots_adjust(wspace=0.32, top=0.82, bottom=0.14, left=0.08, right=0.975)

    GC = {"0–III": "#9DB7C9", "IV": "#E0A458", "V": "#D9743F", "VI": "#B5341F"}
    a = ax[0]; ctrl = d[d.braak_grp4 == "0–III"]["n_expressed"].dropna()
    
    for i, g in enumerate(BORD):
        s = d[d.braak_grp4 == g]["n_expressed"].dropna()
        a.bar(i, s.mean(), 0.62, yerr=sem(s), capsize=4, color=GC[g], edgecolor="white",
              error_kw=dict(lw=1.2, ecolor="#555"), zorder=2)
        a.scatter(np.random.normal(i, 0.07, len(s)), s+np.random.normal(0, 0.04, len(s)),
                  s=6, color="#33333333", zorder=3)
        a.text(i, -0.18, f"n={len(s)}", ha="center", fontsize=8, color="#666")
    _, pV = stats.mannwhitneyu(ctrl, d[d.braak_grp4 == "V"]["n_expressed"].dropna())
    a.plot([0, 0, 2, 2], [3.1, 3.15, 3.15, 3.1], lw=1.1, color="#444")
    a.text(1, 3.17, f"0–III vs V: p={pV:.2f} (ns)", ha="center", fontsize=9)
    a.set_xticks(range(4)); a.set_xticklabels(BORD); a.set_ylim(0, 3.6)
    a.set_ylabel("Cryptic exons detected (n_expressed)")
    a.set_xlabel("Braak stage stratum  (0–III = control)")
    a.set_title("A   Cryptic-exon burden by Braak stage", loc="left"); a.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    a = ax[1]; groups = ["TDP-43 negative", "TDP-43 positive"]; cols = ["#8Fb0c4", "#C44E52"]

    ## additional info
    for i, g in enumerate(groups):
        s = d[d.tdp_pos == g]["n_expressed"].dropna()
        a.bar(i, s.mean(), 0.5, yerr=sem(s), capsize=4, color=cols[i], edgecolor="white",
              error_kw=dict(lw=1.2, ecolor="#555"), zorder=2)
        a.scatter(np.random.normal(i, 0.06, len(s)), s+np.random.normal(0, 0.04, len(s)),
                  s=6, color="#33333330", zorder=3)
        a.text(i, -0.18, f"n={len(s)}\n{(s >= 1).mean()*100:.0f}% ≥1 exon", ha="center", fontsize=8, color="#666")
    pos = d[d.tdp_pos == "TDP-43 positive"]["n_expressed"].dropna()
    neg = d[d.tdp_pos == "TDP-43 negative"]["n_expressed"].dropna()
    _, pw = stats.mannwhitneyu(pos, neg)
    dd = d.dropna(subset=["tdp_pos"]); ct = pd.crosstab(dd.tdp_pos, (dd.n_expressed >= 1))
    chi, pc, _, _ = stats.chi2_contingency(ct)
    a.plot([0, 0, 1, 1], [2.02, 2.07, 2.07, 2.02], lw=1.1, color="#444")
    a.text(0.5, 2.09, f"Wilcoxon p={pw:.3f}  •  χ² p={pc:.3f}", ha="center", fontsize=9)
    a.set_xticks([0, 1]); a.set_xticklabels(["TDP-43\nnegative", "TDP-43\npositive"]); a.set_ylim(0, 2.6)
    a.set_ylabel("Cryptic exons detected (n_expressed)")
    
    a.set_xlabel("Autopsy TDP-43 status (stage 0 vs ≥1)")
    a.set_title("B   Cryptic-exon burden by TDP-43 status", loc="left"); a.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    fig.suptitle("Paper's approach applied to ROSMAP PCC — cryptic-exon burden by Braak & TDP-43 status",
                 fontsize=14, fontweight="bold", y=0.95, color=INK)
    fig.text(0.08, 0.015, "n_expressed (0–3 cryptic-exon count) as a proxy for the paper's Log2(CPM) "
             "— raw reads/CPM are not in the provided files (see pipeline/). Wilcoxon rank-sum & Chi-square.",
             fontsize=8.3, color="#555")
    pu.save(fig, "ROSMAP_TDP43_paper_approach_PCC", OUTDIR)



def main():
    pcc = dp.pcc_master()
    association_battery(pcc)
    pcc_forest(pcc)
    fullcohort_summary(dp.radc_full())
    paper_approach(pcc)


if __name__ == "__main__":
    main()
