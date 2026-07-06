"""
plot_utils.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.patches import Patch

INK = "#1a1a2e"
SC = ["#dfe3e6", "#F6C667", "#E8863C", "#B5341F"]
C_MOL = "#C44E52"     # molecular 
C_AUT = "#2E6E8E"     # autopsy
DX_COLORS = {"NCI": "#6Fae6F", "MCI": "#E0A458", "AD": "#C44E52", "Other": "#9a8fb0"}
BRAAK3_COLORS = {"Braak 0–II": "#8FbC8F", "Braak III–IV": "#E0A458", "Braak V–VI": "#C44E52"}


def setup():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.edgecolor": "#444", "axes.linewidth": 0.8,
        "axes.titlesize": 12.5, "axes.titleweight": "bold",
    })


def save(fig, stem, outdir):
    fig.savefig(f"{outdir}/{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{outdir}/{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _bar_labels(ax, bars, fmt="{:.0f}", dy=1):
    for b in bars:
        h = b.get_height()
        ax.annotate(fmt.format(h), (b.get_x() + b.get_width() / 2, h),
                    ha="center", va="bottom", fontsize=9.5,
                    xytext=(0, dy), textcoords="offset points", color=INK)


def two_panel_by_group(df, group_col, group_order, comp_col, comp_levels,
                       comp_labels, legend_title, rate_thresh, rate_ylabel,
                       suptitle, caption, left_title, right_title, outstem, outdir,
                       group_colors=None, right_ymax=60, xlabel="group"):
    
    setup()
    d = df.dropna(subset=[comp_col])
    tab = pd.crosstab(d[group_col], d[comp_col]).reindex(group_order).fillna(0)
    nper = tab.sum(1).astype(int)
    prop = tab.div(tab.sum(1), 0) * 100
    x = np.arange(len(group_order))
    if group_colors is None:
        group_colors = {g: list(DX_COLORS.values())[i % 4] for i, g in enumerate(group_order)}

    fig, ax = plt.subplots(1, 2, figsize=(13.6, 6.1))
    fig.subplots_adjust(wspace=0.30, top=0.80, bottom=0.16, left=0.07, right=0.985)

    # stacked composition
    a = ax[0]; bottom = np.zeros(len(group_order))
    for s, lev in enumerate(comp_levels):
        vals = prop[lev].values if lev in prop else np.zeros(len(group_order))
        a.bar(x, vals, bottom=bottom, color=SC[s], edgecolor="white", width=0.68)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 7:
                a.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                       fontsize=9.5, color="white" if s >= 3 else INK,
                       fontweight="bold" if s >= 2 else "normal")
        bottom += vals
    a.set_xticks(x); a.set_xticklabels([f"{g}\n(n={nper[g]})" for g in group_order])
    a.set_ylabel("% of participants"); a.set_ylim(0, 100); a.set_xlabel(xlabel)
    a.set_title(left_title, loc="left")
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    # rate + mean
    a = ax[1]
    rate = [(d[d[group_col] == g][comp_col] >= rate_thresh).mean() * 100 for g in group_order]
    meanv = [d[d[group_col] == g][comp_col].mean() for g in group_order]
    bars = a.bar(x, rate, color=[group_colors[g] for g in group_order],
                 edgecolor="white", width=0.62)
    for b, m in zip(bars, meanv):
        h = b.get_height()
        a.annotate(f"{h:.0f}%", (b.get_x() + b.get_width() / 2, h), ha="center",
                   va="bottom", fontsize=10.5, fontweight="bold",
                   xytext=(0, 2), textcoords="offset points", color=INK)
        a.annotate(f"mean {m:.2f}", (b.get_x() + b.get_width() / 2, 1.2),
                   ha="center", va="bottom", fontsize=8.3, color="#555")
    a.set_xticks(x); a.set_xticklabels([f"{g}\n(n={nper[g]})" for g in group_order])
    a.set_ylabel(rate_ylabel); a.set_ylim(0, max(right_ymax, max(rate) * 1.25))
    a.set_xlabel(xlabel); a.set_title(right_title, loc="left")
    a.grid(axis="y", alpha=.25)
    for s in ["top", "right"]:
        a.spines[s].set_visible(False)

    handles = [Patch(fc=SC[i], ec="white", label=comp_labels[i]) for i in range(len(comp_levels))]
    fig.legend(handles=handles, title=legend_title, frameon=False, ncol=len(comp_levels),
               loc="upper center", bbox_to_anchor=(0.5, 0.925), fontsize=9.5,
               title_fontsize=10, columnspacing=1.4, handlelength=1.3)
    fig.suptitle(suptitle, fontsize=14.5, fontweight="bold", y=0.985, color=INK)
    fig.text(0.07, 0.015, caption, fontsize=8.6, color="#555")
    save(fig, outstem, outdir)
    return dict(zip(group_order, rate))


def forest(ax, model, terms, title, sub, xmax, xticks):
    """Forest plot of exp(coef) with 95% CI from a fitted statsmodels Logit """
    ors = np.exp(model.params); ci = np.exp(model.conf_int()); ps = model.pvalues
    y = np.arange(len(terms))[::-1]
    for yi, (k, lab) in zip(y, terms):
        o = ors[k]; lo, hi = ci.loc[k]; p = ps[k]; sig = p < 0.05
        col = C_MOL if (sig and o > 1) else (C_AUT if (sig and o < 1) else "#9aa0a6")
        ax.plot([lo, min(hi, xmax * 0.999)], [yi, yi], color=col, lw=2.2,
                solid_capstyle="round", zorder=2)
        ax.plot(o, yi, "o", color=col, ms=8.5, mec="white", mew=1, zorder=3)
        star = "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"
        ax.text(xmax * 0.99, yi + 0.30, f"{o:.2f} [{lo:.2f}–{hi:.2f}] {star}",
                va="center", ha="right", fontsize=8, color="#333")
    ax.axvline(1, ls="--", lw=1, color="#888", zorder=1)
    ax.set_yticks(y); ax.set_yticklabels([l for _, l in terms], fontsize=9.3)
    ax.set_xscale("log"); ax.set_xlim(0.6, xmax); ax.set_xticks(xticks)
    ax.set_xticklabels([str(t) for t in xticks]); ax.set_xlabel("Odds ratio (log)")
    ax.set_title(title, loc="left", fontsize=11.5, fontweight="bold")
    ax.text(0, 1.10, sub, transform=ax.transAxes, fontsize=8.3, color="#666")
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(left=False); ax.set_ylim(-0.6, len(terms) - 0.4)
