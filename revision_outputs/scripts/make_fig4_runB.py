#!/usr/bin/env python3
"""Regenerate Figure 4 panels A, B and D from the Run B 80,000-cell obs.csv."""
import numpy as np, pandas as pd, matplotlib as mpl
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportion_confint

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "axes.linewidth": 0.8, "xtick.direction": "out", "ytick.direction": "out",
    "savefig.bbox": "tight", "savefig.dpi": 600,
})

REG = ["NT", "TP", "TC", "mLN"]
COMP = ["Epithelial-like", "Stromal-like", "Immune-like",
        "Endothelial-like", "Unclassified"]

df = pd.read_csv("obs.csv")
a, b = df["expr_WBP5"], df["expr_EGFR"]
df["W"] = a >= a[a > 0].quantile(.90)
df["E"] = b >= b[b > 0].quantile(.90)
df["D"] = df["W"] & df["E"]
df["Wo"] = df["W"] & ~df["D"]
df["Eo"] = df["E"] & ~df["D"]

# ---------------------------------------------------------------- panel A
groups = [("WBP5-high", "W"), ("EGFR-high", "E"), ("Double-high", "D")]
fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.6))
for ax, (lab, flag) in zip(axes, groups):
    sub = df[df[flag]]
    t = pd.crosstab(sub["region"], sub["compartment"]).reindex(
        index=REG, columns=COMP).fillna(0)
    pct = t.div(t.sum(axis=1), axis=0) * 100
    im = ax.imshow(pct.values, cmap="magma", vmin=0, vmax=75, aspect="auto")
    for i in range(len(REG)):
        for j in range(len(COMP)):
            v = pct.values[i, j]
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7.5,
                    color="white" if v < 45 else "black")
    ax.set_xticks(range(len(COMP)))
    ax.set_xticklabels(COMP, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(REG)))
    ax.set_yticklabels(REG, fontsize=9, fontweight="bold")
    ax.set_title(f"{lab} cells", fontsize=11, fontweight="bold", pad=14)
    ax.text(0.5, 1.02, f"n={len(sub):,}, row-normalized % within each region",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7.5,
            color="#555555")
    ax.set_xlabel("Proxy compartment", fontsize=9)
    if ax is axes[0]:
        ax.set_ylabel("Region", fontsize=9, fontweight="bold")
fig.colorbar(im, ax=axes, fraction=0.014, pad=0.015)
fig.savefig("Fig4A_runB.png")
fig.savefig("Fig4A_runB.pdf")
plt.close(fig)

# ---------------------------------------------------------------- panel B
colors = {"WBP5-high": "#2FBECD", "EGFR-high": "#E8503A", "Double-high": "#2C3E70"}
fig, ax = plt.subplots(figsize=(7.2, 4.4))
w = 0.26
for k, (lab, flag) in enumerate(groups):
    ys, los, his = [], [], []
    for r in REG:
        sub = df[df["region"] == r]
        n_pos, n = int(sub[flag].sum()), len(sub)
        lo, hi = proportion_confint(n_pos, n, method="wilson")
        ys.append(n_pos / n * 100); los.append(lo * 100); his.append(hi * 100)
    ys = np.array(ys)
    err = np.vstack([ys - np.array(los), np.array(his) - ys])
    ax.bar(np.arange(len(REG)) + (k - 1) * w, ys, w, label=lab,
           color=colors[lab], yerr=err, capsize=3,
           error_kw=dict(lw=0.9, capthick=0.9, ecolor="#333333"))
ax.set_xticks(range(len(REG)))
ax.set_xticklabels(REG, fontsize=10)
ax.set_xlabel("Region", fontsize=10)
ax.set_ylabel("Prevalence among all cells in region (%)", fontsize=10)
ax.legend(frameon=True, fontsize=9, edgecolor="#999999")
ax.spines[["top", "right"]].set_visible(False)
fig.savefig("Fig4B_runB.png"); fig.savefig("Fig4B_runB.pdf")
plt.close(fig)

# ---------------------------------------------------------------- panel D
grp = np.where(df["D"], "Double-high",
      np.where(df["Wo"], "WBP5-only",
      np.where(df["Eo"], "EGFR-only", None)))
s = df.assign(g=grp).dropna(subset=["g"])
order = ["WBP5-only", "EGFR-only", "Double-high"]
ns = {g: int((s["g"] == g).sum()) for g in order}
labels = [f"{g}\n(n={ns[g]:,})" for g in order]

reg_c = ["#4CAF50", "#2C5AA0", "#F0A030", "#B37FD0"]
comp_c = ["#E8503A", "#7EC8E3", "#3CB878", "#2C3E70", "#C8C8C8"]

fig, axes = plt.subplots(1, 2, figsize=(13.6, 4.6))
fig.subplots_adjust(wspace=0.62)
for ax, (var, levels, cols, title) in zip(axes, [
        ("region", REG, reg_c, "Region composition of three subsets"),
        ("compartment", COMP, comp_c, "Proxy compartment composition of three subsets")]):
    t = pd.crosstab(s["g"], s[var]).reindex(index=order, columns=levels).fillna(0)
    pct = t.div(t.sum(axis=1), axis=0) * 100
    bottom = np.zeros(len(order))
    for lv, c in zip(levels, cols):
        ax.bar(range(len(order)), pct[lv].values, 0.62, bottom=bottom,
               label=lv, color=c)
        bottom += pct[lv].values
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Composition within subset (%)", fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=16)
    ax.text(0.5, 1.015, "(80,000 quality-filtered cells)", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8, style="italic", color="#555555")
    ax.legend(bbox_to_anchor=(1.01, 0.5), loc="center left", frameon=False,
              fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
fig.savefig("Fig4D_runB.png", bbox_inches="tight"); fig.savefig("Fig4D_runB.pdf", bbox_inches="tight")
plt.close(fig)
print("[done] Fig4A/B/D regenerated")
