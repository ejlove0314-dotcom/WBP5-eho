"""Fig. 4 - one figure from the Run B 80,000-cell obs.csv / umap_coords.csv.
Same computations as make_fig4_runB.py (A, B, D) and make_fig4c_runB.py (C); panels are
redrawn together at final print size so text is uniform. UMAP point layers rasterised (600 dpi)."""
import sys, numpy as np, pandas as pd
from matplotlib.lines import Line2D
from statsmodels.stats.proportion import proportion_confint
from scipy.stats import chi2_contingency
sys.path.insert(0, sys.argv[1]); import _style as st; plt = st.plt
D, OUT = sys.argv[2], sys.argv[3]
df = pd.read_csv(f"{D}/obs.csv"); um = pd.read_csv(f"{D}/umap_coords.csv")
assert (um.barcode.values == df.barcode.values).all()
df["U1"], df["U2"] = um.UMAP1.values, um.UMAP2.values
REG = ["NT", "TP", "TC", "mLN"]; COMP = ["Epithelial-like", "Stromal-like", "Immune-like", "Endothelial-like", "Unclassified"]
a, b = df.expr_WBP5, df.expr_EGFR
tw, te = a[a > 0].quantile(.90), b[b > 0].quantile(.90)
df["W"] = a >= tw; df["E"] = b >= te; df["D"] = df.W & df.E; df["Wo"] = df.W & ~df.D; df["Eo"] = df.E & ~df.D
grp = np.where(df.D, "Double-high", np.where(df.Wo, "WBP5-only", np.where(df.Eo, "EGFR-only", None)))
s = df.assign(g=grp).dropna(subset=["g"]); order = ["WBP5-only", "EGFR-only", "Double-high"]
for var, lv in (("region", REG), ("compartment", COMP)):
    t = pd.crosstab(s.g, s[var]).reindex(index=order, columns=lv).fillna(0)
    t = t.loc[:, t.sum() > 0]; chi, p, dof, _ = chi2_contingency(t.values)
    print(f"chi2 {var}: {chi:.1f}, df={dof}, p={p:.2e}")
print(f"thresholds W {tw:.4f} E {te:.4f}; n W={int(df.W.sum())} E={int(df.E.sum())} D={int(df.D.sum())} "
      f"Wo={int(df.Wo.sum())} Eo={int(df.Eo.sum())} other={int((~df.W & ~df.E).sum())}")
print("region n", df.region.value_counts().reindex(REG).to_dict())

fig = plt.figure(figsize=(174*st.MM, 225*st.MM))
gs = fig.add_gridspec(3, 1, height_ratios=[0.85, 1.55, 0.95], hspace=0.42, left=0.08, right=0.97, top=0.965, bottom=0.04)
# ---- A
gA = gs[0].subgridspec(1, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.12)
groups = [("WBP5-high", "W"), ("EGFR-high", "E"), ("Double-high", "D")]
for k, (lab, flag) in enumerate(groups):
    ax = fig.add_subplot(gA[k]); sub = df[df[flag]]
    t = pd.crosstab(sub.region, sub.compartment).reindex(index=REG, columns=COMP).fillna(0)
    pct = t.div(t.sum(axis=1), axis=0)*100
    im = ax.imshow(pct.values, cmap="magma", vmin=0, vmax=75, aspect="auto")
    for i in range(4):
        for j in range(5):
            v = pct.values[i, j]; ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5.5,
                                         color="white" if v < 45 else "black")
    ax.set_xticks(range(5)); ax.set_xticklabels(COMP, rotation=45, ha="right", rotation_mode="anchor", fontsize=6)
    ax.set_yticks(range(4)); ax.set_yticklabels(REG if k == 0 else [], fontsize=6.5)
    ax.set_title(f"{lab} cells (n = {len(sub):,})", fontweight="bold", fontsize=7.5)
    ax.tick_params(length=2)
    if k == 0: ax.set_ylabel("Region"); st.letter(ax, "(A)", x=-0.32, y=1.08)
cb = fig.colorbar(im, cax=fig.add_subplot(gA[3])); cb.set_label("Row-normalized %", fontsize=6.5)
cb.ax.tick_params(labelsize=6, length=2); cb.outline.set_linewidth(0.5)
# ---- B and D
g2 = gs[1].subgridspec(1, 2, width_ratios=[0.8, 1.5], wspace=0.25)
gD = gs[2].subgridspec(1, 2, wspace=0.75)
ax = fig.add_subplot(g2[0]); col = {"WBP5-high": "#2FBECD", "EGFR-high": "#E8503A", "Double-high": "#2C3E70"}
for k, (lab, flag) in enumerate(groups):
    ys, lo_, hi_ = [], [], []
    for r in REG:
        sub = df[df.region == r]; npos, n = int(sub[flag].sum()), len(sub)
        lo, hi = proportion_confint(npos, n, method="wilson"); ys.append(npos/n*100); lo_.append(lo*100); hi_.append(hi*100)
    ys = np.array(ys)
    ax.bar(np.arange(4) + (k - 1)*0.26, ys, 0.26, label=lab, color=col[lab],
           yerr=np.vstack([ys - lo_, np.array(hi_) - ys]), capsize=1.5, error_kw=dict(lw=0.6, capthick=0.6, ecolor="#333"))
ax.set_xticks(range(4)); ax.set_xticklabels(REG); ax.set_xlabel("Region"); ax.set_ylabel("Prevalence among cells in region (%)", fontsize=7)
ax.set_ylim(0, 5.6); ax.legend(frameon=False, fontsize=6, loc="upper left", handlelength=1, ncol=1)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
st.letter(ax, "(B)", x=-0.3, y=1.05)
labels = [f"{g}\n(n={int((s.g == g).sum()):,})" for g in order]
for kk, (var, lv, cols, title) in enumerate([("region", REG, ["#4CAF50", "#2C5AA0", "#F0A030", "#B37FD0"], "Region"),
        ("compartment", COMP, ["#E8503A", "#7EC8E3", "#3CB878", "#2C3E70", "#C8C8C8"], "Proxy compartment")]):
    ax = fig.add_subplot(gD[kk]); t = pd.crosstab(s.g, s[var]).reindex(index=order, columns=lv).fillna(0)
    pct = t.div(t.sum(axis=1), axis=0)*100; bot = np.zeros(3)
    for l, c in zip(lv, cols): ax.bar(range(3), pct[l].values, 0.62, bottom=bot, label=l, color=c); bot += pct[l].values
    ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=5.5); ax.set_ylim(0, 100)
    ax.set_ylabel("Composition within subset (%)", fontsize=7); ax.set_title(f"{title} composition", fontweight="bold", fontsize=7.5)
    ax.legend(bbox_to_anchor=(1.0, 0.5), loc="center left", frameon=False, fontsize=5.5, handlelength=0.9)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    if kk == 0: st.letter(ax, "(D)", x=-0.3, y=1.05)
# ---- C
GREY, TEAL, RED, NAVY = "#D9D9D9", "#2AA198", "#E8503A", "#1F3A93"
gC = g2[1].subgridspec(2, 3, width_ratios=[1, 1, 0.62], wspace=0.08, hspace=0.32)
panels = [("WBP5-high", [("W", TEAL, 0.6)], f"≥ {tw:.2f}; n = {int(df.W.sum()):,}"),
          ("EGFR-high", [("E", RED, 0.6)], f"≥ {te:.2f}; n = {int(df.E.sum()):,}"),
          ("Double-high", [("D", NAVY, 1.6)], f"n = {int(df.D.sum()):,}"),
          ("Combined", [("Eo", RED, 0.6), ("Wo", TEAL, 0.6), ("D", NAVY, 1.6)], "mutually exclusive")]
for k, (t, layers, sub) in enumerate(panels):
    ax = fig.add_subplot(gC[k // 2, k % 2])
    ax.scatter(df.U1, df.U2, s=0.15, c=GREY, lw=0, rasterized=True)
    for f, c, sz in layers:
        q = df[df[f]]; ax.scatter(q.U1, q.U2, s=sz, c=c, lw=0, rasterized=True)
    ax.set_title(t, fontweight="bold", fontsize=7.5, pad=9); ax.text(0.5, 1.01, sub, transform=ax.transAxes, ha="center",
                                                              va="bottom", fontsize=5.5, color="#555")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1", fontsize=6.5)
    if k % 2 == 0: ax.set_ylabel("UMAP2", fontsize=6.5)
    if k == 0: st.letter(ax, "(C)", x=-0.14, y=1.13)
    for sp in ax.spines.values(): sp.set_linewidth(0.5)
lax = fig.add_subplot(gC[:, 2]); lax.axis("off")
lax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, mfc=c, mec="none", label=l) for c, l in
                    [(GREY, f"Other (n={int((~df.W & ~df.E).sum()):,})"), (RED, f"EGFR-high only (n={int(df.Eo.sum()):,})"),
                     (TEAL, f"WBP5-high only (n={int(df.Wo.sum()):,})"), (NAVY, f"Double-high (n={int(df.D.sum()):,})")]],
           loc="center left", frameon=False, fontsize=5.5, handletextpad=0.3, borderaxespad=0)
st.save(fig, "Figure_4", OUT)
