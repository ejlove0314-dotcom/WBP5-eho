"""Shared style for CSBJ submission figures (vector PDF/EPS)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    # Liberation Sans is metric-compatible with Arial; Arial is used if installed.
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,   # embed TrueType -> text stays editable
    "svg.fonttype": "none",
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.titlesize": 8.5,
})
MM = 1 / 25.4
def letter(ax, s, x=-0.14, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="left")
def save(fig, stem, outdir):
    for ext in ("pdf", "eps"):
        fig.savefig(f"{outdir}/{stem}.{ext}", bbox_inches="tight", pad_inches=0.03, dpi=600)
    fig.savefig(f"{outdir}/{stem}_preview.png", dpi=300, bbox_inches="tight",
                pad_inches=0.03, facecolor="white")
