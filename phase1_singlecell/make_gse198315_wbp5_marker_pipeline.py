#!/usr/bin/env python3
"""
GSE198315 (multiregional OSCC) — WBP5 marker co-expression pipeline.

Computes WBP5 correlations with p-EMT / invasion markers (ITGA5, LAMC2, VIM,
KRT14, KRT17) at cell level and patient-mean level, stratified by a marker-based
proxy tumor-epithelial definition, and renders the summary and scatter figures.

Proxy tumor epithelial cells are defined as cells with non-zero expression of
any of EPCAM, KRT14, KRT17 or LAMC2 that originate from a tumour-associated
region (TP, TC or mLN).

Usage
-----
    python make_gse198315_wbp5_marker_pipeline.py \
        --data /path/to/GSE198315_OSCC_UMI_count_matrix.txt.gz \
        --outdir /path/to/output

Requires: numpy, pandas, scipy, matplotlib, seaborn, scienceplots
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
import scienceplots
from matplotlib import gridspec

plt.style.use(['science', 'nature', 'no-latex'])
sns.set_palette('colorblind')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Inter', 'Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'legend.fontsize': 8,
    'figure.dpi': 150,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

import argparse

_ap = argparse.ArgumentParser(
    description='GSE198315 WBP5 marker co-expression pipeline.')
_ap.add_argument('--data', required=True, type=Path,
                 help='GSE198315 UMI count matrix, gene x cell, tab-separated, '
                      'gzip-compressed (e.g. GSE198315_OSCC_UMI_count_matrix.txt.gz)')
_ap.add_argument('--outdir', required=True, type=Path,
                 help='Output directory (created if absent)')
_args = _ap.parse_args()

DATA = _args.data
OUTDIR = _args.outdir
OUTDIR.mkdir(parents=True, exist_ok=True)

MARKERS = ['ITGA5', 'LAMC2', 'VIM', 'KRT14', 'KRT17']
NEED = ['WBP5', 'EPCAM', 'KRT14', 'KRT17', 'LAMC2', 'ITGA5', 'VIM']
SITE_PALETTE = {'Primary': '#4DBBD5', 'Lymph node': '#E64B35'}
AXIS_CAPS = {
    'ITGA5': (0.0, 8.0),
    'LAMC2': (0.0, 10.0),
    'VIM': (0.0, 8.0),
    'KRT14': (0.0, 12.0),
    'KRT17': (0.0, 12.0),
}
MAX_PLOT_POINTS = 25000
RNG = np.random.default_rng(42)


def region_from_cell(cell_id: str) -> str:
    prefix = cell_id.split('_', 1)[0]
    if '-' in prefix:
        return prefix.split('-')[-1]
    return 'UNK'


def patient_from_cell(cell_id: str) -> str:
    prefix = cell_id.split('_', 1)[0]
    return prefix.split('-')[0]


def load_selected_matrix() -> pd.DataFrame:
    found: dict[str, np.ndarray] = {}
    with gzip.open(DATA, 'rt', encoding='utf-8', errors='replace') as f:
        header = f.readline().rstrip('\n').split('\t')
        cell_ids = header
        for line_idx, raw in enumerate(f, start=1):
            if '\t' not in raw:
                continue
            gene, rest = raw.rstrip('\n').split('\t', 1)
            if gene in NEED:
                found[gene] = np.fromstring(rest, sep='\t', dtype=np.float32)
                if len(found) == len(NEED):
                    break
            if line_idx % 5000 == 0:
                print(f'scanned {line_idx} rows; found {len(found)}/{len(NEED)} genes', flush=True)
    missing = [g for g in NEED if g not in found]
    if missing:
        raise SystemExit(f'Missing genes in GSE198315 matrix: {missing}')

    obs = pd.DataFrame({'cell_id': cell_ids})
    obs['patient'] = obs['cell_id'].map(patient_from_cell)
    obs['region'] = obs['cell_id'].map(region_from_cell)
    obs['site'] = np.where(obs['region'].eq('mLN'), 'Lymph node', 'Primary')
    for g in NEED:
        obs[g] = found[g]

    epithelial_like = (
        (obs['EPCAM'] > 0) |
        (obs['KRT14'] > 0) |
        (obs['KRT17'] > 0) |
        (obs['LAMC2'] > 0)
    )
    obs['is_epithelial_like'] = epithelial_like.astype(int)
    obs['is_tumor_region'] = obs['region'].isin(['TP', 'TC', 'mLN']).astype(int)
    obs['is_proxy_tumor_epithelial'] = (obs['is_epithelial_like'].eq(1) & obs['is_tumor_region'].eq(1)).astype(int)
    return obs


def safe_corr(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan, np.nan, np.nan
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_rho, spearman_p = spearmanr(x, y)
    return float(pearson_r), float(pearson_p), float(spearman_rho), float(spearman_p)


def append_row(rows: list[dict], df: pd.DataFrame, gene: str, subset: str, analysis: str):
    if analysis == 'cell_level':
        x = df['WBP5'].values
        y = df[gene].values
        pearson_r, pearson_p, spearman_rho, spearman_p = safe_corr(x, y)
        row = {
            'comparison': f'WBP5_vs_{gene}',
            'gene': gene,
            'analysis': analysis,
            'subset': subset,
            'n': int(len(df)),
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_rho': spearman_rho,
            'spearman_p': spearman_p,
            'WBP5_positive_frac': float((df['WBP5'] > 0).mean()) if len(df) else np.nan,
            f'{gene}_positive_frac': float((df[gene] > 0).mean()) if len(df) else np.nan,
            'double_positive_frac': float(((df['WBP5'] > 0) & (df[gene] > 0)).mean()) if len(df) else np.nan,
        }
    else:
        grouped = df.groupby('patient', observed=True)[['WBP5', gene]].mean().reset_index()
        x = grouped['WBP5'].values
        y = grouped[gene].values
        pearson_r, pearson_p, spearman_rho, spearman_p = safe_corr(x, y)
        row = {
            'comparison': f'WBP5_vs_{gene}',
            'gene': gene,
            'analysis': analysis,
            'subset': subset,
            'n': int(len(grouped)),
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_rho': spearman_rho,
            'spearman_p': spearman_p,
        }
    rows.append(row)


def build_correlation_tables(obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []
    proxy_tumor = obs[obs['is_proxy_tumor_epithelial'] == 1].copy()
    primary = proxy_tumor[proxy_tumor['region'].isin(['TP', 'TC'])].copy()
    mln = proxy_tumor[proxy_tumor['region'] == 'mLN'].copy()
    non_tumor = obs[obs['is_proxy_tumor_epithelial'] == 0].copy()

    subset_defs = [
        ('All cells', obs),
        ('Proxy tumor epithelial cells', proxy_tumor),
        ('Non-proxy tumor cells', non_tumor),
        ('Primary proxy tumor epithelial cells', primary),
        ('Lymph node proxy tumor epithelial cells', mln),
    ]
    patient_subset_defs = [
        ('Proxy tumor epithelial cells', proxy_tumor),
        ('Primary proxy tumor epithelial cells', primary),
        ('Lymph node proxy tumor epithelial cells', mln),
    ]

    for gene in ['LAMC2', 'ITGA5', 'VIM', 'KRT17', 'KRT14']:
        for subset_name, df in subset_defs:
            append_row(rows, df, gene, subset_name, 'cell_level')
        for subset_name, df in patient_subset_defs:
            append_row(rows, df, gene, subset_name, 'patient_mean')

    corr_df = pd.DataFrame(rows)
    summary = corr_df[(corr_df['analysis'] == 'cell_level') & (corr_df['subset'] == 'Proxy tumor epithelial cells')][
        ['gene', 'spearman_rho', 'spearman_p', 'pearson_r']
    ].copy()
    patient = corr_df[(corr_df['analysis'] == 'patient_mean') & (corr_df['subset'] == 'Proxy tumor epithelial cells')][
        ['gene', 'spearman_rho', 'spearman_p']
    ].copy().rename(columns={
        'spearman_rho': 'patient_mean_spearman_rho',
        'spearman_p': 'patient_mean_spearman_p',
    })
    non_tumor = corr_df[(corr_df['analysis'] == 'cell_level') & (corr_df['subset'] == 'Non-proxy tumor cells')][
        ['gene', 'spearman_rho']
    ].copy().rename(columns={'spearman_rho': 'non_proxy_spearman_rho'})
    summary = summary.rename(columns={
        'spearman_rho': 'proxy_tumor_spearman_rho',
        'spearman_p': 'proxy_tumor_spearman_p',
        'pearson_r': 'proxy_tumor_pearson_r',
    })
    summary = summary.merge(non_tumor, on='gene', how='left').merge(patient, on='gene', how='left')
    return corr_df, summary


def make_summary_figure(summary_df: pd.DataFrame, proxy_tumor_n: int, patient_n: int):
    order = ['ITGA5', 'LAMC2', 'KRT14', 'VIM', 'KRT17']
    cell = summary_df.set_index('gene').loc[order].reset_index()
    patient = cell.copy()
    patient['spearman_rho'] = patient['patient_mean_spearman_rho']
    patient['spearman_p'] = patient['patient_mean_spearman_p']
    patient['label'] = [f"ρ={v:.3f}" if pd.notna(v) else 'ρ=NA' for v in patient['spearman_rho']]
    patient['sig'] = patient['spearman_p'] < 0.05

    fig = plt.figure(figsize=(7.2, 4.8), constrained_layout=True)
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 1])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])

    colors = ['#4DBBD5' if g != 'KRT17' else '#B0B0B0' for g in order]
    barsA = axA.bar(order, cell['proxy_tumor_spearman_rho'], color=colors)
    axA.axhline(0, color='#777777', lw=0.8, ls='--')
    axA.set_ylim(-0.12, 0.30)
    axA.set_ylabel('Spearman rho')
    axA.set_title('Proxy tumor epithelial cell-level correlation', pad=18)
    axA.text(0.5, 1.01, f'GSE198315 proxy tumor epithelial cells, n={proxy_tumor_n:,}', transform=axA.transAxes,
             ha='center', va='bottom', fontsize=7, style='italic', color='#51616D')
    for bar, rho in zip(barsA, cell['proxy_tumor_spearman_rho']):
        y = bar.get_height()
        dy = 0.010 if y >= 0 else -0.018
        axA.text(bar.get_x() + bar.get_width() / 2, y + dy, f'ρ={rho:.3f}', ha='center', va='bottom', fontsize=7)

    barsB = axB.bar(order, patient['spearman_rho'], color=colors)
    axB.axhline(0, color='#777777', lw=0.8, ls='--')
    axB.set_ylim(-0.50, 1.00)
    axB.set_ylabel('Spearman rho')
    axB.set_title('Patient-mean correlation', pad=18)
    axB.text(0.5, 1.01, f'Mean proxy tumor epithelial expression across {patient_n} patients', transform=axB.transAxes,
             ha='center', va='bottom', fontsize=7, style='italic', color='#51616D')
    for bar, txt, sig in zip(barsB, patient['label'], patient['sig']):
        y = bar.get_height()
        dy = 0.03 if y >= 0 else -0.05
        axB.text(bar.get_x() + bar.get_width() / 2, y + dy, txt + ('*' if sig else ''), ha='center', va='bottom', fontsize=7)

    for ax, lbl in zip([axA, axB], 'AB'):
        ax.text(-0.14, 1.05, lbl, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')

    fig.suptitle('GSE198315: WBP5 correlation with EMT/epithelial markers', fontsize=12, fontweight='bold')
    fig.savefig(OUTDIR / 'Figure_GSE198315_WBP5_marker_summary.pdf')
    fig.savefig(OUTDIR / 'Figure_GSE198315_WBP5_marker_summary.png', dpi=600)
    plt.close(fig)


def make_scatter_panels(obs: pd.DataFrame) -> pd.DataFrame:
    proxy_tumor = obs[obs['is_proxy_tumor_epithelial'] == 1].copy()

    fig = plt.figure(figsize=(15.0, 11.4), constrained_layout=False)
    gs = gridspec.GridSpec(2, 3, figure=fig)
    axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(6)]

    stats_rows = []
    plot_total_n = 0
    filtered_total_n = 0

    for ax, gene, panel in zip(axes[:5], MARKERS, 'ABCDE'):
        pair_df = proxy_tumor.loc[(proxy_tumor['WBP5'] > 0) & (proxy_tumor[gene] > 0), ['WBP5', gene, 'site']].copy()
        filtered_total_n += len(pair_df)

        rho, p = spearmanr(pair_df['WBP5'], pair_df[gene]) if len(pair_df) >= 3 else (np.nan, np.nan)
        stats_rows.append({
            'gene': gene,
            'n_filtered_cells': int(len(pair_df)),
            'spearman_rho': float(rho) if pd.notna(rho) else np.nan,
            'spearman_p': float(p) if pd.notna(p) else np.nan,
            'primary_n': int((pair_df['site'] == 'Primary').sum()),
            'lymph_node_n': int((pair_df['site'] == 'Lymph node').sum()),
        })

        if len(pair_df) > MAX_PLOT_POINTS:
            plot_df = pair_df.sample(n=MAX_PLOT_POINTS, random_state=42).copy()
        else:
            plot_df = pair_df.sample(frac=1.0, random_state=42).copy()
        plot_total_n += len(plot_df)

        sns.scatterplot(
            data=plot_df, x='WBP5', y=gene, hue='site', palette=SITE_PALETTE,
            s=4, alpha=0.12, linewidth=0, ax=ax, legend=False
        )
        for site in ['Primary', 'Lymph node']:
            site_df = pair_df[pair_df['site'] == site]
            if len(site_df) >= 3 and site_df['WBP5'].nunique() > 1 and site_df[gene].nunique() > 1:
                sns.regplot(
                    data=site_df, x='WBP5', y=gene, scatter=False, ci=95, ax=ax,
                    color=SITE_PALETTE[site], line_kws={'lw': 2.2, 'alpha': 0.98}
                )

        x_q1, x_q99 = pair_df['WBP5'].quantile([0.05, 0.97])
        y_q1, y_q99 = pair_df[gene].quantile([0.05, 0.97])
        x_lower = max(0.0, float(x_q1) - 0.15)
        x_upper = min(8.0, float(x_q99) + 0.20)
        _, cap_y_upper = AXIS_CAPS[gene]
        y_lower = max(0.0, float(y_q1) - 0.20)
        y_upper = min(float(cap_y_upper), float(y_q99) + 0.30)
        if x_upper <= x_lower + 0.5:
            x_lower, x_upper = 0.0, 8.0
        if y_upper <= y_lower + 0.5:
            y_lower, y_upper = 0.0, float(cap_y_upper)
        ax.set_xlim(x_lower, x_upper)
        ax.set_ylim(y_lower, y_upper)

        ptxt = 'NS' if pd.isna(p) or p >= 0.05 else f'p={p:.1e}'
        rho_txt = 'NA' if pd.isna(rho) else f'{rho:.3f}'
        ax.set_title(f'{gene}\nSpearman ρ={rho_txt}, {ptxt}', pad=16, fontsize=10)
        ax.set_xlabel('WBP5 (UMI count)', fontsize=10)
        ax.set_ylabel('Expr (UMI count)', fontsize=10)
        ax.tick_params(axis='both', labelsize=8, pad=3)
        ax.text(-0.12, 1.04, panel, transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')

    axL = axes[5]
    axL.axis('off')
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=SITE_PALETTE['Primary'], markersize=6, alpha=0.8, label='Primary (TP+TC)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=SITE_PALETTE['Lymph node'], markersize=6, alpha=0.8, label='Lymph node (mLN)'),
        plt.Line2D([0], [0], color='#2F2F2F', lw=1.4, label='Linear fit (all filtered cells)'),
    ]
    axL.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.88),
               frameon=True, framealpha=0.95, facecolor='white', edgecolor='#E1E5EA',
               handletextpad=0.8, labelspacing=0.8, borderpad=0.8)
    axL.text(0.5, 0.28,
             'Proxy tumor epithelial cells only\n'
             'Filter: WBP5 > 0 and matched marker > 0\n'
             'x = WBP5, y = matched marker\n'
             f'full proxy tumor epithelial n={len(proxy_tumor):,}\n'
             f'filtered total={filtered_total_n:,}, plotted total={plot_total_n:,}\n'
             f'max plotted per panel={MAX_PLOT_POINTS:,}',
             ha='center', va='center', fontsize=9.2,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E1E5EA', alpha=0.95),
             transform=axL.transAxes)
    axL.text(0.02, 1.02, 'F', transform=axL.transAxes, fontsize=12, fontweight='bold', va='top')

    fig.suptitle('GSE198315: WBP5 vs markers (proxy tumor epithelial cells)', fontsize=12, fontweight='bold', y=0.97)
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.11, top=0.92, wspace=0.32, hspace=0.40)
    fig.savefig(OUTDIR / 'Figure_GSE198315_WBP5_marker_scatter_panels_filtered.pdf')
    fig.savefig(OUTDIR / 'Figure_GSE198315_WBP5_marker_scatter_panels_filtered.png', dpi=600)
    plt.close(fig)
    return pd.DataFrame(stats_rows)


def main():
    obs = load_selected_matrix()
    corr_df, summary_df = build_correlation_tables(obs)
    scatter_stats = make_scatter_panels(obs)

    proxy_tumor = obs[obs['is_proxy_tumor_epithelial'] == 1].copy()
    patient_n = proxy_tumor['patient'].nunique()
    make_summary_figure(summary_df, proxy_tumor_n=len(proxy_tumor), patient_n=patient_n)

    combined = summary_df.merge(scatter_stats, on='gene', how='left')

    corr_df.to_csv(OUTDIR / 'gse198315_wbp5_marker_correlations.tsv', sep='\t', index=False)
    summary_df.to_csv(OUTDIR / 'gse198315_wbp5_marker_correlations_summary.tsv', sep='\t', index=False)
    scatter_stats.to_csv(OUTDIR / 'gse198315_wbp5_marker_scatter_filtered_stats.tsv', sep='\t', index=False)
    combined.to_csv(OUTDIR / 'gse198315_wbp5_combined_summary.tsv', sep='\t', index=False)

    report = {
        'dataset': 'GSE198315 multiregional OSCC scRNA-seq',
        'n_total_cells': int(len(obs)),
        'n_proxy_tumor_epithelial_cells': int(len(proxy_tumor)),
        'n_patients_proxy_tumor_epithelial': int(patient_n),
        'regions_proxy_tumor_epithelial': proxy_tumor['region'].value_counts().to_dict(),
        'note': 'Because the deposited GEO matrix lacks author-provided malignant-cell annotation, the GSE103322 pipeline was adapted using marker-based proxy tumor epithelial cells (EPCAM/KRT14/KRT17/LAMC2-positive cells in TP/TC/mLN regions).',
    }
    (OUTDIR / 'gse198315_wbp5_marker_pipeline_report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
