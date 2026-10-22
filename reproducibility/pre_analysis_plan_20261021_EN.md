# Pre-analysis plan — four credibility analyses

**Recorded before the results were seen. 2026-10-21, KST.**

> **Note on this file.** This is an English translation of the Korean working
> document `pre_analysis_plan_20261021_KO.md`, which is deposited alongside it
> and is the authoritative record. The translation adds nothing to the original
> and omits nothing from it. The plan was written before any of the four
> analyses below was run, and the decision rules were not altered afterwards.
> Analysis outputs are in `results_cred/`.

Decision: the TCEAL/WBP gene-family specificity analysis is **included**. All
four analyses are to be reported in the response letter and in the
Supplementary Information **regardless of which way the results go**.

This document was written before execution. The criteria are not to be changed
after the results are seen. If a criterion must be changed, the reason for the
change is to be recorded separately.

---

## Common definitions

- **WBP5 correlation profile**: Spearman correlation on within-patient rank
  z-scores, taken as a partial correlation with respect to the technical axis.
  - **GSE198315**: tumour-region epithelium (Epithelial ∩ TC/TP); technical axis
    = log1p(total_counts). Same formulation and same gene set as step7b (D3)
    (detection rate ≥ 1%). Three subsamples (seeds 42 / 2026 / 7).
  - **GSE103322**: malignant cells (**2,200 cells, excluding the pooled "combo"
    plate**); technical axis = log1p(number of detected genes). Gene set as in
    P3 (variance > 0 in malignant cells).
- **Signature**: Puram p-EMT 100 genes (primary), top 30 genes (secondary).
- **Top-list size**: 200 (primary); 100 and 500 (sensitivity).
- **Reproduction check**: the new script must first recompute the existing D3
  (GSE198315) and P3 (GSE103322) profiles and agree with them. If they do not
  agree, none of the downstream results is to be used.

---

## 1. Cross-cohort overlap of gene identities

- **Primary result**: the list of names of Puram-100 genes that fall in the top
  200 of the adjusted GSE103322 profile **and** in the top 200 of all three
  GSE198315 seeds.
- **Significance**: within the set of genes common to both cohorts, a
  hypergeometric test of the overlap between each cohort's set of Puram genes in
  its top list. Secondarily, the whole-transcriptome top-200 overlap
  (hypergeometric) and the between-cohort Spearman correlation of the full
  profiles over the common genes.
- **Decision rule**: if the hypergeometric p for the Puram-gene overlap is
  < 0.05, state that "the same p-EMT genes are recovered in both cohorts".
  Otherwise report the names only and make no claim of reproducibility.

## 2. Leave-one-patient-out of the gene-level enrichment

- **Scope**: GSE198315, 10 runs for each of the three seeds; GSE103322
  (adjusted and unadjusted), as many runs as there are patients with ≥ 20 cells.
- **Statistics**: Puram-100 permutation p (10,000 permutations, expression-bin
  matched, as in D4); number of signature genes among the top 200.
- **Decision rule**: if p < 0.05 for every omission, state that "no single
  patient produces the result"; if p < 0.01 for every omission, state "strongly
  robust". If any single omission gives p ≥ 0.05, name that patient and report
  the result as it stands.

## 3. TCEAL/WBP gene-family specificity

- **Core family**: TCEAL1–8. **Extended** (Xq22 neighbours, specified in
  advance): BEX1–5 (BEX3 = NGFRAP1). **Outgroup**: WBP1, WBP2, WBP4, WBP11.
- **Evaluability threshold**: detection rate ≥ 10% in the analysed cells.
- **Statistics**: for each gene, the expression-bin-matched analytical z for the
  Puram 100 within that gene's own technical-axis-adjusted profile. The primary
  comparison is the percentile within **the 1,000 reference genes closest to
  WBP5 in expression level** (excluding Puram-100, family, outgroup and control
  genes). The p-value of the analytical z is overestimated because of
  between-gene correlation and is **not** used for the decision.
- **Decision rules**:
  - WBP5 at or above the 95th percentile against the reference genes in
    GSE103322 **and** in all three GSE198315 seeds → "p-EMT co-expression is
    pronounced relative to expression-matched genes".
  - WBP5's z higher than that of **every** evaluable core TCEAL member in both
    cohorts (all three seeds included) → "pronounced within the family as well".
  - At least half of the evaluable core members also at or above the 95th
    percentile in both cohorts → "a property shared across the family". To be
    stated as such in the response letter.
  - Otherwise → to be described as it stands, as an intermediate result.
- **3b, locus test**: adjust the WBP5 profile for the technical axis **plus** an
  Xq22 module score (the within-patient rank mean of the evaluable core and
  extended members, excluding WBP5), then recompute the Puram-100 permutation p.
  If p < 0.05 in GSE103322 and in at least two of the three seeds → "association
  beyond the Xq22 co-expression module". Otherwise → "not separable from the
  module", stating at the same time that the adjustment may be excessive because
  the whole module may respond to the same program.

## 4. Negative and positive controls

- **Negative controls** (specified in advance): the 11 uniformly expressed
  housekeeping genes of Eisenberg & Levanon, *Trends Genet* 2013 — C1orf43,
  CHMP2A, EMC7, GPI, PSMB2, PSMB4, RAB7A, REEP5, SNRPD3, VCP, VPS29.
- **Positive controls**: LAMC2, ITGA5 (the markers used in the manuscript). The
  percentile distribution of all Puram-100 members is also to be reported (each
  gene computed with itself excluded from the signature).
- **Decision rule**: if the median percentile of the 11 housekeeping genes lies
  between 25 and 75 and at most 1 of the 11 reaches the 95th percentile → "the
  analytical framework does not return positives for arbitrary genes". If the
  positive controls reach the 95th percentile → "the framework recovers known
  p-EMT genes".
- The proportion of the 1,000 reference genes with analytical z > 1.96 is also
  to be reported, as the baseline level of co-expression structure in the
  transcriptome.

---

## Outcome against these criteria

Recorded after execution, from `results_cred/` (see the response letter and
Supplementary Tables S14–S16 for the values):

| Criterion | Outcome |
|---|---|
| 1. Cross-cohort gene identity | **Not met** — only OCIAD2 and MAGED1 shared, p = 0.20. Names reported; no claim of reproducibility made. |
| 2. Leave-one-patient-out | **Met, strongly robust** — permutation p ≤ 1×10⁻⁴ in all runs. The one patient whose omission left no 30-gene overlap (HN25) is named and reported. |
| 3. Family specificity, ≥ 95th percentile in all four runs | **Not met** — WBP5 at the 89.7th–92.8th percentile. |
| 3. Family specificity, highest within the family | **Not met** — TCEAL3 exceeded WBP5 in all four runs. |
| 3. Family-shared property | **Not met** — 0/4, 1/4, 1/4, 1/4 of core members at ≥ 95th percentile. |
| 3b. Xq22 module adjustment | **Met** — association beyond the Xq22 module. |
| 4. Negative and positive controls | **Met** — housekeeping median percentile 44.9–56.8, none ≥ 95th; LAMC2 and ITGA5 ≥ 99.8th. |

The three unmet criteria under item 3 are reported in the manuscript and the
Supplementary Information as stated above, and the gene-family result is
presented as a limitation of the nomination rather than as support for it.
