# Addendum to the pre-analysis plan — published CAF signatures (CSBJ-D-26-00923R1)

Written 16 November 2026, after the pre-specified specificity analysis (pre_analysis_plan_R2_20261116.md) had been
run and read, and before the analysis described here was run. It is therefore a sensitivity analysis and is
labeled as such wherever it is reported.

## Why it was added

The reviewer asked for cancer-associated fibroblast (CAF) signatures. The pre-specified analysis used a
fibroblast-identity signature derived within GSE103322, which on inspection includes immediate-early genes not
specific to fibroblast identity. Published CAF signatures are added to test whether the fibroblast result depends
on that derived signature.

## Gene sets

The five pan-cancer CAF subtype signatures of Galbo, Zang and Zheng (Clin Cancer Res 2021), Supplementary Table S3,
used exactly as published: pan-myCAF (53 genes), pan-dCAF (69), pan-iCAF (42), pan-iCAF-2 (31), pan-pCAF (17).
No signature is published for the sixth subtype (pan-nCAF). Each signature is tested separately and all five are
reported. To be confirmed from that paper before any statement about independence is made: whether its head and
neck cohort is GSE103322.

## Analysis

Unchanged from the pre-specified analysis (script revision c, option --extra): same cells, profile, nulls,
mutually exclusive sets and two-way conditional analysis. Null draws exclude the p-EMT signature, the four
pre-specified programs and the five pan-CAF signatures. The pre-specified outputs are not recomputed.

## Reading rules (expression-matched null; discovery cohort, and at least two of three secondary subsamples)

- A pan-CAF signature is "an independent correlate of WBP5" if its unique genes remain enriched after adjustment
  for the p-EMT-unique module (p ≤ 0.05).
- The p-EMT signature is "not explained by" a pan-CAF signature if the p-EMT-unique genes remain enriched after
  adjustment for that signature's module (p ≤ 0.05).

| Outcome | Wording |
|---|---|
| At least one pan-CAF signature is an independent correlate in the discovery cohort | the fibroblast result of the discovery cohort is corroborated by a published CAF signature |
| None is | the discovery-cohort fibroblast result is specific to the derived signature and is reported as not corroborated |
| p-EMT is explained by a pan-CAF signature in both cohorts | the association is described as not separable from that CAF program |

No outcome of this analysis restores a claim of p-EMT specificity: the general-EMT result of the pre-specified
analysis stands whatever is found here.
