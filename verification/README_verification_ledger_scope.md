# Scope of `verification_ledger.csv` and `VERIFICATION.md`

Both files in this directory were compiled between July and August 2026, while
the manuscript was under review. `VERIFICATION.md` is the narrative audit report
and records the manuscript it audited (`CSBJ_WBP5_v10_20260825`, audited
2026-08-27); `verification_ledger.csv` is its machine-readable companion. Both
use the section, table and figure numbering of the **submitted** version. The
values they record are current; the locations are not. One caveat applies.

## Locations that no longer exist

The revision moved the structural and generative-chemistry material into the
Supplementary Information. The revised manuscript has Tables 1 and 2 and
Figs. 1 to 5, and its Results end at Section 3.3.

Thirteen ledger rows refer to `Section 3.4`, `Section 3.5`, `Section 3.6`,
`Section 3.7`, `Table 3`, `Table 4` or `Fig. 6A`, and `VERIFICATION.md` uses the
same labels in its narrative. Those labels belong to the submitted version. The quantities themselves are unchanged and are now reported
in the Supplementary Information, chiefly in Supplementary Note S1 and Tables
S13 to S18. For example, the four rows giving the P0 pocket descriptors
(Drug Score 0.837, volume 229.6, surface 811.5, depth 24.9; the twelve
pocket-lining residues; mean pLDDT 79.4; eight detected pockets) are now
Table S14, and the rows for the docking replicates and the compound profiles
are now Tables S13 and S15.

## Rows revised with the manuscript

Two Section 3.3 rows recorded values that the revision replaced. Both have been
updated in the file itself, using the `supersedes` notation the ledger already
uses elsewhere, so that the superseded value stays visible:

| Quantity | Submitted | Revised |
|---|---|---|
| gene-set permutation empirical p | 0.005, from 1,000 unmatched 22-gene sets | **0.0002**, from 10,000 sets matched to the YAP/TAZ signature for expression level and detection rate |
| partial rho after p-EMT adjustment | +0.064, p = 2.6e-03 | **0.040**, p = 0.06 |

Both revised values come from
`revision_outputs/gse103322_puram/P5_YAP_TAZ_permutation.csv`, which also gives
the matched-null mean (0.0399), the 95th percentile of |rho| (0.0966), the
permutation resolution (1e-04) and the unmatched-null result (0.0003). The
change is described in the response to reviewers, which quotes the submitted
and revised passages side by side.

The related row `WBP5 vs YAP/TAZ target signature, rho=0.17 p=4.2e-16` is
unchanged and remains correct; the revised manuscript reports it as rho = 0.173.

## Everything else

The remaining rows were checked against the revised manuscript during the
November 2026 verification sweep and still hold. That sweep also compared every
table cell in the manuscript and the Supplementary Information against the
deposited source files, and every statistic quoted in the running text of the
manuscript, the Supplementary Information and the response to reviewers. Where a row and the revised
manuscript differ in rounding only, the manuscript value is the one to use.
