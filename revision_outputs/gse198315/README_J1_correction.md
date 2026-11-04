# Correction to `J1_Fig4A_proxy_compartment_composition.csv`

**File:** `revision_outputs/gse198315/J1_Fig4A_proxy_compartment_composition.csv`

## What was wrong

The deposited version listed 2,344 WBP5-high cells. The manuscript, the
Supplementary Information and Fig. 4 all use 2,345. One cell was lost when this
verification table was exported on 2026-10-24; the EGFR-high (878) and
double-high (84) totals were unaffected.

The missing cell is a WBP5-only cell in the tumor-core region, assigned to the
Endothelial-like proxy compartment. Four independent lines of evidence place it
there, and no other position is consistent with all four:

1. **The two chi-square tests of Fig. 4(D).** Computed over the three mutually
   exclusive groups named in Methods Section 2.2 and Supplementary Table S3
   (WBP5-only 2,261; EGFR-only 794; double-high 84), restoring this one cell
   gives chi-square 26.79 on 6 df with 1 of 12 expected counts below five, and
   692.47 on 8 df with 2 of 15 — matching the published values in every field.
   Without it the statistics are 26.87 and 691.95.
2. **Fig. 4A.** The printed WBP5-high tumor-core row reads 67.4 and 9.3 for the
   Stromal-like and Endothelial-like compartments, which require 1,056 cells
   with 98 endothelial. The uncorrected file gives 67.5 and 9.2.
3. **Supplementary Table S3(e).** The WBP5-high Endothelial-like share of
   11.51% requires 270 of 2,345. The uncorrected file gives 269 of 2,344,
   which is 11.48%.
4. **The group total itself**, which the Fig. 4 legend and
   `make_fig4c_runB.py` both state as 2,345.

A fifth line of reasoning was considered and discarded: the WBP5-only
tumor-core share of 44.6% in Table S3(a) does not discriminate, because
1,007/2,260 and 1,008/2,261 both round to 44.6%.

## What has been changed

One row. `WBP5-high,TC`:

| Field | Before | After |
|---|---|---|
| n_cells_in_group_and_region | 1055 | **1056** |
| n_Endothelial-like | 97 | **98** |
| pct_Stromal-like | 67.5 | **67.4** |
| pct_Endothelial-like | 9.2 | **9.3** |

The percentage columns were recomputed from the counts in the same row; the
other three (19.3, 3.0, 0.9) are unchanged by the correction. Group totals are
now 2,345 / 878 / 84. No other row, column or value was touched.

## Reproduce

```
python -c "
import csv
d=list(csv.DictReader(open('J1_Fig4A_proxy_compartment_composition.csv')))
t={}
for r in d: t[r['group']]=t.get(r['group'],0)+int(r['n_cells_in_group_and_region'])
print(t)
"
```

should print `{'WBP5-high': 2345, 'EGFR-high': 878, 'Double-high': 84}`.

The chi-square check in point 1 above uses these counts with the double-high
row subtracted from the two "-high" rows, which is how Methods Section 2.2
defines the mutually exclusive groups. That subtraction matters: carried out on
the overlapping WBP5-high and EGFR-high groups instead, the compartment test
returns 652.85 rather than the 692.47 the Fig. 4(D) legend reports, on the same
8 degrees of freedom and with the same 2 of 15 expected counts below five. The
degrees of freedom and the expected-count check therefore do not distinguish
the two groupings, and only the statistic does.
