# Correction to the per-residue RSA deposit

**File affected:** `revision_outputs/structure/WBP5_RSA_per_residue.csv` and
`revision_outputs/structure/compute_rsa_shrake_rupley.py`

## What was wrong

The deposited script computed solvent-accessible surface area with
`ShrakeRupley(n_points=960)`, and the deposited CSV was the output of that
960-point run. Supplementary Table S18 reports the values obtained with
Biopython's default 100-point sphere sampling, which is what its legend states.
A reader who ran the deposited script therefore obtained 23 relative solvent
accessibility values that differed from Table S18 by up to about 0.03.

No value in the manuscript or the Supplementary Information was incorrect. The
discrepancy was confined to the deposited script and CSV.

## What has been changed

- `compute_rsa_shrake_rupley.py` now uses the default 100-point sampling, takes
  the sampling density as an optional third argument, and prints the P0 summary
  so that agreement with Table S18 can be checked directly.
- `WBP5_RSA_per_residue.csv` has been regenerated at 100 points. All 23
  residues listed in Table S18 now reproduce exactly, as do the pLDDT values.
- `WBP5_RSA_per_residue_960pt.csv` is retained as a companion file so that the
  earlier numbers remain available.

## Why the choice of sampling does not affect the conclusion

Over the P0 region (residues 38-58):

| Sphere sampling | RSA range | Median |
|---|---|---|
| 100 points (reported in Table S18) | 0.266 to 0.593 | 0.516 |
| 960 points | 0.290 to 0.624 | 0.521 |

Every residue of the P0 region is solvent-exposed under either setting, so the
statement made in the Table S18 legend holds in both cases.

## Reproduce

```
python compute_rsa_shrake_rupley.py AF-Q9UHQ7-F1-model_v6.pdb WBP5_RSA_per_residue.csv
python compute_rsa_shrake_rupley.py AF-Q9UHQ7-F1-model_v6.pdb WBP5_RSA_per_residue_960pt.csv 960
```

Requires Biopython. Verified with Biopython 1.88.
