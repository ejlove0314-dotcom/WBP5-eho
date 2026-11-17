#!/usr/bin/env python3
"""
make_full_from_annotated.py — rebuild "<name>_full_annotated.h5ad" from "<name>_annotated.h5ad".

Step 2 of GSE198315_reannotation_20261015.py writes two files from the same object:
  <name>_annotated.h5ad       highly variable genes, with the full log-normalized matrix kept in .raw
  <name>_full_annotated.h5ad  = adata.raw.to_adata() plus the cell_type labels
so the second can be recovered from the first without re-running anything.

Usage:
  python make_full_from_annotated.py gse198315_checkpoint_annotated.h5ad
"""
import sys
import anndata as ad

src = sys.argv[1]
if not src.endswith("_annotated.h5ad") or src.endswith("_full_annotated.h5ad"):
    raise SystemExit("give the <name>_annotated.h5ad file")
dst = src.replace("_annotated.h5ad", "_full_annotated.h5ad")
a = ad.read_h5ad(src)
if a.raw is None:
    raise SystemExit("ERROR: this file has no .raw matrix; it cannot be used")
full = a.raw.to_adata()
for c in ("cell_type", "region", "patient", "total_counts"):
    if c not in full.obs:
        if c not in a.obs:
            raise SystemExit(f"ERROR: column '{c}' is missing")
        full.obs[c] = a.obs[c].values
print(f"{src}: {a.shape} -> full matrix {full.shape}")
print(full.obs["cell_type"].value_counts().to_string())
full.write_h5ad(dst, compression="gzip")
print("written:", dst)
