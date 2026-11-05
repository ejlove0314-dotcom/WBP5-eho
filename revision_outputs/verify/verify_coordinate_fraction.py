#!/usr/bin/env python3
"""Re-derive the genomic-coordinate assignment rate reported for the GSE198315
inferCNV run (Methods Section 2.6; response to reviewers).

The step-6 script logs this fraction to stderr at run time but does not save it.
This script recomputes it from the same two inputs, using the identical gene-name
matching rule, and prints:

    genes_with_coordinates 20788 / 33694 (61.7%)

Usage
-----
    python verify_coordinate_fraction.py \
        --h5ad gse198315_checkpoint_full_annotated.h5ad \
        --gtf  gencode.v50.basic.annotation.gtf.gz

Requires only h5py; neither scanpy nor infercnvpy is needed.
"""
import argparse
import gzip

import h5py


def gtf_gene_names(gtf_path):
    """Gene symbols carrying usable coordinates, by the rule used in step 6."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    names = set()
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            i = cols[8].find('gene_name "')
            if i < 0:
                continue
            try:
                int(cols[3]), int(cols[4])
            except ValueError:
                continue
            names.add(cols[8][i + 11:cols[8].find('"', i + 11)])
    if not names:
        raise SystemExit("error: no 'gene' rows with gene_name found in the GTF")
    return names


def var_names(h5ad_path):
    """The var index of an .h5ad, for plain, categorical or nullable-string storage."""
    with h5py.File(h5ad_path, "r") as fh:
        var = fh["var"]
        obj = var[var.attrs["_index"]]
        if not isinstance(obj, h5py.Group):
            return [x.decode() if isinstance(x, bytes) else str(x) for x in obj[:]]
        if "values" in obj:                      # nullable-string-array
            values = obj["values"][:]
            mask = obj["mask"][:] if "mask" in obj else None
            out = []
            for i, x in enumerate(values):
                if mask is not None and mask[i]:
                    out.append(None)
                else:
                    out.append(x.decode() if isinstance(x, bytes) else str(x))
            return out
        categories = [                            # categorical
            x.decode() if isinstance(x, bytes) else str(x)
            for x in obj["categories"][:]
        ]
        return [categories[c] if c >= 0 else None for c in obj["codes"][:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="step-2 output (…_full_annotated.h5ad)")
    ap.add_argument("--gtf", required=True, help="GENCODE annotation used by step 6")
    args = ap.parse_args()

    coords = gtf_gene_names(args.gtf)
    genes = var_names(args.h5ad)
    n = len(genes)
    assigned = sum(1 for g in genes if g in coords)
    print(f"gtf_unique_gene_names {len(coords)}")
    print(f"genes_with_coordinates {assigned} / {n} ({100.0 * assigned / n:.1f}%)")
    print(f"  exact {100.0 * assigned / n:.4f}%")


if __name__ == "__main__":
    main()
