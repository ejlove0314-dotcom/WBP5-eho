#!/usr/bin/env python3
"""
Repair reproducibility/GSE103322_WBP5_cell_level_raw.csv for re-deposition.

Three file-side defects, none of which affects any value reported in the
manuscript or the supplementary information:

(1) 18 rows carry the literal string "#NAME?" in `non_cancer_type` and
    `compartment`.  The underlying GSE103322 annotation for this category is
    "-Fibroblast"; because the label begins with a hyphen, a spreadsheet
    application parses it as a formula and writes back #NAME?.  The category is
    the one tabulated in Table S6 as "Fibroblast (excl.)" (n = 18 before the
    pooled-plate exclusion, n = 17 after).

(2) 173 rows carry an empty `patient`.  Their cell identifiers use the form
    HNSCC_17_P10_G09_S273_comb, with an underscore after the cohort prefix,
    whereas the remaining cells of the same patient use HNSCC17_...  The
    identifier parser matched only the second form, so these cells -- all
    T cells of patient 17 -- lost their patient assignment when the file was
    written.  The analysis itself retained them: with these cells included the
    malignant-versus-pooled-non-malignant paired comparison returns the
    published mean difference of 0.79 log2 units, d_z = 0.73 and Wilcoxon
    p = 0.049, whereas without them it returns 0.72.

(3) A further 109 rows carry an empty `patient` by design: they are the pooled
    "combo" plate cells that Methods Section 2.2 excludes.  These are left
    unassigned and are now flagged explicitly in a new `pooled_plate` column so
    that the exclusion is reproducible from the file alone.

`compartment` and `non_cancer_type` are additionally normalized to the labels
used in Table S6, so that the table can be reproduced from this file without an
external mapping and so that no label begins with a character a spreadsheet
treats as an operator.  The one GEO annotation string that this renames,
"-Fibroblast" -> "Fibroblast (excl.)", is recorded in the README rather than in
a column of its own, because a column holding the original string would
reintroduce the same failure mode.
"""
import csv
import re
import sys

SRC = sys.argv[1]
DST = sys.argv[2]

# GEO annotation string -> Table S6 compartment label
NORM = {
    "fibroblast":  "Fibroblast",
    "-fibroblast": "Fibroblast (excl.)",
    "endothelial": "Endothelial",
    "macrophage":  "Macrophage",
    "mast":        "Mast",
    "b cell":      "B cell",
    "dendritic":   "Dendritic",
    "t cell":      "T cell",
    "myocyte":     "Myocyte",
}
UNANNOTATED = {"", "0", "nan", "na"}

PT17 = re.compile(r"^HNSCC_17_")
COMBO = re.compile(r"^HNSCC_combo", re.I)

rows = list(csv.DictReader(open(SRC, newline="")))
fields = list(rows[0])
out_fields = fields + ["pooled_plate"]

counts = dict(name_fixed=0, pt17_filled=0, combo_flagged=0, unannotated=0, myocyte=0)

for r in rows:
    # --- (1) restore the GEO label that the spreadsheet corrupted ----------
    raw = r["non_cancer_type"]
    if raw == "#NAME?":
        raw = "-Fibroblast"
        counts["name_fixed"] += 1

    # --- compartment, normalized to the Table S6 labels -------------------
    if r["is_cancer"] == "1":
        r["compartment"] = "Cancer cells"
        r["non_cancer_type"] = ""
    else:
        key = raw.strip().lower()
        if key in UNANNOTATED:
            label = "Unannotated"
            counts["unannotated"] += 1
        else:
            label = NORM.get(key, raw.strip())
            if key == "myocyte":
                counts["myocyte"] += 1
        r["compartment"] = label
        r["non_cancer_type"] = label

    # --- (2) and (3) the two kinds of empty patient field -----------------
    combo = bool(COMBO.match(r["cell_id"]))
    r["pooled_plate"] = "TRUE" if combo else "FALSE"
    if combo:
        counts["combo_flagged"] += 1
    elif not r["patient"].strip() and PT17.match(r["cell_id"]):
        r["patient"] = "17"
        counts["pt17_filled"] += 1

still_blank = [r for r in rows
               if not r["patient"].strip() and r["pooled_plate"] == "FALSE"]
assert not still_blank, f"{len(still_blank)} rows still have no patient"
assert counts["name_fixed"] == 18, counts
assert counts["pt17_filled"] == 173, counts
assert counts["combo_flagged"] == 109, counts

with open(DST, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=out_fields)
    w.writeheader()
    w.writerows(rows)

for k, v in counts.items():
    print(f"  {k:16s} {v}")
print(f"\nwrote {DST}  ({len(rows)} rows, {len(out_fields)} columns)")
