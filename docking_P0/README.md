# MTiOpenScreen library screen at pocket P0

Library-scale docking of the MTiOpenScreen Drugs-lib collection against the same pocket,
receptor file, grid definition and exhaustiveness used in the local pipeline. This screen
provides the external reference scale reported in Section 3.5 of the manuscript.

## Job record

| Parameter | Value |
|-----------|-------|
| Server | MTiOpenScreen (RPBS, Universite Paris Cite) |
| Job ID | V06821248773098 |
| Executed | 13 August 2026 |
| Receptor | `WBP5_AlphaFold.pdb` (unmodified AF-Q9UHQ7-F1-model_v6) |
| Library | Drugs-lib; 7,173 stereoisomers representing 4,574 unique drugs |
| Ligand limit | 10,000 (exceeds library size, so the whole library was screened) |
| Search space | custom mode; center -23.835, -0.451, 23.910; box 25 x 25 x 25 A |
| Grid spacing | 1 A (server default) |
| Exhaustiveness | 8 |

## Files

| File | Contents |
|------|----------|
| `output.table.csv` | Ranked output table returned by the server (1,500 top-scoring entries). |
| `index.xml` | Server job record, including the command line and all input parameters. |

## Summary statistics

Computed directly from `output.table.csv`:

| Quantity | Value |
|----------|-------|
| Entries returned | 1,500 (20.9% of the 7,173-entry library) |
| Distinct drugs among them | 1,037 |
| Energy range | -8.7 to -5.7 kcal/mol |
| Mean +/- SD | -6.17 +/- 0.38 kcal/mol |
| Pazopanib | -5.8 kcal/mol |
| Entries scoring more favorably than pazopanib | 1,180 |
| Entries tied with pazopanib | 257 |
| Pazopanib rank band | 1,181st-1,437th of 7,173 (top 16.5-20.0%) |

## Comparability with the local pipeline

Server and local energies are not interchangeable. Pazopanib scored -5.8 kcal/mol here
against -6.86 +/- 0.05 kcal/mol under the local protocol, a 1.06 kcal/mol offset arising
from the independent ligand preparation used by each pipeline; the receptor file, grid
definition and exhaustiveness were identical. Rankings from this screen are interpreted
only within the server scale.

The server returns only its top-scoring entries, so all library members scoring below
-5.7 kcal/mol - 79.1% of the collection - are unobserved.

## Superseded run

An earlier screen (job P12689576379061, 10 August 2026) was run at grid center
`2.402, -3.010, 1.561`, which is pocket P_3 rather than P0. Its output is not included
here and none of its values appear in the manuscript.
