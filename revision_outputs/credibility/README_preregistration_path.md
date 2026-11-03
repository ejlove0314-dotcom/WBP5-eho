# The `preregistration` line in `CR_run_info.txt`

Each of the four run records in this directory carries the line

    preregistration claude/CSBJ_신뢰도분석4종_사전규칙_20261021.md

That path does not resolve in this repository. It is the working location of
the pre-analysis plan on the machine the analyses were run on, written into the
record by `revision_outputs/scripts/CRED_profiles_20261021.py` at run time.

**The document it names is deposited here as
`reproducibility/pre_analysis_plan_20261021_KO.md`**, with an English
translation beside it as `pre_analysis_plan_20261021_EN.md`. The Korean file is
the original and the authoritative record; the response letter cites the
English translation. The deposited Korean file is that same document: its
common definitions, its four sets of decision criteria, its pre-specified
gene lists and its negative-control list are the plan the runs recorded here
were judged against.

The hash is computed over the file with LF line endings, as it is stored
in the repository; a checkout that converts line endings will not reproduce
it. The repository carries a `.gitattributes` that fixes this.

The line itself has deliberately not been rewritten, for the reason given under
"Running the deposited scripts" in the top-level `README.md`: these records are
the record of what was executed. Line 1 of each record states the SHA-256 of
the script that produced it, and that hash matches the deposited script as it
stands; editing either the script or the record would break a check a reader
can run. The reference is therefore corrected here rather than in place.

No value reported in the manuscript, the Supplementary Information or the
response letter depends on this line.

## Where each file is

| Named in the record | Deposited as |
|---|---|
| `claude/CSBJ_신뢰도분석4종_사전규칙_20261021.md` | `reproducibility/pre_analysis_plan_20261021_KO.md` |
| `results_cred/` (the run-time output directory) | `revision_outputs/credibility/` |

The second row applies to the `--outdir results_cred\...` arguments in the
usage examples of `CRED_profiles_20261021.py` and
`CRED_crosscohort_20261021.py`, which likewise record the layout used at run
time. The two mentions of `results_cred/` in the pre-analysis plan itself were
in the note added when the plan was deposited and in the section recorded after
execution, not in the pre-registered criteria, and both now name this
directory.
