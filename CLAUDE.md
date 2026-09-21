# somatic-wes-benchmark — working agreement

## What this repository is
A reorganisation of a single large exploratory notebook (`legacy/StabilityAnalysis.ipynb`)
into one notebook per analysis, so that a reader can follow the study from raw data to
published figure. The science is finished; this is a restructuring and verification task.

## Hard rules

1. **`legacy/` is frozen.** Never edit, never execute, never import from it. It is the
   reference implementation used to verify migrated code.

2. **Never invent a number.** Every value in a table or figure must come from executed
   code in this repository. If a legacy value cannot be reproduced, say so in the
   notebook's audit output and stop — do not substitute an approximation.

3. **Parity before replacement.** A migrated notebook is not finished until its outputs
   match the legacy values. Each notebook ends with a parity cell that asserts this and
   writes `results/NN_name/audit/parity.txt`. If parity fails, the notebook must fail
   loudly rather than write results.

4. **No silent changes to analysis logic.** Variant filtering, key representation
   (`CHROM_POS_REF_ALT1`), truth-set construction, pipeline inclusion, and metric
   definitions are fixed. If a legacy behaviour looks like a bug, write it in
   `docs/legacy-deviations.md` with the evidence and ask — do not fix it silently.

5. **One notebook, one analysis, one results folder.** `results/NN_name/` holds
   `tables/`, `figures/`, `audit/`. Nothing writes outside its own folder.

6. **Shared code goes in `src/swb/`.** If two notebooks need the same function, it moves
   to `src/` and both import it. No copy-paste between notebooks.

7. **Expensive work is cached.** Notebook 01 parses the 480 VCFs once and writes
   `data/derived/sets_dict.parquet`. Downstream notebooks read the cache and must never
   re-parse VCFs unless the notebook's stated purpose requires raw VCF fields.

8. **Deterministic output.** Fixed random seeds. Explicit `sort_values` before every
   `to_csv`. `index=False`. Running a notebook twice produces byte-identical files.

## Known data hazards

- **macOS AppleDouble files.** `._filename` files sit beside real files on the external
  drive and parse as empty. Every glob must exclude names starting with `._`.
  Use `swb.io.real_files()`.
- **Non-UTF-8 bytes** appear in some VCF headers. Open with `errors="replace"`.
- **Tumour/normal order.** Read `##tumor_sample=` from the VCF header. The FILENAME LISTS
  NORMAL FIRST — any filename-based inference is inverted.
- **FILTER is not plain `PASS`.** The truth set uses `HighConf;PASS` and `MedConf;PASS`.
  Test with `"PASS" in f[6].split(";")`.
- **SomaticSniper has no PASS.** Its FILTER is `.` on every record; the legacy filter
  exempts it. Preserve this and state the asymmetry in any output that uses it.
- **`TestCaseNo` format differs between files** — `"TestCase 269"` in metadata,
  `"269"` in cached sets. Normalise via `swb.io.canon()`.
- **`TestCases_Indels.csv` has a different schema** from `TestCases.csv`: 13 columns, no
  `Environment`, 320 runs rather than 480.

## Notebook template
Every notebook has these sections, in order:
1. Purpose — what question, which reviewer comment if any, what it does NOT do
2. Imports — `from swb import io, metrics, viz`
3. Config — paths from `swb.config`, no hard-coded absolute paths
4. Load + validate — assert row counts, factor levels, no missing values
5. Analysis
6. Save tables
7. Figures
8. **Parity check** — assert against legacy values
9. Audit summary written to `results/NN_name/audit/`

## Style
- Terse, targeted edits. No explanatory prose in commit messages or code comments beyond
  what a reader needs.
- Do not refactor code you were not asked to touch.
- Ask before deviating from these rules.

## Data and paths — THE REPOSITORY IS SELF-CONTAINED

The repository must run on a reviewer's machine with nothing outside it.

- **Never write an absolute path.** No `/Volumes/`, no `/Users/`, no home directory.
- Every path is relative to the repository root, resolved in `src/swb/config.py`:

      REPO_ROOT = Path(__file__).resolve().parents[2]
      DATA      = REPO_ROOT / "data"
      RAW_SNV   = DATA / "raw" / "snv"        # 480 per-run SNV VCFs  (gitignored)
      RAW_INDEL = DATA / "raw" / "indel"      # 320 per-run indel VCFs (gitignored)
      TRUTH     = DATA / "truth"
      REFERENCE = DATA / "reference"
      DERIVED   = DATA / "derived"
      RESULTS   = REPO_ROOT / "results"

- `SWB_DATA_ROOT` is REMOVED. Delete every reference to it.
- The BED and stratification files are gzipped (`.bed.gz`). Read them with gzip.
- Two tiers:
  * Committed: `data/truth/`, `data/reference/`, small CSVs in `data/derived/`.
  * Local only: `data/raw/` (raw VCFs) and `data/derived/*.parquet` caches.
- **Notebooks 02 onward must run from the committed cache alone.** Only notebook 01
  (which builds the cache) and analyses that genuinely need raw VCF fields may read
  `data/raw/`. Any notebook reading `data/raw/` states this in its Purpose section and
  checks that the directory exists, failing with a clear message pointing to
  `data/README.md` if it does not.

## Execution
Notebooks are executed by the user in Jupyter. Write code and parity assertions; do not
attempt to run notebooks that read `data/raw/`.

## Environment
Python 3.8, pandas 1.x, conda env `mywork2_rebuild_env`. No pandas 2.x-only idioms.