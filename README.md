# somatic-wes-benchmark
Analysis code for the SEQC2 somatic WES study: 480 pipeline configurations (SARIYER/ALTAY; trimming, base recalibration, duplicate handling; BWA-MEM/Bowtie2; MuTect2/Strelka2/SomaticSniper). Computes precision/recall/F1-score, execution times, ANOVA, TMB, and consensus/ensemble results.

## Layout
- `notebooks/` — one notebook per analysis (`NN_name.ipynb`)
- `results/NN_name/` — `tables/`, `figures/`, `audit/` for the matching notebook
- `src/swb/` — shared code (`audit`, `config`, `io`, `metrics`, `stats`, `viz`)
- `data/` — `raw/` (per-run VCFs, local only), `truth/`, `reference/` (metadata, BED and stratification files, `.bed.gz`), `derived/` (caches); see `data/README.md`
- `tests/` — pytest for `src/swb/` (`pytest` from the repo root)
- `docs/` — `legacy-deviations.md`
- `legacy/` — frozen original notebooks

Working rules are in `CLAUDE.md`. The repository is self-contained: every path is relative to the repository root and set in `src/swb/config.py`.

## Running the notebooks

Run them in the order below. Notebooks 00 and 01 need the raw VCFs in `data/raw/snv/` (see `data/README.md`).

1. **Environment.** Create it once from `environment.yml` (or reuse `mywork2_rebuild_env` if it already exists), then activate it:
   ```bash
   conda env create -f environment.yml
   conda activate mywork2_rebuild_env
   ```
2. **Data.** Place the raw SNV VCFs in `data/raw/snv/` as described in `data/README.md`. Everything else the notebooks read is committed.
3. **Kernel.** Register the environment as a Jupyter kernel (safe to repeat):
   ```bash
   python -m ipykernel install --user --name mywork2_rebuild_env --display-name "Python (mywork2_rebuild_env)"
   ```
   The notebooks are saved with the generic `python3` kernel. Select `Python (mywork2_rebuild_env)` when Jupyter asks, or start Jupyter from the activated environment so that `python3` resolves to it.
4. **Start Jupyter from `notebooks/`.** The first cell of each notebook falls back to `../src` if `swb` is not installed, which needs that working directory:
   ```bash
   cd notebooks
   jupyter lab
   ```
5. **Run the notebooks top to bottom, in this order:**

   | Order | Notebook | Needs | Writes |
   |---|---|---|---|
   | 1 | `00_data_preparation` | `data/raw/snv/`, `data/reference/TestCases.csv`, `data/truth/hc_bed_filtered.recode.vcf`, `data/reference/sorted_exome_hc.bed.gz` | `results/00_data_preparation/` (`run_manifest.csv`) |
   | 2 | `01_variant_sets` | notebook 00's manifest; parses all 480 VCFs on the drive | `data/derived/sets_dict.parquet`, `results/01_variant_sets/` |
   | 3 | `07_stratification` | notebooks 00 and 01, `data/truth/hc_bed_filtered.recode.vcf`, `data/reference/stratification/*.bed.gz` | `results/07_stratification/` |

   Notebooks 02 to 06 are not written yet. Each notebook needs the one before it in the table, and 07 reads the cache written by 01, so a later notebook never re-parses the run VCFs. Every notebook also reads the legacy caches in `data/derived/` for its parity check.

**Parity.** The last code cells check the notebook's outputs against the legacy values. If a check fails the notebook raises `AssertionError`, only `audit/parity.txt` (marked FAIL) is written, and no tables, figures or cache appear. On success the outputs are moved into place and `audit/summary.txt` is written. Running a notebook again overwrites its outputs with identical files.

**Tests.** From the repository root, with the environment active: `pytest`.
