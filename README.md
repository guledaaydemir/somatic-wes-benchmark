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

Run them in the order below. Notebook 01 needs the raw VCFs in `data/raw/snv/` (see `data/README.md`); notebook 00 runs without them.

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
   | 1 | `00_data_preparation` | the committed files in `data/truth/`, `data/reference/`, `data/derived/`; `data/raw/` is optional (listed and hashed if present) | `results/00_data_preparation/` (file inventory with SHA-256, BED territory, truth-VCF census) |
   | 2 | `01_variant_sets` | `data/raw/snv/` (480 run folders; parses every VCF), `data/reference/TestCases.csv` | `data/derived/sets_dict.parquet`, `results/01_variant_sets/` |
   | 3 | `02_pairwise_iou` | notebook 01's cache, `data/reference/TestCases.csv` | `results/02_pairwise_iou/` |
   | 4 | `02_pairwise_iou_indels` | `data/derived/sets_dict_indels.csv`, `data/reference/TestCases_Indels.csv`, `data/derived/legacy_indel_iou_recorded.csv` | `results/02_pairwise_iou_indels/` |
   | 5 | `03_execution_time` | notebook 01's cache (run keys only), `data/reference/TestCases.csv` | `results/03_execution_time/` |
   | 6 | `04_tmb` | notebook 01's cache, `data/reference/TestCases.csv`, `data/reference/sorted_exome_hc.bed.gz`, `data/reference/sequencing_yield.csv` | `results/04_tmb/` |
   | 7 | `05_metrics_anova` | notebook 01's cache, `data/truth/hc_bed_filtered.recode.vcf`, `data/reference/TestCases.csv`, `data/reference/sequencing_yield.csv` | `results/05_metrics_anova/` |
   | 8 | `06_validated_vcfs` | the committed validated VCFs in `data/reference/validated_vcfs/` (45 files, parsed in about ten seconds), `data/truth/hc_bed_filtered.recode.vcf` | `results/06_validated_vcfs/` |
   | 9 | `08_consensus` | notebook 01's cache, `data/truth/hc_bed_filtered.recode.vcf`, `data/reference/TestCases.csv` | `results/08_consensus/` |
   | 10 | `09_region_complexity` | notebook 01's cache, `data/truth/hc_bed_filtered.recode.vcf`, `data/reference/stratification/*.bed.gz`, `data/reference/TestCases.csv` | `results/09_region_complexity/` |
   | 11 | `07_stratification` | notebook 01's cache, `results/00_data_preparation/tables/run_manifest.csv` (committed; no notebook regenerates it any more), `data/truth/hc_bed_filtered.recode.vcf`, `data/reference/stratification/*.bed.gz` | `results/07_stratification/` |
   | 12 | `10_ensemble` | notebook 01's cache, `data/truth/hc_bed_filtered.recode.vcf`, `data/reference/TestCases.csv` | `results/10_ensemble/` |

   Notebooks 02 to 05, 07, 08 and 09 read the cache written by 01, so a later notebook never re-parses the run VCFs. Notebook 06 does not use the cache: it parses the committed validated VCFs. Every notebook also reads the legacy caches in `data/derived/` for its parity check.

**Parity.** The last code cells check the notebook's outputs against the legacy values. If a check fails the notebook raises `AssertionError`, only `audit/parity.txt` (marked FAIL) is written, and no tables, figures or cache appear. On success the outputs are moved into place and `audit/summary.txt` is written. Running a notebook again overwrites its outputs with identical files.

**Tests.** From the repository root, with the environment active: `pytest`.
