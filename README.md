# somatic-wes-benchmark
Analysis code for the SEQC2 somatic WES study: 480 pipeline configurations (SARIYER/ALTAY; trimming, base recalibration, duplicate handling; BWA-MEM/Bowtie2; MuTect2/Strelka2/SomaticSniper). Computes precision/recall/F1-score, execution times, ANOVA, TMB, and consensus/ensemble results.

## Layout
- `notebooks/` — one notebook per analysis (`NN_name.ipynb`)
- `results/NN_name/` — `tables/`, `figures/`, `audit/` for the matching notebook
- `src/swb/` — shared code (`config`, `io`, `metrics`, `viz`)
- `data/derived/` — cached intermediates; `data/reference/` — truth-set and region files
- `docs/` — `legacy-deviations.md`
- `legacy/` — frozen original notebooks

Working rules are in `CLAUDE.md`. Raw VCFs live on an external drive; set `SWB_DATA_ROOT` to override the default location.
