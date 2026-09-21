# Legacy deviations

Legacy behaviours that look like bugs. Each entry: what the legacy code does, the evidence, and the decision. `src/swb/` reproduces the legacy behaviour unless an entry says it was decided otherwise (D1, D4, D5, D6, D7); nothing else is changed without an explicit decision.
Cell numbers are 0-based indices over all cells of `legacy/StabilityAnalysis.ipynb`.

## D1. BED filter ignores chromosomes
Status: **decided.** Both filters are implemented. **The corrected filter changes the region-stratified results substantially; published region-stratified values require revision.**

`filter_vcf_with_bed_single_chrom` only restricts by chromosome when `chrom` is given. Every legacy call (cells 215, 226, 237) omits it, so each variant's POS is compared against the regions of *all* chromosomes: a chr2 variant at POS 15 passes a chr1-only region `10-20`.

Implementation:
- `swb.io.bed_filter_legacy` reproduces the chromosome-blind behaviour exactly. It is used only in parity cells, and in notebook 07 to fill the comparison table. It gives the same rows as the legacy two-pointer sweep (kept as the reference in `tests/test_bed_filter.py`; equal on 25 random cases and on real data) without walking every region, which took 36 s per call on the largest BED.
- `swb.io.bed_filter` matches a variant only against regions of its own chromosome. It is the default for all new analysis. Its region boundary is still the legacy one (D2).
- Notebook 07 computes both and writes `results/07_stratification/tables/legacy_vs_corrected.csv` (per-region counts under each) and `stratified_metrics.csv` (per run and region, both filters).

Reproduction of the published values: the legacy filter reproduces all 2880 (run, region) rows of `data/derived/union_metrics.csv` (IoU, precision, recall, F1) to within 2.3e-16, so the published region-stratified numbers came from the chromosome-blind logic.

Effect, 480 runs pooled (692,587 variants), from `data/derived/sets_dict.csv` and the six GIAB BEDs, computed with the two filters above. Notebook 07's `legacy_vs_corrected.csv` is the authoritative version of this table:

| region | kept, legacy | kept, corrected | fraction, legacy | fraction, corrected |
|---|---:|---:|---:|---:|
| alldifficultregions | 640,231 | 206,951 | 0.924 | 0.299 |
| notinalldifficultregions | 690,451 | 485,636 | 0.997 | 0.701 |
| alllowmapandsegdupregions | 523,227 | 27,551 | 0.755 | 0.040 |
| notinalllowmapandsegdupregions | 692,563 | 665,036 | 1.000 | 0.960 |
| segdups | 419,479 | 25,867 | 0.606 | 0.037 |
| notinsegdups | 692,571 | 666,720 | 1.000 | 0.963 |

Under the legacy filter 61% of all variants fall in segmental duplications; under the corrected filter 3.7% do. Mean F1 per run in `segdups` goes from 0.519 to 0.033 and in `alllowmapandsegdupregions` from 0.590 to 0.037 (D9 applies to both). Under the corrected filter 4 runs (140, 286, 428, 574) keep no variant in `segdups` and in `alllowmapandsegdupregions`, so precision and F1 are undefined there; legacy code would raise `ZeroDivisionError`, and notebook 07 leaves them blank.

Tests: `test_legacy_ignores_chromosome_when_chrom_is_none` pins the legacy behaviour; `test_corrected_respects_chromosome` and `test_corrected_matches_one_legacy_call_per_chromosome` pin the corrected one.

## D2. BED boundary convention
Status: **pending decision** (found by reading code; effect not measured)

The filter keeps `start <= POS < end`. BED is 0-based half-open and VCF POS is 1-based, so the matching test is `start < POS <= end`. Legacy includes the base just before each region and drops the last base of it (cells 213, 226, 238). Both `bed_filter_legacy` and `bed_filter` keep this convention, so the D1 comparison isolates the chromosome handling.

## D3. `filtering_df` has two identical columns
Status: **decided**: ported as-is, documented, never read as a filtering cascade.

Cell 32 appends `(key, total_len, len(temp), without_indels_len)` after both filters have run, so `After FilterPass` and `After RemoveIndels` are both the count after every filter (PASS, `is_snp`, single-base REF and ALT_1). `After FilterPass` is mislabelled: it is not an intermediate count. Verified in `data/derived/filtering_df.csv`: equal in 480 of 480 rows. The "Changes in VCF Lengths (SNPS)" figure shows the same box twice.

Handling: `swb.io.parse_vcf` returns records read and records kept, and `01_variant_sets` writes exactly those two (`records_read`, `records_kept`). Its parity cell checks `records_kept` against both legacy columns and labels them as not a cascade. Nothing new should treat those two columns as a step-by-step filtering count.

## D4. `anova_lm(model, type=7)` silently fell back to Type I
Status: **decided**: Type II in all new code.

`anova_lm` reads `typ`, not `type` (statsmodels 0.14.1: `typ = kwargs.get('typ', 1)`), so `type=7` was ignored and the legacy table (cell 133) is Type I (sequential) sums of squares. The design is fully balanced (5 × 2 × 2 × 2 × 2 × 3 × 2 factor levels, one run per combination), where Type I, II and III coincide, so the published variance percentages are unaffected. Checked on the real data (run-level F1, 480 runs, 127 terms): Type II and Type I sums of squares agree to 1.2e-14, and `variance_pct` returns the same 10 sources in the same order, differing by at most 7e-14 percentage points (VariantCaller 53.9%, Center 19.4%, Mapper:VariantCaller 8.0%).

Implementation: `swb.stats.fit_anova` defaults to `typ=2`; `typ=1` is kept for reproducing legacy. With one run per cell the model is saturated (residual df 0), and statsmodels cannot form Type II there (it divides 0 by 0), so `swb.stats.type2_sum_sq` computes Type II by nested-model residual sums of squares. It matches statsmodels' `typ=2` to 2.6e-14 on unbalanced data with every cell populated (where Type I differs by up to 2.1), and `fit_anova` raises if a design cell is empty. The legacy formula equals `full_factorial_formula("f1_score")` term for term, in the same order (127 terms, checked).

## D5. Tumour sample: two legacy functions disagreed
Status: **decided**: the `##tumor_sample=` header is authoritative; the legacy functions are deleted, not ported.

`extract_tumor_sample_name` returned the first sample column and `_get_disease_sample` the last (cell 269). Neither is right, and nothing may infer the tumour from the file name, whose accession order lists the normal first.

`swb.io.tumor_sample` is the one resolver. It reads `##tumor_sample=` (which must name a sample column); else returns `TUMOR` when the sample columns are exactly `NORMAL` and `TUMOR`; else raises `ValueError`. Nothing in `src/swb/` infers the tumour from the file name or column order.

Evidence, headers of the 476 run VCFs held in `data/raw/vcf_snps` (4 of the 480 are not in the repository copy and were not checked): 156 Mutect VCFs resolve from the header, 320 Strelka and SomaticSniper VCFs from their `NORMAL`/`TUMOR` columns, none raise. In all 156 Mutect files the file name lists the normal first and the tumour second, matching the header. Mutect column order is not fixed: tumour first in 124 files, normal first in 32, so column position is unreliable even within one caller.

## D6. "Mapper" t-test line prints baseRecalibration numbers
Status: **decided**: legacy ported as `_legacy`, corrected version written for new analysis.

Cell 128 prints `Mapper :: T-Test: t-statistic = {t_stat_br}, p-value = {p_value_ttest_br}`, which are the baseRecalibration values; the Mapper values are `t_stat_ma` / `p_value_ttest_ma`. **Any Mapper t-statistic or p-value read from that printed line is a baseRecalibration comparison and must be re-run** with `swb.stats.factor_tests`; this includes any such value in the manuscript. The repository does not show which of the two the manuscript quotes: the p-value annotated under "Mapper" in the box-plot of cell 129 comes from `p_value_ttest_ma` and is the Mapper one, so a Mapper value taken from the plot is not affected. Check each manuscript Mapper value against its source.

Implementation: `swb.stats.factor_tests_legacy` reproduces the printed text as-is, mislabel included; `swb.stats.factor_tests` returns every comparison under its own label and is the one to use. Rows are still paired by position within each group, as in legacy; whether the row order aligns was not checked, and this is unchanged.

## D7. Clustermap colours depend on `PYTHONHASHSEED`
Status: **decided**: the determinism fix is kept.

`process_dataframe` (cells 198, 200; the same pattern is in 151, 160, 184, 187) builds `lut = dict(zip(set(sample), palette))`. String set order changes between Python processes, so centres get different colours on each run. `swb.viz.cluster_pivot` orders centres alphabetically.

## D8. `clean_bed_name` has four definitions
Status: **decided**: the last definition (cell 238) is ported as `swb.viz.clean_bed_name`.

Python keeps the last definition, so cell 238 is the one in force in the notebook. The other three, kept here for reference (`\n` is a line break in the label):

| BED (`GRCh38_` prefix and `.bed` dropped) | cell 224 | cell 229 | cell 230 | cell 238 (ported) |
|---|---|---|---|---|
| segdups | `Segdups` | `Segmental Duplications` | `Segmental Duplications` | `Segdups` |
| notinsegdups | `Notinsegdups` | `Notinsegdups` (its key is misspelled `notinsegdeups`) | `Not In Segmental Duplications` | `Notinsegdups` |
| alllowmapandsegdupregions | `All Low Map \nAnd Segdup Regions` | `All Low Map\nAnd Segdup Regions` | `All Low Mapping And \nSegmental Duplication Regions` | `All Low Map\nAnd Segdup Regions` |
| notinalllowmapandsegdupregions | `Not In All Low Map \nAnd Segdup Regions` | `Not In All Low Map\nAnd Segdup Regions` | `Not In All Low Mapping And \nSegmental Duplication Regions` | `Not In All Low Map\nAnd Segdup Regions` |
| alldifficultregions | `All Difficult \nRegions` | `All Difficult\nRegions` | `All Difficult Regions` | `All Difficult\nRegions` |
| notinalldifficultregions | `Not In All \nDifficult Regions` | `Not In All\nDifficult Regions` | `Not In All\nDifficult Regions` | `Not In All\nDifficult Regions` |
| `findings from this study` | `Findings From This Study` | `Findings from\nThis Study` | `Total findings from this study` | `Findings from\nThis Study` |

Cell 224 works differently from the others: it strips `GRCh38_` and a trailing `.bed` (regex `.bed$`), replaces `_` with spaces and title-cases, then maps a few title-cased names. Cell 238 differs from 229 only by lacking the two segdups entries, so its segdups regions are labelled `Segdups` and `Notinsegdups`. It is not known which definition produced the published figures.

## D9. Region-stratified metrics use the whole truth set
Status: **pending decision** (metric definitions are fixed; not changed here)

Legacy scores each region-filtered call set against the *unfiltered* truth set (`get_precision/get_recall/get_f1score(filtered_variant_set, high_confidence)`, cells 215, 226, 237). Recall inside a region is therefore TP-in-region divided by all truth variants, and falls with the size of the region whatever the caller does. Confirmed by reproduction: computing the metrics against the unfiltered truth set (1,161 variants) reproduces all 2880 rows of `union_metrics.csv` to within 2.3e-16. Notebook 07 keeps this definition for both filters and states it in its purpose.

## D10. Support fraction to support count
Status: **decided.** Numbered D10 because D2 is the BED boundary convention (pending); this entry has no earlier number.

`swb.metrics.fraction_to_count(q, n_lists)` transfers a consensus support fraction `q` to `n_lists` call sets, for `swb.metrics.ensemble_call`: `k = ceil(q * n_lists - 1e-9)`, clipped to `[1, n_lists]`. It raises `ValueError` if `q` is not in `(0, 1]` or `n_lists < 1`. Ceil keeps the transferred rule at least as strict as the rule it was chosen under. The epsilon absorbs float error on exact ratios: `7/25 * 25` is `7.000000000000001`, and a bare ceil gives 8 instead of 7. (`34/120 * 120` is exactly `34.0` in floating point and needs no guard, but other exact ratios do.)

Values: `(34/120, 120) -> 34`, `(34/120, 96) -> 28`, `(34/120, 24) -> 7`, `(0.5, 24) -> 12`, `(1/96, 24) -> 1`, `(1.0, 24) -> 24`. `k` is non-decreasing in `q` and in `n_lists`.

Tests in `tests/test_metrics.py`: `test_fraction_to_count_known_values`, `test_fraction_to_count_non_decreasing_in_q_and_in_n_lists`, `test_fraction_to_count_equals_exact_integer_ceil_on_rational_grid`, and the rejection tests. No legacy cell computing this conversion has been traced, so there is no legacy behaviour to reproduce.

## D11. TMB territory counts overlapping BED intervals twice
Status: **pending decision** (metric definitions are fixed; the legacy value is reproduced).

`calculate_region_size` (cell 38) sums `end - start` over every line of the BED, so a base covered by two intervals counts twice. `sorted_exome_hc.bed.gz` has 272,014 intervals: raw sum 80,762,432 bp, merged unique length 80,711,709 bp, so 50,723 bp (0.063%) are counted twice. Notebook 00's `tables/bed_territory.csv` is the authoritative version of these numbers.

The legacy TMB (`TMB` in `vcfcomparison_df_full.csv`, and `tmb_list_snp.csv`) equals variants per run divided by the raw sum in millions, for all 480 runs, to within 7.1e-15. The merged length does not reproduce it (largest difference 2.6e-2). Using the merged length would raise every TMB by 0.063%. Notebook 00 asserts both statements in its parity cell.

Nothing computes TMB in `src/swb/` or the notebooks yet except `swb.metrics.tmb(n_variants, region_size_bp)`, which takes the territory as an argument. Open question: raw sum (legacy) or merged length for new analysis.
