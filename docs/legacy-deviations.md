# Legacy deviations

Legacy behaviours that look like bugs. Each entry: what the legacy code does, the evidence, and the decision. `src/swb/` reproduces the legacy behaviour unless an entry says it was decided otherwise (D1, D4, D5, D6, D7); nothing else is changed without an explicit decision.
Cell numbers are 0-based indices over all cells of `legacy/StabilityAnalysis.ipynb`.

## D1. BED filter ignores chromosomes
Status: **decided.** Both filters are implemented. **The corrected filter changes the region-stratified results substantially; published region-stratified values require revision.**

`filter_vcf_with_bed_single_chrom` only restricts by chromosome when `chrom` is given. Every legacy call (cells 215, 226, 237) omits it, so each variant's POS is compared against the regions of *all* chromosomes: a chr2 variant at POS 15 passes a chr1-only region `10-20`.

Implementation:
- `swb.io.bed_filter_legacy` reproduces the chromosome-blind behaviour exactly. It is used only in parity cells, and in notebook 07 to fill the comparison table. It gives the same rows as the legacy two-pointer sweep (kept as the reference in `tests/test_bed_filter.py`; equal on 25 random cases and on real data) without walking every region, which took 36 s per call on the largest BED.
- `swb.io.bed_filter` matches a variant only against regions of its own chromosome. It is the default for all new analysis. Its region boundary is still the legacy one (D2).
- The truth counts of the manuscript text (38, 1,123, 43, 1,118, 294 and 867 variants in the six regions) are the corrected filter's: the chromosome-blind filter keeps 685, 1,161, 859, 1,161, 1,058 and 1,156 of the 1,161 (notebook 09's `tables/region_truth.csv`). Notebook 07 computes both and writes `results/07_stratification/tables/legacy_vs_corrected.csv` (per-region counts under each) and `stratified_metrics.csv` (per run and region, both filters).

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

Implementation: `swb.stats.fit_anova` defaults to `typ=2`; `typ=1` is kept for reproducing legacy. With one run per cell the model is saturated (residual df 0), and statsmodels cannot form Type II there (it divides 0 by 0), so `swb.stats.type2_sum_sq` computes Type II by nested-model residual sums of squares. Notebook 05 also fits Type III (`typ=3`, sum-to-zero contrasts) on five factors, for the figure of D13: it equals Type II to 9e-14 on this balanced design. Notebook 05 reproduces the ten sources above 1% that legacy cell 134 printed (VariantCaller 53.94%, Center 19.39%, ...), in the same order, to 5e-7, the precision of the printed values. It matches statsmodels' `typ=2` to 2.6e-14 on unbalanced data with every cell populated (where Type I differs by up to 2.1), and `fit_anova` raises if a design cell is empty. The legacy formula equals `full_factorial_formula("f1_score")` term for term, in the same order (127 terms, checked).

## D5. Tumour sample: two legacy functions disagreed
Status: **decided**: the `##tumor_sample=` header is authoritative; the legacy functions are deleted, not ported.

`extract_tumor_sample_name` returned the first sample column and `_get_disease_sample` the last (cell 269). Neither is right, and nothing may infer the tumour from the file name, whose accession order lists the normal first.

`swb.io.tumor_sample` is the one resolver. It reads `##tumor_sample=` (which must name a sample column); else returns `TUMOR` when the sample columns are exactly `NORMAL` and `TUMOR`; else raises `ValueError`. Nothing in `src/swb/` infers the tumour from the file name or column order.

Evidence, headers of the 476 run VCFs held in `data/raw/vcf_snps` (4 of the 480 are not in the repository copy and were not checked): 156 Mutect VCFs resolve from the header, 320 Strelka and SomaticSniper VCFs from their `NORMAL`/`TUMOR` columns, none raise. In all 156 Mutect files the file name lists the normal first and the tumour second, matching the header. Mutect column order is not fixed: tumour first in 124 files, normal first in 32, so column position is unreliable even within one caller.

## D6. "Mapper" t-test line prints baseRecalibration numbers
Status: **decided**: legacy ported as `_legacy`, corrected version written for new analysis.

Cell 128 prints `Mapper :: T-Test: t-statistic = {t_stat_br}, p-value = {p_value_ttest_br}`, which are the baseRecalibration values; the Mapper values are `t_stat_ma` / `p_value_ttest_ma`. **Any Mapper t-statistic or p-value read from that printed line is a baseRecalibration comparison and must be re-run** with `swb.stats.factor_tests`; this includes any such value in the manuscript. The repository does not show which of the two the manuscript quotes: the p-value annotated under "Mapper" in the box-plot of cell 129 comes from `p_value_ttest_ma` and is the Mapper one, so a Mapper value taken from the plot is not affected. Check each manuscript Mapper value against its source.

Implementation: `swb.stats.factor_tests_legacy` reproduces the printed text as-is, mislabel included; `swb.stats.factor_tests` returns every comparison under its own label and is the one to use. Rows are still paired by position within each group, as in legacy. Notebook 05 reproduces every printed precision test with the runs in ascending `TestCaseNo` order, which is the legacy row order. The corrected mapper test for precision has t = -13.72 (p = 5.6e-32); the printed "Mapper" line, t = 13.60 (p = 1.4e-31), is the base recalibration test.

**Pairing (open).** `ttest_rel` pairs the k-th run of one level with the k-th run of the other. For the mapper (240 pairs) and for each caller comparison (160 pairs) every pair differs in the tested setting alone. For base recalibration and for trimming only 72 of 240 pairs do; in the other 168 the two runs also differ in trimming, or in base recalibration, respectively. So the legacy base recalibration and trimming t-tests compare unlike runs in 70% of their pairs. Notebook 05 keeps the legacy pairing so that the printed values reproduce, and writes the counts to `tables/pairing.csv`. Open question: keep it, or pair runs that differ in the tested setting alone.

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

Legacy scores each region-filtered call set against the *unfiltered* truth set (`get_precision/get_recall/get_f1score(filtered_variant_set, high_confidence)`, cells 215, 226, 237). Recall inside a region is therefore TP-in-region divided by all truth variants, and falls with the size of the region whatever the caller does. Confirmed by reproduction: computing the metrics against the unfiltered truth set (1,161 variants) reproduces all 2880 rows of `union_metrics.csv` to within 2.3e-16. Notebook 07 keeps this definition for both filters and states it in its purpose. Notebook 09 filters the truth set with the same BED as the calls, as the manuscript text describes, so recall is TP over the truth variants in the region; with it the truth counts, shares and medians of the manuscript text reproduce. Which definition to publish is still your decision.

## D10. Support fraction to support count
Status: **decided.** Numbered D10 because D2 is the BED boundary convention (pending); this entry has no earlier number.

`swb.metrics.fraction_to_count(q, n_lists)` transfers a consensus support fraction `q` to `n_lists` call sets, for `swb.metrics.ensemble_call`: `k = ceil(q * n_lists - 1e-9)`, clipped to `[1, n_lists]`. It raises `ValueError` if `q` is not in `(0, 1]` or `n_lists < 1`. Ceil keeps the transferred rule at least as strict as the rule it was chosen under. The epsilon absorbs float error on exact ratios: `7/25 * 25` is `7.000000000000001`, and a bare ceil gives 8 instead of 7. (`34/120 * 120` is exactly `34.0` in floating point and needs no guard, but other exact ratios do.)

Values: `(34/120, 120) -> 34`, `(34/120, 96) -> 28`, `(34/120, 24) -> 7`, `(0.5, 24) -> 12`, `(1/96, 24) -> 1`, `(1.0, 24) -> 24`. `k` is non-decreasing in `q` and in `n_lists`.

Tests in `tests/test_metrics.py`: `test_fraction_to_count_known_values`, `test_fraction_to_count_non_decreasing_in_q_and_in_n_lists`, `test_fraction_to_count_equals_exact_integer_ceil_on_rational_grid`, and the rejection tests. No legacy cell computing this conversion has been traced, so there is no legacy behaviour to reproduce.

## D11. TMB territory counts overlapping BED intervals twice
Status: **pending decision** (metric definitions are fixed; the legacy value is reproduced).

`calculate_region_size` (cell 38) sums `end - start` over every line of the BED, so a base covered by two intervals counts twice. `sorted_exome_hc.bed.gz` has 272,014 intervals: raw sum 80,762,432 bp, merged unique length 80,711,709 bp, so 50,723 bp (0.063%) are counted twice. Notebook 00's `tables/bed_territory.csv` is the authoritative version of these numbers.

The legacy TMB (`TMB` in `vcfcomparison_df_full.csv`, and `tmb_list_snp.csv`) equals variants per run divided by the raw sum in millions, for all 480 runs, to within 7.1e-15. The merged length does not reproduce it (largest difference 2.6e-2). Using the merged length would raise every TMB by 0.063%. Notebook 00 asserts both statements in its parity cell.

Notebook 04 (`04_tmb`) computes TMB with both territories: `tmb` (raw sum, the legacy value; its parity target) and `tmb_merged`, in `tables/tmb.csv`. It asserts that the merged territory does not reproduce the legacy TMB. `swb.metrics.tmb(n_variants, region_size_bp)` takes the territory as an argument. Open question: raw sum (legacy) or merged length for new analysis.

## D12. Pairwise IoU is rounded to 2 decimals
Status: **pending decision** (metric definitions are fixed; the legacy values are reproduced).

The legacy pair tables use `inter_union` (cell 8), which returns `round(intersection / union, 2)`: `df_generateds_r_p.csv` (all ordered pairs of the 480 runs, self-pairs included), the per-sample `df_ea.csv` to `df_nc.csv` (cell 68), and `grouped_inter_union.csv` (cell 101), whose means are taken over the rounded values, as is the bar chart drawn from it. The IoU in `union_metrics.csv` (against the truth set) uses `intersection_over_union` and is not rounded.

Effect, from the cached sets: a rounded pair differs from the exact value by at most 0.005, and the 30 grouped means (6 parameters x 5 samples) by at most 0.0014. Rounding also hides small differences: the computing-environment contrast for sample IL is 1.00 rounded and 0.999986 unrounded. Notebook 02's `tables/grouped_iou.csv` (`mean_iou` against `mean_iou_2dp`) and its audit summary are the authoritative version of these numbers.

For the 320 indel runs (notebook `02_pairwise_iou_indels`, from `sets_dict_indels.csv`) the same rounding changes a pair by at most 0.005 and the grouped means by at most 0.0015. Notebook 02 writes both values. Its parity cell and figure use the rounded one, which reproduces the legacy tables to 1e-12; `iou` is unrounded and is not compared with anything. Open question: rounded (legacy) or unrounded for new analysis.

## D13. The ANOVA p-values rest on near-identical replicates
Status: **pending decision** (the legacy variance percentages are unaffected).

The legacy 7-factor ANOVA has one run per design cell, so it is saturated: it has no residual, hence no F statistics and no p-values (D4). The Type III ANOVA figure of notebook 05 needs p-values. The only full-factorial model with residual degrees of freedom treats environment and duplicate handling as replicates: five factors (mapper, caller, center, trimming, base recalibration), 120 cells of 4 runs, 360 residual degrees of freedom. The design is balanced, so each term's share of the total sum of squares is the same as in the 7-factor table (the ten sources legacy printed agree to 5e-7), and Type III equals Type II.

The four runs of a cell differ only in environment and duplicate handling, and they are practically identical (D12: their IoU is 1.00 rounded). The residual is therefore 3.4e-9 of the total sum of squares. Every F statistic is between 9.6e5 and 2.9e10, and every term has a p-value far below 1e-300 (scipy's p-value underflows to 0; -log10 p, computed in log space, is 646 to 1477). The smallest term (trimming and base recalibration, 0.0013% of the total) is as significant as the largest (variant caller, 53.9%), so these p-values do not say which terms matter. Notebook 05's `tables/anova_f1_type3.csv` is the authoritative version of these numbers.

The figure `anova_contribution_significance` is drawn as specified, with this caveat in its footnote. Open question: whether the ANOVA should carry p-values at all, or use another error term or a permutation test.

## D14. A printed count carries the wrong center's label
Status: **decided**: reproduced as it is and documented; nothing depends on the label.

Legacy cell 23 ends with `print_set_count("EA", set.intersection(*ll_sets_dict.values()))` and prints `EA Validated Set Count: 2388`. The sets are the LL ones: 2388 is the size of the LL intersection, and EA's is 2464. Notebook 06 asserts both. Cells 28 and 29 build the Venn from the correct set of each center.

## D15. The validated sets keep their indels
Status: **pending decision**

The section of the legacy notebook that reads the validated VCFs is titled "(without indels)", but `process_files` (cell 23) removes neither indels nor non-PASS records: it keeps the first ALT of every record. Of the 575,670 distinct variants in the 45 validated VCFs, 26,778 are indels (REF or ALT_1 not a single base). Notebook 06 counts them per file in `tables/validated_sets.csv` (`n_snv`, `n_indel`). The truth set has no indels, and the intersection of each center's nine sets, which the Venn uses, contains none, so the figure is not affected. `swb.io.validated_sets` keeps the legacy behaviour. Open question: whether the validated sets should be reduced to SNVs.

## D16. The manuscript's consensus optimum and gain do not match the legacy algorithm
Status: **pending decision** (the legacy algorithm is reproduced as it is).

The legacy consensus cells (203 to 211) were never executed in the legacy notebook, so no output was recorded. Notebook 08 runs them: 120 pipelines (the first run of each group of four runs that differ in environment and duplicate handling), a support count per variant, and the precision, recall and F1 of the variants called by at least n pipelines against the truth set, for n = 1 to 120. The curve equals `ensemble_call` scored with `precision`, `recall` and `f1` at every n.

Against the manuscript text (F1 0.94 with variants called by at least 34 pipelines, +1.97% over the best single pipeline):
- The maximum F1 is 0.94092 at **n = 33**, not 34. n = 34 gives 0.94086, 5.2e-5 lower; both round to 0.94. The maximum is unique. F1 stays within 1% of it from n = 26 to 47.
- The gain over the best single pipeline (run 776, F1 0.92182) is **+2.07%** (0.01910 absolute) from the unrounded F1. The manuscript's +1.97% is what the rounded 0.94 gives against 0.92182.
- Which of the four runs of a group represents its pipeline changes the F1 curve by at most 1.8e-5, and not the optimum.
- Another paragraph of the manuscript text says the maximum is reached with "approximately 30" pipelines. The plateau contains 30, 33 and 34.

Notebook 08's `tables/consensus_curve.csv` and `tables/consensus_best.csv` are the authoritative numbers. Open question: which of n = 33 and n = 34 the manuscript should report, and whether the gain should be quoted from the unrounded F1.

## D17. The ensemble analysis has no legacy implementation
Status: **decided**: reproduced from the manuscript text; documented, nothing to fix.

The legacy notebook has no code for the exhaustive search of pipeline combinations or the leave-one-center-out (LOCO) validation (its cells with `combinations` belong to a vcfeval search, notebook cells 242 to 252). Notebook 10 implements the method as the manuscript describes it, scored with the ported `ensemble_call`, `precision`, `recall` and `f1`; parity checks the vectorised search against those functions exactly and every number in the manuscript text against the results. All of them reproduce. Choices the text leaves open, fixed as follows:
- Ties in the maximum: the first one wins (smallest n, then lowest combination index, then lowest threshold). No tie occurs in the reported numbers.
- The SD of the five held-out F1 values is the sample SD (0.0431; 0.0386 with n in the denominator); only the sample SD gives the manuscript's 0.043.
- The gain in LOCO is the held-out F1 of the selected ensemble minus the held-out F1 of the individual pipeline with the highest mean training F1 (mean +0.0176; the text's 0.018).
- The manuscript's 0.0045 gap of the Illumina fold to the reported configuration is the difference of the 4-decimal values 0.9374 and 0.9329; unrounded it is 0.00445 (0.93736 and 0.93291).
- The individual pipelines are the first run of each group of four, as in notebook 08 (D16); the choice is immaterial there (at most 1.8e-5 in F1).

