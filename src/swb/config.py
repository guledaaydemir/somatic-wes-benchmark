"""Paths and constants. Every path is relative to the repository root; nothing is read from outside it."""
from pathlib import Path

SEED = 0

# --- repository layout ----------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
RAW_SNV = DATA / "raw" / "snv"        # 480 per-run SNV VCFs   (gitignored)
RAW_INDEL = DATA / "raw" / "indel"    # 320 per-run indel VCFs (gitignored)
TRUTH = DATA / "truth"
REFERENCE = DATA / "reference"
DERIVED = DATA / "derived"
RESULTS = REPO_ROOT / "results"

# --- inputs ---------------------------------------------------------------
METADATA_CSV = REFERENCE / "TestCases.csv"
TRUTH_VCF = TRUTH / "hc_bed_filtered.recode.vcf"
SEQUENCING_YIELD_CSV = REFERENCE / "sequencing_yield.csv"    # manuscript Table 1: Gb and tumour/normal coverage per sample
EXOME_BED = REFERENCE / "sorted_exome_hc.bed.gz"
STRATIFICATION_DIR = REFERENCE / "stratification"    # *.bed.gz

SETS_PARQUET = DERIVED / "sets_dict.parquet"
MANIFEST = RESULTS / "00_data_preparation" / "tables" / "run_manifest.csv"

# --- legacy caches (copies of the legacy notebook's outputs, used for parity) --
LEGACY_METADATA = DERIVED / "vcfcomparison_df_full.csv"
LEGACY_FILENAMES = DERIVED / "vcf_filenames.csv"
LEGACY_FILENAMES_INDELS = DERIVED / "vcf_filenames_indels.csv"
LEGACY_SETS = DERIVED / "sets_dict.csv"
LEGACY_FILTERING = DERIVED / "filtering_df.csv"
LEGACY_UNION_METRICS = DERIVED / "union_metrics.csv"
LEGACY_TMB = DERIVED / "tmb_list_snp.csv"
LEGACY_TIME_TABLE = DERIVED / "df_timeTable.csv"
LEGACY_SETS_INDELS = DERIVED / "sets_dict_indels.csv"
LEGACY_META_INDELS = DERIVED / "vcfcomparison_full_indels.csv"
LEGACY_FILTERING_INDELS = DERIVED / "filtering_df_indels.csv"
# Values printed in the outputs of legacy/StabilityAnalysis_Indels.ipynb (the pair tables themselves are not kept)
LEGACY_INDEL_IOU_RECORDED = DERIVED / "legacy_indel_iou_recorded.csv"
# Values printed in the outputs of legacy/StabilityAnalysis.ipynb: F1 by sample (cell 126), factor tests (128), F1 of the first runs (132), variance table (134)
LEGACY_METRICS_RECORDED = DERIVED / "legacy_metrics_recorded.csv"

# --- design of the 480 SNV runs: factor levels, one run per combination ------
FACTOR_LEVELS = {
    "Sample": ["EA", "FD", "IL", "LL", "NC"],
    "Environment": ["Altay + COSAP", "Uhem + COSAP"],
    "isTrimmed": ["NO", "YES"],
    "baseRecalibration": ["NO", "YES"],
    "Mapper": ["BOWTIE", "BWA"],
    "VariantCaller": ["Mutect", "SomaticSniper", "Strelka"],
    "Duplicates": ["DELETE", "MARK"],
}

# --- design of the 320 indel runs (TestCases_Indels.csv spells the levels differently from TestCases.csv) ---
METADATA_INDEL_CSV = REFERENCE / "TestCases_Indels.csv"
FACTOR_LEVELS_INDEL = {
    "Sample": ["EA", "FD", "IL", "LL", "NC"],
    "Environment": ["Altay + COSAP", "Uhem + COSAP"],
    "isTrimmed": ["no", "yes"],
    "baseRecalibration": ["no", "yes"],
    "Mapper": ["BWA", "Bowtie"],
    "VariantCaller": ["Mutect", "Strelka"],
    "Duplicates": ["DELETE", "MARK"],
}

# --- paper-validated call sets: 3 aligners x 3 callers per centre (the folder NV is not one of the five centres) ---
VALIDATED_VCFS = REFERENCE / "validated_vcfs"
VALIDATED_CENTERS = ["LL", "NC", "EA", "IL", "FD"]  # the order of the legacy notebook
LEGACY_VALIDATED_FILES = {c: DERIVED / "validated_vcfs_{}_files.csv".format(c) for c in VALIDATED_CENTERS}
LEGACY_VALIDATED_SETS = {c: DERIVED / "validated_vcfs_{}_sets.csv".format(c) for c in VALIDATED_CENTERS}

# A pipeline: the settings that matter (environment and duplicate handling are dropped, notebooks 08 and 10)
PIPELINE_FACTORS = ["Sample", "Mapper", "VariantCaller", "isTrimmed", "baseRecalibration"]

# --- legacy pairwise IoU tables (rounded to 2 decimals), parity targets of notebook 02 ---
LEGACY_PAIRS = DERIVED / "df_generateds_r_p.csv"
LEGACY_PAIRS_BY_SAMPLE = {s: DERIVED / "df_{}.csv".format(s.lower()) for s in FACTOR_LEVELS["Sample"]}
LEGACY_GROUPED_IOU = DERIVED / "grouped_inter_union.csv"

# Parameter name (as in the legacy grouped IoU table) -> the one setting that differs within a pair of runs
IOU_CONTRASTS = {
    "Computing Environment": "Environment", "Mapper": "Mapper", "Variant Caller": "VariantCaller",
    "Base Recalibration": "baseRecalibration", "Mark Duplicates": "Duplicates", "Trimming": "isTrimmed",
}

# --- TestCases.csv schema (24 columns; the indel file has a different one) --
METADATA_COLUMNS = [
    "TestCaseNo", "phase", "Sample", "operatingSystem", "HOW", "isTrimmed", "baseRecalibration", "Mapper",
    "M_Version", "Variant Caller", "VC_Version", "ReferenceGenom", "Duplicates", "Variant Filter", "FilterPass",
    "SOMATICfilter", "remove_indels", "Sort", "File Link", "Result VCF Count", "Kim yaptı :)", "Slurm_id",
    "Donanım", "ElapsedTime",
]
# Names the legacy notebook gave these columns (cell 40).
METADATA_RENAME = {
    "HOW": "Environment",
    "Variant Caller": "VariantCaller",
    "Variant Filter": "VariantFilter",
    "File Link": "FileLink",
    "Kim yaptı :)": "Author",
    "Donanım": "Node",
}
# Empty in TestCases.csv; the legacy notebook filled them with computed values (TMB, variant count).
METADATA_COMPUTED_COLUMNS = ["Sort", "Result VCF Count"]


def results_dir(name):
    """results/<name>/ for a notebook, e.g. results_dir("01_variant_sets")."""
    return RESULTS / name
