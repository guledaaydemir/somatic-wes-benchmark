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
