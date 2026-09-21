"""Paths and constants. Reads SWB_DATA_ROOT (default: /Volumes/E4 Pro/Bioinformatics-StabilityAnalysis).

DATA_ROOT holds the raw inputs (vcf/...); everything else lives inside the repository.
"""
import os
from pathlib import Path

SEED = 0

# --- external drive -------------------------------------------------------
DATA_ROOT = Path(os.environ.get("SWB_DATA_ROOT", "/Volumes/E4 Pro/Bioinformatics-StabilityAnalysis"))
VCF_DIR = DATA_ROOT / "vcf" / "TESTCASES_bedded"
METADATA_CSV = DATA_ROOT / "vcf" / "TestCases.csv"
TRUTH_VCF = DATA_ROOT / "vcf" / "hc_bed_filtered.recode.vcf"
EXOME_BED = DATA_ROOT / "vcf" / "sorted_exome_hc.bed"
STRATIFICATION_DIR = DATA_ROOT / "vcf" / "genome_stratifications"

# --- repository -----------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DERIVED = REPO_ROOT / "data" / "derived"
RESULTS = REPO_ROOT / "results"

SETS_PARQUET = DERIVED / "sets_dict.parquet"
MANIFEST = RESULTS / "00_data_preparation" / "tables" / "run_manifest.csv"

# --- legacy caches (copies of the legacy notebook's outputs, used for parity) --
LEGACY_METADATA = DERIVED / "vcfcomparison_df_full.csv"
LEGACY_FILENAMES = DERIVED / "vcf_filenames.csv"
LEGACY_SETS = DERIVED / "sets_dict.csv"
LEGACY_FILTERING = DERIVED / "filtering_df.csv"
LEGACY_UNION_METRICS = DERIVED / "union_metrics.csv"

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
