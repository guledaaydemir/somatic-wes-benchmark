"""File discovery, VCF reading and filtering, variant keys, BED helpers and set caches.

Ported from legacy/StabilityAnalysis.ipynb; each docstring names the legacy function it replaces.
Variant keys are CHROM_POS_REF_ALT1 strings, e.g. "chr1_12345_A_T".
"""
import os

import numpy as np
import pandas as pd

SOMATICSNIPER = "SomaticSniper"
KEY_COLUMNS = ["CHROM", "POS", "REF", "ALT_1"]
EXCLUDED_NAMES = frozenset({".DS_Store", "fixed.recode.vcf", "gatk.recode.vcf"})


# --- file discovery -------------------------------------------------------

def real_files(root, suffix=".vcf", exclude=EXCLUDED_NAMES):
    """Sorted paths under root ending in suffix. Skips macOS AppleDouble '._' files and excluded names."""
    found = []
    for dirpath, _, names in os.walk(root):
        for fn in names:
            if fn.startswith("._") or fn in exclude or not fn.endswith(suffix):
                continue
            found.append(os.path.join(dirpath, fn))
    return sorted(found)


def canon(x):
    """Normalise a TestCaseNo ("TestCase 269", "269", 269, 269.0) to "269"."""
    if isinstance(x, float) and x.is_integer():
        x = int(x)
    s = str(x).strip()
    return s[len("TestCase "):].strip() if s.lower().startswith("testcase ") else s


def vcf_filenames(base):
    """{canon(folder name): VCF path} (legacy: vcf_filenames, cell 15; last file per folder wins)."""
    return {canon(os.path.basename(os.path.dirname(p))): p for p in real_files(base)}


def validated_vcf_files(prefix):
    """{folder relative to prefix: path} for validated-set files (legacy: cell 23 os.walk)."""
    return {
        os.path.relpath(os.path.dirname(p), prefix): p
        for p in real_files(prefix, suffix="", exclude={".DS_Store"})
    }


def bed_files(directory):
    """Sorted '.bed' paths directly in directory, skipping '._' files (legacy: get_bed_file_names)."""
    return sorted(
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".bed") and not f.startswith("._")
    )


# --- VCF reading ----------------------------------------------------------

def read_variants(path):
    """scikit-allel dataframe of variants/* (CHROM, POS, REF, ALT_1, FILTER_PASS, is_snp, ...); the legacy reader."""
    import allel

    return allel.vcf_to_dataframe(path, fields="variants/*")


def vcf_header(path):
    """Header lines ('##...' and '#CHROM'). Decodes with errors='replace': some VCF headers hold non-UTF-8 bytes."""
    lines = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("#"):
                break
            lines.append(line.rstrip("\r\n"))
    return lines


def tumor_sample(path):
    """Name of the tumour sample of a VCF: the one resolver, never inferred from file name or column order.

    1. the '##tumor_sample=' header (Mutect VCFs), which must name one of the sample columns;
    2. else the column literally named TUMOR, when the sample columns are exactly NORMAL and TUMOR
       (Strelka and SomaticSniper VCFs);
    3. else ValueError.
    File names list the normal sample first and the column order differs by caller, so neither says
    which sample is the tumour (docs/legacy-deviations.md D5).
    """
    header = vcf_header(path)
    samples = []
    for line in header:
        if line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
    for line in header:
        if line.startswith("##tumor_sample="):
            name = line.split("=", 1)[1].strip()
            if name not in samples:
                raise ValueError("{}: ##tumor_sample={!r} is not a sample column {}".format(path, name, samples))
            return name
    if sorted(samples) == ["NORMAL", "TUMOR"]:
        return "TUMOR"
    raise ValueError("{}: no ##tumor_sample= header and sample columns {} are not NORMAL/TUMOR".format(path, samples))


# --- filters and keys -----------------------------------------------------

def apply_pass_filter(df, caller):
    """Keep FILTER_PASS rows. SomaticSniper (exact name) is exempt: its FILTER is '.' on every record."""
    df = df.copy()
    return df if caller == SOMATICSNIPER else df[df["FILTER_PASS"].astype("bool")]


def filters_control(df, caller):
    """PASS filter (SomaticSniper exempt), then is_snp (legacy: filtersControl)."""
    df = apply_pass_filter(df, caller)
    return df[df["is_snp"].astype("bool")]


def remove_indels(df):
    """Keep rows whose REF and ALT_1 are single bases."""
    df = df[df["REF"].apply(len) == 1]
    return df[df["ALT_1"].apply(len) == 1]


def variant_set(df):
    """{CHROM_POS_REF_ALT1} for every row."""
    if df.empty:
        return set()
    key = df["CHROM"].astype(str)
    for col in KEY_COLUMNS[1:]:
        key = key + "_" + df[col].astype(str)
    return set(key)


def parse_vcf(path, caller, indels=False):
    """(variant set, records read, records kept).

    SNP VCFs: PASS (SomaticSniper exempt) + is_snp + 1-bp REF/ALT (legacy cell 32).
    Indel VCFs: PASS only (legacy cell 34).
    """
    df = read_variants(path)
    total = len(df)
    df = apply_pass_filter(df, caller) if indels else remove_indels(filters_control(df, caller))
    return variant_set(df), total, len(df)


def truth_set(path):
    """High-confidence truth set: no FILTER test, indels removed (legacy: high_confidence, cell 26)."""
    return variant_set(remove_indels(read_variants(path)))


def validated_sets(files):
    """{key: variant set} for validated VCFs; no FILTER test, no indel removal, empty VCFs dropped (legacy: process_files)."""
    sets = {}
    for key in files:
        df = read_variants(files[key])
        if len(df):
            sets[key] = variant_set(df)
    return sets


# --- caches and tables ----------------------------------------------------

def _makedirs_for(path):
    os.makedirs(os.path.dirname(os.fspath(path)) or ".", exist_ok=True)


def save_sets(sets, path):
    """Write {key: set} as (Key, Identifier) rows sorted by both; .parquet or .csv by suffix.

    Keys with an empty set are not written (same as the legacy CSV cache).
    """
    rows = sorted((str(k), v) for k, s in sets.items() for v in s)
    df = pd.DataFrame(rows, columns=["Key", "Identifier"])
    _makedirs_for(path)
    if os.fspath(path).endswith(".parquet"):
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, index=False)


def load_sets(path):
    """{str(key): set} from a save_sets cache, or the legacy sets_dict.csv (Key, Identifier)."""
    read = pd.read_parquet if os.fspath(path).endswith(".parquet") else pd.read_csv
    sets = read(path).groupby("Key")["Identifier"].apply(set).to_dict()
    return {str(k): v for k, v in sets.items()}


def write_csv(df, path, sort_by):
    """Sorted, index-free CSV so reruns are byte-identical."""
    _makedirs_for(path)
    df.sort_values(sort_by, kind="mergesort").to_csv(path, index=False)


# --- BED regions ----------------------------------------------------------

def read_bed(path):
    """BED as chrom/start/end columns (legacy: pd.read_csv(..., sep='\\t', header=None, names=[...]))."""
    return pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end"])


def region_size(bed_file):
    """Sum of (end - start) over BED lines, overlaps counted twice (legacy: calculate_region_size)."""
    size = 0
    with open(bed_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            size += int(fields[2]) - int(fields[1])
    return size


def merge_intervals(intervals):
    """Merge overlapping (start, end) pairs."""
    merged = []
    for s, e in sorted(intervals):
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return merged


def unique_region_size(bed_df):
    """Bases covered per chromosome after merging overlaps (legacy: calculate_unique_region_size)."""
    chrom_intervals = {}
    for _, r in bed_df.iterrows():
        chrom_intervals.setdefault(r["chrom"], []).append((int(r["start"]), int(r["end"])))
    total = 0
    for ivs in chrom_intervals.values():
        for s, e in merge_intervals(ivs):
            total += e - s
    return total


def variants_to_df(variant_set_):
    """DataFrame (#CHROM, POS, REF, ALT) from CHROM_POS_REF_ALT1 keys (legacy: convert_to_vcf_df)."""
    df = pd.DataFrame([v.split("_") for v in variant_set_], columns=["#CHROM", "POS", "REF", "ALT"])
    df["POS"] = df["POS"].astype(int)
    return df


def df_to_variants(vcf_df):
    """CHROM_POS_REF_ALT1 keys from a (#CHROM, POS, REF, ALT) frame (legacy: convert_to_variant_set)."""
    return set(
        "{}_{}_{}_{}".format(c, p, r, a)
        for c, p, r, a in zip(vcf_df["#CHROM"], vcf_df["POS"], vcf_df["REF"], vcf_df["ALT"])
    )


class BedRegions:
    """BED intervals prepared once for repeated membership tests.

    A region is start <= POS < end (legacy convention, docs/legacy-deviations.md D2). Regions are kept
    pooled over all chromosomes (for bed_filter_legacy) and per chromosome (for bed_filter).
    """

    def __init__(self, bed_df, bed_chrom_col_name="chrom"):
        self.pooled = _prepare_intervals(bed_df["start"], bed_df["end"])
        self.by_chrom = {
            c: _prepare_intervals(g["start"], g["end"]) for c, g in bed_df.groupby(bed_chrom_col_name)
        }


_NO_INTERVALS = (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))


def _prepare_intervals(starts, ends):
    """Starts sorted, with the running maximum of the ends in that order."""
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    order = np.lexsort((ends, starts))
    return starts[order], np.maximum.accumulate(ends[order])


def _inside(prepared, pos):
    """True where some interval has start <= pos < end."""
    starts, running_max_end = prepared
    pos = np.asarray(pos, dtype=np.int64)
    if len(starts) == 0:
        return np.zeros(len(pos), dtype=bool)
    last = np.searchsorted(starts, pos, side="right") - 1  # last interval starting at or before pos
    return (last >= 0) & (running_max_end[np.maximum(last, 0)] > pos)


def bed_filter_legacy(vcf_df, bed, invert_stringency, chrom=None,
                      vcf_chrom_col_name="#CHROM", bed_chrom_col_name="chrom"):
    """Legacy BED filter, reproduced exactly (legacy: filter_vcf_with_bed_single_chrom). For parity cells only.

    With chrom=None (how every legacy call used it) POS is compared against the regions of ALL
    chromosomes, so a chr2 variant is kept by a chr1 region; docs/legacy-deviations.md D1. Use bed_filter
    for analysis.

    Keeps rows of vcf_df inside (or, if invert_stringency, outside) the regions. bed is a DataFrame
    (chrom, start, end) or a BedRegions. Result is sorted by POS; index labels are row positions in the
    sorted frame (reset before masking). Same result as the legacy two-pointer sweep, which
    tests/test_bed_filter.py keeps as the reference; this version needs no walk over every region.
    """
    regions = bed if isinstance(bed, BedRegions) else BedRegions(bed, bed_chrom_col_name)
    if chrom:
        temp_vcf_df = vcf_df[vcf_df[vcf_chrom_col_name] == chrom]
        prepared = regions.by_chrom.get(chrom, _NO_INTERVALS)
    else:
        temp_vcf_df = vcf_df.copy()
        prepared = regions.pooled

    sorted_vcf_df = temp_vcf_df.sort_values(by=["POS"], kind="mergesort").reset_index(drop=True)
    mask = _inside(prepared, sorted_vcf_df["POS"])
    final_mask = ~mask if invert_stringency else mask
    return sorted_vcf_df[final_mask]


def bed_filter(vcf_df, bed, invert_stringency=False,
               vcf_chrom_col_name="#CHROM", bed_chrom_col_name="chrom"):
    """Chromosome-aware BED filter: the default for all analysis (corrects docs/legacy-deviations.md D1).

    Keeps rows of vcf_df whose POS lies in a region of the SAME chromosome (or, if invert_stringency,
    in none). A chromosome absent from the BED has no regions. Boundaries are unchanged from legacy
    (start <= POS < end, D2 is still open). bed is a DataFrame (chrom, start, end) or a BedRegions.
    Result is sorted by (chrom, POS) with a fresh index. Raises ValueError if the VCF and the BED share
    no chromosome name (e.g. 'chr1' vs '1'), instead of silently dropping everything.
    """
    regions = bed if isinstance(bed, BedRegions) else BedRegions(bed, bed_chrom_col_name)
    chroms = vcf_df[vcf_chrom_col_name]
    if chroms.isna().any() or (chroms == "").any():
        raise ValueError("variants with a missing chromosome name")
    if len(vcf_df) and not set(chroms) & set(regions.by_chrom):
        raise ValueError("VCF and BED share no chromosome name: VCF {}, BED {}".format(
            sorted(set(chroms))[:3], sorted(regions.by_chrom)[:3]))

    ordered = vcf_df.sort_values(by=[vcf_chrom_col_name, "POS"], kind="mergesort").reset_index(drop=True)
    pos = ordered["POS"].to_numpy()
    mask = np.zeros(len(ordered), dtype=bool)
    for chrom, rows in ordered.groupby(vcf_chrom_col_name).indices.items():
        mask[rows] = _inside(regions.by_chrom.get(chrom, _NO_INTERVALS), pos[rows])
    final_mask = ~mask if invert_stringency else mask
    return ordered[final_mask].reset_index(drop=True)
