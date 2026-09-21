"""bed_filter_legacy (chromosome-blind, for parity) and bed_filter (chromosome-aware, the default).

docs/legacy-deviations.md D1. Both use a start <= POS < end region (D2, still open).
"""
import numpy as np
import pandas as pd
import pytest

from swb import io


def bed(rows):
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def calls(pos, chroms=None):
    chroms = chroms or ["chr1"] * len(pos)
    return pd.DataFrame({"#CHROM": chroms, "POS": pos, "REF": "A", "ALT": "T"})


def legacy_sweep(vcf_df, bed_df, invert_stringency, chrom=None,
                 vcf_chrom_col_name="#CHROM", bed_chrom_col_name="chrom"):
    """The legacy two-pointer sweep (cell 226), verbatim: the reference bed_filter_legacy must agree with."""
    if chrom:
        temp_vcf_df = vcf_df[vcf_df[vcf_chrom_col_name] == chrom]
        temp_bed_df = bed_df[bed_df[bed_chrom_col_name] == chrom]
    else:
        temp_vcf_df = vcf_df.copy()
        temp_bed_df = bed_df.copy()

    mask = np.zeros(len(temp_vcf_df), dtype=bool)
    sorted_vcf_df = temp_vcf_df.sort_values(by=["POS"]).reset_index(drop=True)
    sorted_bed_df = temp_bed_df.sort_values(by=["start", "end"]).reset_index(drop=True)

    vi, ri = 0, 0
    while vi < len(sorted_vcf_df) and ri < len(sorted_bed_df):
        pos = sorted_vcf_df.at[vi, "POS"]
        start = sorted_bed_df.at[ri, "start"]
        end = sorted_bed_df.at[ri, "end"]
        if pos < start:
            vi += 1
        elif start <= pos < end:
            mask[vi] = True
            vi += 1
        else:
            ri += 1

    final_mask = ~mask if invert_stringency else mask
    return sorted_vcf_df[final_mask]


def rows(df):
    return sorted(map(tuple, df[["#CHROM", "POS", "REF", "ALT"]].values.tolist()))


# --- bed_filter_legacy -------------------------------------------------------

def test_legacy_region_is_start_inclusive_end_exclusive():
    kept = io.bed_filter_legacy(calls([9, 10, 15, 19, 20]), bed([("chr1", 10, 20)]), invert_stringency=False)
    assert list(kept["POS"]) == [10, 15, 19]


def test_legacy_invert_and_sorted_output():
    out = io.bed_filter_legacy(calls([50, 12, 5]), bed([("chr1", 10, 20)]), invert_stringency=True)
    assert list(out["POS"]) == [5, 50]
    assert list(out.index) == [0, 2]  # positions in the sorted frame: 5 -> 0, 12 -> 1, 50 -> 2


def test_legacy_chrom_restricts_when_given():
    vcf = calls([15, 15], ["chr1", "chr2"])
    out = io.bed_filter_legacy(vcf, bed([("chr1", 10, 20)]), invert_stringency=False, chrom="chr1")
    assert list(out["#CHROM"]) == ["chr1"]


def test_legacy_ignores_chromosome_when_chrom_is_none():
    # Characterises legacy behaviour (D1): the chr2 call is kept by a chr1-only region.
    vcf = calls([15, 15], ["chr1", "chr2"])
    out = io.bed_filter_legacy(vcf, bed([("chr1", 10, 20)]), invert_stringency=False)
    assert sorted(out["#CHROM"]) == ["chr1", "chr2"]


def test_legacy_overlapping_and_nested_regions():
    regions = bed([("chr1", 0, 100), ("chr1", 10, 20), ("chr1", 150, 160)])
    out = io.bed_filter_legacy(calls([5, 15, 99, 100, 155, 200]), regions, invert_stringency=False)
    assert list(out["POS"]) == [5, 15, 99, 155]


@pytest.mark.parametrize("seed", range(25))
def test_legacy_agrees_with_verbatim_sweep_on_random_data(seed):
    rng = np.random.RandomState(seed)
    chroms = ["chr1", "chr2", "chr3"]
    n_regions = rng.randint(0, 40)
    starts = rng.randint(0, 200, n_regions)
    regions = bed(list(zip(rng.choice(chroms[:2], n_regions), starts, starts + rng.randint(1, 60, n_regions))))
    n_calls = rng.randint(0, 60)  # positions repeat and cross chromosomes on purpose
    vcf = pd.DataFrame({"#CHROM": rng.choice(chroms, n_calls), "POS": rng.randint(0, 260, n_calls),
                        "REF": "A", "ALT": "T"})
    for invert in (False, True):
        for chrom in (None, "chr1", "chr3"):
            got = io.bed_filter_legacy(vcf, regions, invert, chrom=chrom)
            want = legacy_sweep(vcf, regions, invert, chrom=chrom)
            assert rows(got) == rows(want), (seed, invert, chrom)


def test_prepared_regions_give_the_same_result_as_a_dataframe():
    regions = bed([("chr1", 10, 20), ("chr2", 0, 5)])
    vcf = calls([12, 3, 30], ["chr1", "chr2", "chr1"])
    prepared = io.BedRegions(regions)
    assert rows(io.bed_filter_legacy(vcf, prepared, False)) == rows(io.bed_filter_legacy(vcf, regions, False))
    assert rows(io.bed_filter(vcf, prepared)) == rows(io.bed_filter(vcf, regions))


# --- bed_filter (corrected) --------------------------------------------------

def test_corrected_respects_chromosome():
    # Same input as test_legacy_ignores_chromosome_when_chrom_is_none: only the chr1 call is inside a chr1 region.
    vcf = calls([15, 15], ["chr1", "chr2"])
    out = io.bed_filter(vcf, bed([("chr1", 10, 20)]), invert_stringency=False)
    assert list(out["#CHROM"]) == ["chr1"]
    inv = io.bed_filter(vcf, bed([("chr1", 10, 20)]), invert_stringency=True)
    assert list(inv["#CHROM"]) == ["chr2"]


def test_corrected_matches_one_legacy_call_per_chromosome():
    regions = bed([("chr1", 10, 20), ("chr2", 100, 200), ("chr2", 150, 300)])
    vcf = calls([5, 15, 19, 20, 120, 250, 400, 15], ["chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr3"])
    want = pd.concat([io.bed_filter_legacy(vcf, regions, False, chrom=c) for c in ["chr1", "chr2", "chr3"]])
    assert rows(io.bed_filter(vcf, regions)) == rows(want)
    assert list(io.bed_filter(vcf, regions)["POS"]) == [15, 19, 120, 250]


def test_corrected_keeps_the_legacy_boundary_convention():
    # D2 is undecided: only the chromosome handling differs from legacy.
    kept = io.bed_filter(calls([9, 10, 15, 19, 20]), bed([("chr1", 10, 20)]))
    assert list(kept["POS"]) == [10, 15, 19]


def test_corrected_agrees_with_legacy_on_a_single_chromosome():
    regions = bed([("chr1", 0, 100), ("chr1", 10, 20), ("chr1", 150, 160)])
    vcf = calls([5, 15, 99, 100, 155, 200])
    for invert in (False, True):
        assert rows(io.bed_filter(vcf, regions, invert)) == rows(io.bed_filter_legacy(vcf, regions, invert))


def test_corrected_chromosome_missing_from_bed_has_no_regions():
    vcf = calls([15, 15], ["chr1", "chrM"])
    regions = bed([("chr1", 10, 20)])
    assert list(io.bed_filter(vcf, regions)["#CHROM"]) == ["chr1"]
    assert list(io.bed_filter(vcf, regions, invert_stringency=True)["#CHROM"]) == ["chrM"]


def test_corrected_raises_when_no_chromosome_name_is_shared():
    with pytest.raises(ValueError, match="no chromosome name"):
        io.bed_filter(calls([15], ["chr1"]), bed([("1", 10, 20)]))


def test_corrected_rejects_missing_chromosome_names():
    with pytest.raises(ValueError, match="missing chromosome"):
        io.bed_filter(calls([15], [None]), bed([("chr1", 10, 20)]))


def test_corrected_empty_input():
    out = io.bed_filter(calls([]), bed([("chr1", 10, 20)]))
    assert len(out) == 0 and list(out.columns) == ["#CHROM", "POS", "REF", "ALT"]


def test_corrected_output_is_sorted_with_fresh_index():
    vcf = calls([50, 12, 5, 15], ["chr2", "chr1", "chr1", "chr2"])
    out = io.bed_filter(vcf, bed([("chr1", 0, 100), ("chr2", 0, 100)]))
    assert list(zip(out["#CHROM"], out["POS"])) == [("chr1", 5), ("chr1", 12), ("chr2", 15), ("chr2", 50)]
    assert list(out.index) == [0, 1, 2, 3]
