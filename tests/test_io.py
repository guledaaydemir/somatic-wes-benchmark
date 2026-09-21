"""Toy-data checks for swb.io: hazards from CLAUDE.md, filters, keys, caches, BED helpers."""
import pandas as pd
import pytest

from swb import io


def test_canon():
    assert io.canon("TestCase 269") == "269"
    assert io.canon("testcase 269") == "269"
    assert io.canon("269") == "269"
    assert io.canon(269) == "269"
    assert io.canon(269.0) == "269"
    assert io.canon(" TestCase 5 ") == "5"


def test_real_files_skips_appledouble_and_excluded(tmp_path):
    (tmp_path / "TestCase 1").mkdir()
    (tmp_path / "TestCase 2").mkdir()
    keep = [tmp_path / "TestCase 1" / "a.vcf", tmp_path / "TestCase 2" / "b.vcf"]
    for p in keep:
        p.write_text("x")
    for name in ["._a.vcf", ".DS_Store", "fixed.recode.vcf", "gatk.recode.vcf", "c.vcf.gz", "d.txt"]:
        (tmp_path / "TestCase 1" / name).write_text("x")
    assert io.real_files(str(tmp_path)) == sorted(str(p) for p in keep)


def test_vcf_filenames_keys_are_canon(tmp_path):
    (tmp_path / "TestCase 269").mkdir()
    (tmp_path / "TestCase 269" / "x.vcf").write_text("x")
    (tmp_path / "TestCase 269" / "._x.vcf").write_text("x")
    assert list(io.vcf_filenames(str(tmp_path))) == ["269"]


def test_validated_vcf_files_keys_relative_folder(tmp_path):
    (tmp_path / "bowtie" / "mutect").mkdir(parents=True)
    (tmp_path / "bowtie" / "mutect" / "WES.vcf").write_text("x")
    (tmp_path / "bowtie" / "mutect" / ".DS_Store").write_text("x")
    assert list(io.validated_vcf_files(str(tmp_path))) == ["bowtie/mutect"]


def test_bed_files(tmp_path):
    for name in ["b.bed", "a.bed", "._a.bed", "c.txt"]:
        (tmp_path / name).write_text("x")
    assert [p.split("/")[-1] for p in io.bed_files(str(tmp_path))] == ["a.bed", "b.bed"]


def variants(**cols):
    base = {
        "CHROM": ["chr1", "chr1", "chr2", "chr2"],
        "POS": [100, 200, 300, 400],
        "REF": ["A", "C", "G", "AT"],
        "ALT_1": ["T", "G", "A", "A"],
        "FILTER_PASS": [True, False, True, True],
        "is_snp": [True, True, False, False],
    }
    base.update(cols)
    return pd.DataFrame(base)


def test_pass_filter_drops_non_pass_for_mutect_and_strelka():
    for caller in ["Mutect", "Strelka"]:
        assert list(io.apply_pass_filter(variants(), caller)["POS"]) == [100, 300, 400]


def test_pass_filter_exempts_somaticsniper():
    df = variants(FILTER_PASS=[False, False, False, False])  # FILTER is '.' on every SomaticSniper record
    assert len(io.apply_pass_filter(df, "SomaticSniper")) == 4
    assert len(io.apply_pass_filter(df, "Mutect")) == 0


def test_pass_filter_does_not_mutate_input():
    df = variants()
    io.apply_pass_filter(df, "Mutect")
    assert len(df) == 4


def test_filters_control_is_pass_then_snp():
    assert list(io.filters_control(variants(), "Mutect")["POS"]) == [100]
    assert list(io.filters_control(variants(), "SomaticSniper")["POS"]) == [100, 200]


def test_remove_indels():
    df = variants(REF=["A", "AT", "G", "C"], ALT_1=["T", "A", "GG", "T"])
    assert list(io.remove_indels(df)["POS"]) == [100, 400]


def test_variant_set_key_format_is_chrom_pos_ref_alt1():
    assert io.variant_set(variants().iloc[:2]) == {"chr1_100_A_T", "chr1_200_C_G"}
    assert io.variant_set(variants().iloc[0:0]) == set()


def test_parse_vcf_snp_and_indel_paths(monkeypatch):
    monkeypatch.setattr(io, "read_variants", lambda path: variants())
    assert io.parse_vcf("x.vcf", "Mutect") == ({"chr1_100_A_T"}, 4, 1)
    assert io.parse_vcf("x.vcf", "SomaticSniper") == ({"chr1_100_A_T", "chr1_200_C_G"}, 4, 2)
    # indel path: PASS only, no is_snp / 1-bp requirement
    assert io.parse_vcf("x.vcf", "Mutect", indels=True) == (
        {"chr1_100_A_T", "chr2_300_G_A", "chr2_400_AT_A"}, 4, 3)


def test_truth_set_has_no_filter_test_but_drops_indels(monkeypatch):
    monkeypatch.setattr(io, "read_variants", lambda path: variants(FILTER_PASS=[False] * 4))
    assert io.truth_set("hc.vcf") == {"chr1_100_A_T", "chr1_200_C_G", "chr2_300_G_A"}


def test_validated_sets_drops_empty(monkeypatch):
    frames = {"full.vcf": variants(), "empty.vcf": variants().iloc[0:0]}
    monkeypatch.setattr(io, "read_variants", lambda path: frames[path])
    out = io.validated_sets({"k1": "full.vcf", "k2": "empty.vcf"})
    assert list(out) == ["k1"] and len(out["k1"]) == 4


COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"


def vcf_with(tmp_path, name, header_lines, samples):
    p = tmp_path / name
    p.write_text("\n".join(["##fileformat=VCFv4.2"] + header_lines + ["\t".join([COLUMNS] + samples), "1\t5"]) + "\n")
    return str(p)


def test_vcf_header_and_tumor_sample_tolerate_non_utf8(tmp_path):
    p = tmp_path / "m.vcf"
    p.write_bytes(b"##fileformat=VCFv4.2\n##note=caf\xe9\n##normal_sample=N1\n##tumor_sample=T1\n"
                  + COLUMNS.encode() + b"\tT1\tN1\n1\t5\n")
    assert io.tumor_sample(str(p)) == "T1"
    assert io.vcf_header(str(p))[-1].endswith("\tT1\tN1")


def test_tumor_sample_header_wins_whatever_the_column_order_or_file_name(tmp_path):
    header = ["##normal_sample=N1", "##tumor_sample=T1"]
    for name, samples in [("normal-first_N1-T1.vcf", ["N1", "T1"]), ("tumour-first_T1-N1.vcf", ["T1", "N1"])]:
        assert io.tumor_sample(vcf_with(tmp_path, name, header, samples)) == "T1"


def test_tumor_sample_falls_back_to_the_column_named_tumor(tmp_path):
    for samples in (["NORMAL", "TUMOR"], ["TUMOR", "NORMAL"]):
        assert io.tumor_sample(vcf_with(tmp_path, "s.vcf", [], samples)) == "TUMOR"


def test_tumor_sample_raises_when_it_cannot_be_resolved(tmp_path):
    with pytest.raises(ValueError, match="no ##tumor_sample="):
        io.tumor_sample(vcf_with(tmp_path, "a.vcf", [], ["SRR1", "SRR2"]))  # no file-name or column-order guess
    with pytest.raises(ValueError, match="no ##tumor_sample="):
        io.tumor_sample(vcf_with(tmp_path, "b.vcf", [], ["TUMOR"]))  # NORMAL/TUMOR pair incomplete
    with pytest.raises(ValueError, match="no ##tumor_sample="):
        io.tumor_sample(vcf_with(tmp_path, "c.vcf", [], ["NORMAL", "TUMOR", "OTHER"]))
    with pytest.raises(ValueError, match="not a sample column"):
        io.tumor_sample(vcf_with(tmp_path, "d.vcf", ["##tumor_sample=T9"], ["N1", "T1"]))
    with pytest.raises(ValueError, match="not a sample column"):
        io.tumor_sample(vcf_with(tmp_path, "e.vcf", ["##tumor_sample="], ["N1", "T1"]))


@pytest.mark.parametrize("suffix", [".parquet", ".csv"])
def test_sets_cache_round_trip_and_determinism(tmp_path, suffix):
    a = {"10": {"chr1_1_A_T", "chr1_2_C_G"}, "2": {"chr2_5_G_A"}}
    b = {"2": {"chr2_5_G_A"}, "10": {"chr1_2_C_G", "chr1_1_A_T"}}  # same content, different insertion order
    pa, pb = tmp_path / ("a" + suffix), tmp_path / ("b" + suffix)
    io.save_sets(a, pa)
    io.save_sets(b, pb)
    assert pa.read_bytes() == pb.read_bytes()
    assert io.load_sets(pa) == a


def test_load_sets_stringifies_integer_keys(tmp_path):
    p = tmp_path / "legacy.csv"
    pd.DataFrame({"Key": [296, 296, 606], "Identifier": ["x", "y", "z"]}).to_csv(p, index=False)
    assert io.load_sets(p) == {"296": {"x", "y"}, "606": {"z"}}


def test_write_csv_sorted_without_index(tmp_path):
    p = tmp_path / "sub" / "t.csv"
    io.write_csv(pd.DataFrame({"k": ["b", "a"], "v": [1, 2]}), p, sort_by="k")
    assert p.read_text() == "k,v\na,2\nb,1\n"


# --- BED -------------------------------------------------------------------

def bed(rows):
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def test_variants_df_round_trip():
    keys = {"chr1_100_A_T", "chr2_5_G_C"}
    df = io.variants_to_df(keys)
    assert list(df.columns) == ["#CHROM", "POS", "REF", "ALT"]
    assert df["POS"].dtype.kind == "i"
    assert io.df_to_variants(df) == keys
    assert io.df_to_variants(io.variants_to_df(set())) == set()


def test_merge_intervals_and_unique_region_size():
    assert io.merge_intervals([(5, 10), (0, 3), (8, 12), (12, 15)]) == [[0, 3], [5, 15]]
    df = bed([("chr1", 0, 10), ("chr1", 5, 15), ("chr2", 0, 10)])
    assert io.unique_region_size(df) == 25  # chr1: 0-15 = 15, chr2: 10


def test_region_size_sums_lines_without_merging(tmp_path):
    p = tmp_path / "r.bed"
    p.write_text("chr1\t0\t10\nchr1\t5\t15\n")
    assert io.region_size(str(p)) == 20  # overlap counted twice, as in legacy
    assert io.read_bed(str(p)).shape == (2, 3)
