"""Checks for swb.audit: nothing reaches its final path unless parity passes."""
import pytest

from swb import audit


def stage_two_files(tmp_path):
    stage = audit.Stage(tmp_path / "out" / ".staging")
    (tmp_path / "out" / "tables").mkdir(parents=True)
    a = stage.path(tmp_path / "out" / "tables" / "a.csv")
    b = stage.path(tmp_path / "derived" / "sets.parquet")
    a.write_text("a")
    b.write_text("b")
    return stage


def test_stage_path_is_stable_and_keeps_suffix(tmp_path):
    stage = audit.Stage(tmp_path / ".staging")
    p = stage.path(tmp_path / "x" / "sets.parquet")
    assert p.suffix == ".parquet"
    assert stage.path(tmp_path / "x" / "sets.parquet") == p
    assert not (tmp_path / "x" / "sets.parquet").exists()


def test_parity_pass_commits_files_and_writes_report(tmp_path):
    stage = stage_two_files(tmp_path)
    parity = audit.Parity()
    parity.check("counts match", True, "480 runs")
    parity.finish(tmp_path / "out" / "audit" / "parity.txt", stage)
    assert (tmp_path / "out" / "tables" / "a.csv").read_text() == "a"
    assert (tmp_path / "derived" / "sets.parquet").read_text() == "b"
    assert not (tmp_path / "out" / ".staging").exists()
    assert (tmp_path / "out" / "audit" / "parity.txt").read_text() == (
        "parity: PASS (1 checks)\nPASS  counts match  [480 runs]\n")


def test_parity_failure_raises_and_writes_nothing_but_the_report(tmp_path):
    stage = stage_two_files(tmp_path)
    parity = audit.Parity()
    parity.check("ok check", True)
    parity.check("bad check", False, "3 runs differ")
    with pytest.raises(AssertionError, match="bad check"):
        parity.finish(tmp_path / "out" / "audit" / "parity.txt", stage)
    assert not (tmp_path / "out" / "tables" / "a.csv").exists()
    assert not (tmp_path / "derived" / "sets.parquet").exists()
    assert not (tmp_path / "out" / ".staging").exists()
    assert "FAIL  bad check  [3 runs differ]" in (tmp_path / "out" / "audit" / "parity.txt").read_text()


def test_commit_fails_if_a_staged_file_was_never_written(tmp_path):
    stage = audit.Stage(tmp_path / ".staging")
    stage.path(tmp_path / "out" / "never.csv")
    with pytest.raises(FileNotFoundError, match="never.csv"):
        stage.commit()


def test_stem_registers_png_and_pdf(tmp_path):
    stage = audit.Stage(tmp_path / ".staging")
    stem = stage.stem(tmp_path / "figs" / "plot")
    for ext in (".png", ".pdf"):
        open(stem + ext, "w").write("x")
    stage.commit()
    assert (tmp_path / "figs" / "plot.png").exists() and (tmp_path / "figs" / "plot.pdf").exists()


def test_skipped_checks_are_reported_but_do_not_fail_or_count_as_passed():
    parity = audit.Parity()
    parity.check("a", True)
    parity.skip("raw paths", "data/raw absent")
    assert parity.ok
    assert parity.report() == "parity: PASS (1 checks, 1 skipped)\nPASS  a\nSKIP  raw paths  [data/raw absent]\n"
    plain = audit.Parity()
    plain.check("a", True)
    assert plain.report() == "parity: PASS (1 checks)\nPASS  a\n"  # unchanged when nothing is skipped
