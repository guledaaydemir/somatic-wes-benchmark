"""Toy-data checks for swb.stats."""
import itertools

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from scipy.stats import ttest_rel
from statsmodels.formula.api import ols

from swb import stats


def test_formula_default_factors_has_127_terms_in_combination_order():
    formula = stats.full_factorial_formula("f1_score")
    response, rhs = formula.split(" ~ ")
    terms = rhs.split(" + ")
    assert response == "f1_score"
    assert len(terms) == 2 ** 7 - 1
    assert terms[0] == "C(isTrimmed)"
    assert terms[7] == "C(isTrimmed):C(baseRecalibration)"
    assert terms[-1] == ":".join("C({})".format(f) for f in stats.FACTORS)


def test_formula_two_factors():
    assert stats.full_factorial_formula("y", ["a", "b"]) == "y ~ C(a) + C(b) + C(a):C(b)"


def balanced():
    return pd.DataFrame({
        "a": ["x", "x", "y", "y", "x", "x", "y", "y"],
        "b": ["p", "q", "p", "q", "p", "q", "p", "q"],
        "y": [1.0, 2.0, 4.0, 7.0, 2.0, 3.0, 5.0, 6.0],
    })


def unbalanced(n=60, seed=0):
    rng = np.random.RandomState(seed)
    df = pd.DataFrame({"a": rng.choice(list("xyz"), n, p=[0.5, 0.3, 0.2]), "b": rng.choice(list("pq"), n, p=[0.7, 0.3]),
                       "c": rng.choice(list("uv"), n)})
    df["y"] = rng.normal(size=n) + (df["a"] == "x") + 0.5 * ((df["b"] == "p") & (df["c"] == "u"))
    return df


def test_fit_anova_default_is_type_two_and_partitions_total_ss():
    df = balanced()
    table = stats.fit_anova(df, "y", ["a", "b"])
    assert list(table.index) == ["C(a)", "C(b)", "C(a):C(b)", "Residual"]
    assert table["sum_sq"].sum() == pytest.approx(((df["y"] - df["y"].mean()) ** 2).sum())


def test_fit_anova_type_two_matches_statsmodels_and_differs_from_type_one_when_unbalanced():
    df = unbalanced()
    model = ols(stats.full_factorial_formula("y", ["a", "b", "c"]), data=df).fit()
    want = sm.stats.anova_lm(model, typ=2)
    got = stats.fit_anova(df, "y", ["a", "b", "c"])
    assert got["sum_sq"].values == pytest.approx(want["sum_sq"].values)
    type1 = stats.fit_anova(df, "y", ["a", "b", "c"], typ=1)
    assert np.abs(type1["sum_sq"] - got["sum_sq"]).max() > 0.1  # the type matters when unbalanced


def test_types_one_and_two_coincide_when_balanced():
    df = pd.concat([balanced().assign(c=c) for c in ("u", "v")], ignore_index=True)
    df["y"] = df["y"] + np.arange(len(df)) * 0.37 % 1
    t1 = stats.fit_anova(df, "y", ["a", "b", "c"], typ=1)
    t2 = stats.fit_anova(df, "y", ["a", "b", "c"], typ=2)
    assert t1["sum_sq"].values == pytest.approx(t2["sum_sq"].values)


@pytest.mark.parametrize("n, seed", [(60, 0), (100, 2)])
def test_type2_sum_sq_matches_statsmodels_typ2_on_unbalanced_data(n, seed):
    df = unbalanced(n, seed)
    cells = df.groupby(["a", "b", "c"]).size()
    assert len(cells) == 12 and cells.nunique() > 1  # every cell populated, sizes unequal
    model = ols(stats.full_factorial_formula("y", ["a", "b", "c"]), data=df).fit()
    want = sm.stats.anova_lm(model, typ=2)["sum_sq"].drop("Residual")
    got = stats.type2_sum_sq(model)
    assert set(got.index) == set(want.index)
    assert got.reindex(want.index).values == pytest.approx(want.values, abs=1e-9)


def test_fit_anova_saturated_model_gives_type_two_and_matches_type_one_when_balanced():
    # one run per cell of a 2 x 2 x 3 design: residual df 0, where statsmodels' own Type II fails
    cells = pd.DataFrame(list(itertools.product("xy", "pq", "uvw")), columns=["a", "b", "c"])
    cells["y"] = np.random.RandomState(1).normal(size=len(cells))
    model = ols(stats.full_factorial_formula("y", ["a", "b", "c"]), data=cells).fit()
    assert model.df_resid == 0
    t2 = stats.fit_anova(cells, "y", ["a", "b", "c"], typ=2)
    t1 = stats.fit_anova(cells, "y", ["a", "b", "c"], typ=1)
    assert t2["sum_sq"].drop("Residual").values == pytest.approx(t1["sum_sq"].drop("Residual").values, abs=1e-12)
    assert t2["sum_sq"].sum() == pytest.approx(((cells["y"] - cells["y"].mean()) ** 2).sum())


def test_fit_anova_saturated_with_an_empty_cell_raises():
    cells = pd.DataFrame(list(itertools.product("xy", "pq")), columns=["a", "b"]).iloc[:3].copy()  # (y, q) missing
    cells["y"] = [1.0, 2.0, 3.0]
    assert ols(stats.full_factorial_formula("y", ["a", "b"]), data=cells).fit().df_resid == 0
    for typ in (1, 2):
        with pytest.raises(ValueError, match="empty design cells"):
            stats.fit_anova(cells, "y", ["a", "b"], typ=typ)


def test_fit_anova_rejects_other_types():
    with pytest.raises(ValueError, match="typ"):
        stats.fit_anova(balanced(), "y", ["a", "b"], typ=3)


def factor_frame():
    rng = np.random.RandomState(4)
    return pd.DataFrame({
        "Mapper": ["BWA", "BOWTIE"] * 6,
        "baseRecalibration": ["YES", "YES", "NO", "NO"] * 3,
        "isTrimmed": ["YES", "YES", "YES", "NO", "NO", "NO"] * 2,
        "VariantCaller": ["Mutect"] * 4 + ["SomaticSniper"] * 4 + ["Strelka"] * 4,
        "Precision": rng.uniform(0.5, 1.0, 12),
    })


def test_factor_tests_are_labelled_and_match_scipy():
    df = factor_frame()
    got = stats.factor_tests(df, "Precision")
    assert list(got) == ["Mapper", "baseRecalibration", "isTrimmed", "VariantCaller (one-way ANOVA)",
                         "Mutect vs SomaticSniper", "Mutect vs Strelka", "SomaticSniper vs Strelka"]
    bwa = df.loc[df["Mapper"] == "BWA", "Precision"]
    bowtie = df.loc[df["Mapper"] == "BOWTIE", "Precision"]
    assert got["Mapper"] == pytest.approx(tuple(ttest_rel(bwa, bowtie)))
    assert got["Mapper"] != pytest.approx(got["baseRecalibration"])


def test_factor_tests_legacy_prints_baserecalibration_under_the_mapper_label():
    df = factor_frame()
    correct = stats.factor_tests(df, "Precision")
    lines = stats.factor_tests_legacy(df, "Precision").splitlines()
    t_br, p_br = correct["baseRecalibration"]
    assert lines[0] == "Mapper                    :: T-Test: t-statistic = {}, p-value = {}".format(t_br, p_br)
    assert lines[0].split("::")[1] == lines[-2].split("::")[1]  # identical to the baseRecalibration line
    assert repr(correct["Mapper"][0]) not in lines[0]
    assert stats.factor_tests_legacy(df, "Precision").count("\n") == 10  # 7 lines + 3 blank lines, as printed


def test_variance_pct_sorts_and_drops_small_sources():
    table = pd.DataFrame({"sum_sq": [10.0, 60.0, 0.5, 29.5]}, index=["A", "B", "C", "Residual"])
    out = stats.variance_pct(table)
    assert list(out.index) == ["B", "Residual", "A"]  # C = 0.5 % is not above 1
    assert out.index.name == "Source"
    assert out["Result"].tolist() == pytest.approx([60.0, 29.5, 10.0])
    assert list(stats.variance_pct(table, min_pct=0).index) == ["B", "Residual", "A", "C"]
