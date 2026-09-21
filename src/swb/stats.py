"""Factorial ANOVA helpers, ported from legacy/StabilityAnalysis.ipynb (cells 133-134)."""
import itertools

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import f_oneway, ttest_rel
from statsmodels.formula.api import ols

FACTORS = ["isTrimmed", "baseRecalibration", "Duplicates", "Mapper", "VariantCaller", "Center", "Environment"]


def full_factorial_formula(response, factors=FACTORS):
    """'response ~ C(a) + ... + C(a):C(b) + ...': every main effect and interaction, in itertools.combinations order.

    With the default FACTORS this has the same 127 terms in the same order as the hand-typed legacy formula.
    """
    terms = [
        ":".join("C({})".format(f) for f in combo)
        for k in range(1, len(factors) + 1)
        for combo in itertools.combinations(factors, k)
    ]
    return "{} ~ {}".format(response, " + ".join(terms))


def type2_sum_sq(model):
    """Type II sums of squares of a fitted formula OLS model, by nested-model residual sums of squares.

    For each term T: RSS(all terms except T and terms containing T) - RSS(the same terms plus T).
    Agrees with statsmodels' anova_lm(typ=2) wherever that runs (every design cell populated); unlike it,
    this also runs on a saturated model (statsmodels divides 0 by 0 when the residual df is 0).
    """
    info = model.model.data.design_info
    exog = np.asarray(model.model.exog, dtype=float)
    endog = np.asarray(model.model.endog, dtype=float)

    def rss(terms):
        cols = np.concatenate([np.arange(exog.shape[1])[info.term_slices[t]] for t in terms])
        beta = np.linalg.lstsq(exog[:, cols], endog, rcond=None)[0]
        resid = endog - exog[:, cols] @ beta
        return float(resid @ resid)

    out = {}
    for term in info.terms:
        if not term.factors:
            continue
        kept = [u for u in info.terms if not set(u.factors) >= set(term.factors)]
        out[term.name()] = rss(kept) - rss(kept + [term])
    return pd.Series(out)


def fit_anova(df, response, factors=FACTORS, typ=2):
    """ANOVA table (sum_sq, df, ...) of the full-factorial OLS model; Type II by default, typ=1 also allowed.

    Legacy called anova_lm(model, type=7); anova_lm ignores `type` (its keyword is `typ`), so the legacy
    table silently fell back to Type I. The design is balanced, where Type I, II and III coincide, so
    the published variance percentages are unaffected (docs/legacy-deviations.md D4).

    With one run per factor combination the model is saturated (residual df 0): F and p do not exist,
    and the Type II sums of squares come from type2_sum_sq. Takes about 20 s for 7 factors.
    """
    if typ not in (1, 2):
        raise ValueError("typ must be 1 or 2, got {!r}".format(typ))
    model = ols(full_factorial_formula(response, factors), data=df).fit()
    if model.df_resid > 0:
        return sm.stats.anova_lm(model, typ=typ)
    n_cells = int(np.prod([df[f].nunique() for f in factors]))
    if df.groupby(list(factors)).ngroups != n_cells:
        raise ValueError("empty design cells: the saturated model is rank-deficient, sums of squares are not defined")
    table = sm.stats.anova_lm(model, typ=1)
    if typ == 2:
        type2 = type2_sum_sq(model)
        terms = [name for name in table.index if name != "Residual"]
        if set(type2.index) != set(terms):
            raise ValueError("Type II terms differ from the ANOVA table terms")
        table.loc[terms, "sum_sq"] = type2[terms].values
    return table


def variance_pct(table, min_pct=1):
    """Percent of total sum of squares per source, largest first, keeping sources above min_pct."""
    total = sum(table["sum_sq"])
    out = pd.DataFrame({
        "Source": list(table.index),
        "Result": [ss / total * 100 for ss in table["sum_sq"]],
    })
    out = out.sort_values(by="Result", ascending=False, kind="mergesort").set_index("Source")
    return out.loc[out["Result"] > min_pct]


# --- factor comparisons (legacy cell 128) --------------------------------------

def _factor_stats(df, metric):
    """Every statistic of legacy cell 128, keyed by comparison. Rows are paired by position within each group.

    Mapper, baseRecalibration and isTrimmed: paired t-test of the first level against the second, on rows
    where metric is present. Callers: pairwise paired t-tests and a one-way ANOVA, on all rows.
    """
    present = ~df[metric].isnull()

    def levels(column, first, second):
        sub = df[df[column].isin([first, second]) & present]
        return sub[sub[column] == first][metric], sub[sub[column] == second][metric]

    scores = {c: df[df["VariantCaller"] == c][metric] for c in ("Mutect", "SomaticSniper", "Strelka")}
    return {
        "Mapper": ttest_rel(*levels("Mapper", "BWA", "BOWTIE")),
        "baseRecalibration": ttest_rel(*levels("baseRecalibration", "YES", "NO")),
        "isTrimmed": ttest_rel(*levels("isTrimmed", "YES", "NO")),
        "VariantCaller (one-way ANOVA)": f_oneway(scores["Mutect"], scores["SomaticSniper"], scores["Strelka"]),
        "Mutect vs SomaticSniper": ttest_rel(scores["Mutect"], scores["SomaticSniper"]),
        "Mutect vs Strelka": ttest_rel(scores["Mutect"], scores["Strelka"]),
        "SomaticSniper vs Strelka": ttest_rel(scores["SomaticSniper"], scores["Strelka"]),
    }


def factor_tests(df, metric):
    """{comparison: (statistic, p-value)}, each under its own label. Use this for new analysis.

    Same statistics as the legacy cell, including pairing rows by position (docs/legacy-deviations.md D6).
    """
    return {name: (float(res[0]), float(res[1])) for name, res in _factor_stats(df, metric).items()}


def factor_tests_legacy(df, metric):
    """The text legacy cell 128 printed, reproduced as-is. Do not use for new analysis.

    Its "Mapper" line prints the baseRecalibration statistic and p-value (D6): a Mapper value read from
    this text is a baseRecalibration comparison.
    """
    r = _factor_stats(df, metric)
    t_stat_br, p_value_ttest_br = r["baseRecalibration"]
    t_stat_it, p_value_ttest_it = r["isTrimmed"]
    t_stat, p_value = r["VariantCaller (one-way ANOVA)"]
    t_statistic_ms, p_value_ms = r["Mutect vs SomaticSniper"]
    t_statistic_mst, p_value_mst = r["Mutect vs Strelka"]
    t_statistic_st, p_value_st = r["SomaticSniper vs Strelka"]
    lines = [
        f"Mapper                    :: T-Test: t-statistic = {t_stat_br}, p-value = {p_value_ttest_br}\n",
        f"Mutect vs. SomaticSniper  :: T-Test: t-statistic = {t_statistic_ms}, p-value:{p_value_ms}",
        f"Mutect vs. Strelka        :: T-Test: t-statistic = {t_statistic_mst}, p-value:{p_value_mst}",
        f"SomaticSniper vs. Strelka :: T-Test: t-statistic = {t_statistic_st}, p-value:{p_value_st}\n",
        f"variantCaller stats.f_oneway  :: T-Test: t-statistic = {t_stat}, p-value = {p_value}\n",
        f"baseRecalibration         :: T-Test: t-statistic = {t_stat_br}, p-value = {p_value_ttest_br}",
        f"isTrimmed                 :: T-Test: t-statistic = {t_stat_it}, p-value = {p_value_ttest_it}",
    ]
    return "".join(line + "\n" for line in lines)
