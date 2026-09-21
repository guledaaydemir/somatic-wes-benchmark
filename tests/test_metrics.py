"""Toy-data checks for swb.metrics. Numbers are worked out by hand in the comments."""
import itertools

import numpy as np
import pandas as pd
import pytest

from swb import metrics

# pred = {a,b,c,d}, real = {c,d,e}: TP = {c,d} = 2, FP = {a,b} = 2, FN = {e} = 1
PRED = {"a", "b", "c", "d"}
REAL = {"c", "d", "e"}


def test_confusion_counts():
    assert metrics.true_pos(PRED, REAL) == 2
    assert metrics.false_pos(PRED, REAL) == 2
    assert metrics.false_neg(PRED, REAL) == 1


def test_precision():
    assert metrics.precision(PRED, REAL) == 0.5  # 2 / (2 + 2)


def test_recall():
    assert metrics.recall(PRED, REAL) == pytest.approx(2 / 3)  # 2 / (2 + 1)


def test_f1():
    assert metrics.f1(PRED, REAL) == pytest.approx(4 / 7)  # 2TP / (2TP + FP + FN) = 4 / 7


def test_prf1_matches_parts():
    p, r, f = metrics.prf1(PRED, REAL)
    assert (p, r, f) == (metrics.precision(PRED, REAL), metrics.recall(PRED, REAL), metrics.f1(PRED, REAL))


def test_prf1_perfect_and_disjoint():
    assert metrics.prf1({"a", "b"}, {"a", "b"}) == (1.0, 1.0, 1.0)
    assert metrics.prf1({"a"}, {"b"}) == (0.0, 0.0, 0)  # legacy f1 returns int 0 when p + r == 0


def test_precision_empty_prediction_raises():
    with pytest.raises(ZeroDivisionError):
        metrics.precision(set(), REAL)


def test_recall_empty_truth_raises():
    with pytest.raises(ZeroDivisionError):
        metrics.recall(PRED, set())


def test_jaccard():
    assert metrics.jaccard({"a", "b", "c"}, {"b", "c", "d"}) == 0.5  # 2 / 4
    assert metrics.jaccard({"a", "b", "c"}, {"c", "d", "e", "f"}) == pytest.approx(1 / 6)  # unrounded
    assert metrics.jaccard(["a", "a", "b"], ["b"]) == 0.5  # inputs are coerced to sets
    assert metrics.jaccard({"a"}, {"b"}) == 0
    assert metrics.jaccard(set(), set()) == 0


def test_jaccard_2dp_rounds_and_raises_on_empty():
    assert metrics.jaccard_2dp({"a", "b", "c"}, {"c", "d", "e", "f"}) == 0.17  # round(1/6, 2)
    assert metrics.jaccard_2dp({"a", "b"}, {"a", "b"}) == 1.0
    with pytest.raises(ZeroDivisionError):
        metrics.jaccard_2dp(set(), set())


SETS = [{"a", "b"}, {"b", "c"}, {"b"}]  # a: 1 run, b: 3 runs, c: 1 run


def test_variant_counts():
    assert metrics.variant_counts(SETS) == {"a": 1, "b": 3, "c": 1}
    assert metrics.variant_counts([]) == {}


def test_ensemble_call_modes():
    assert metrics.ensemble_call(SETS, 1) == {"a", "b", "c"}
    assert metrics.ensemble_call(SETS, 2) == {"b"}
    assert metrics.ensemble_call(SETS, 3) == {"b"}
    assert metrics.ensemble_call(SETS, 4) == set()
    assert metrics.ensemble_call(SETS, 1, mode="eq") == {"a", "c"}
    assert metrics.ensemble_call(SETS, 3, mode="eq") == {"b"}
    assert metrics.ensemble_call(SETS, 1, mode="le") == {"a", "c"}
    assert metrics.ensemble_call(SETS, 3, mode="le") == {"a", "b", "c"}


def test_ensemble_call_accepts_dict_values_and_generators():
    d = {"1": {"a", "b"}, "2": {"b"}}
    assert metrics.ensemble_call(d.values(), 2) == {"b"}
    assert metrics.ensemble_call((s for s in d.values()), 2) == {"b"}


def test_ensemble_call_unknown_mode():
    with pytest.raises(KeyError):
        metrics.ensemble_call(SETS, 1, mode="gt")


def test_ensemble_scored_against_truth():
    truth = {"b", "c"}
    # >= 2 runs -> {b}: TP 1, FP 0, FN 1
    assert metrics.prf1(metrics.ensemble_call(SETS, 2), truth) == (1.0, 0.5, pytest.approx(2 / 3))


def test_tmb():
    assert metrics.tmb(30, 30_000_000) == 1.0
    assert metrics.tmb(5, 2_500_000) == 2.0
    assert metrics.tmb(0, 1_000_000) == 0.0


# --- fraction_to_count -------------------------------------------------------

@pytest.mark.parametrize("q, n, k", [
    (34 / 120, 120, 34),  # reproduces the published threshold exactly
    (34 / 120, 96, 28),   # 27.2 -> 28
    (34 / 120, 24, 7),    # 6.8 -> 7
    (0.5, 24, 12),        # exact ratio stays 12
    (1 / 96, 24, 1),      # 0.25 -> 1
    (1.0, 24, 24),
])
def test_fraction_to_count_known_values(q, n, k):
    assert metrics.fraction_to_count(q, n) == k
    assert isinstance(metrics.fraction_to_count(q, n), int)


def test_fraction_to_count_epsilon_stops_exact_ratios_rounding_up():
    # 7/25 * 25 is 7.000000000000001 in floating point: a bare ceil would give 8.
    assert 7 / 25 * 25 > 7
    assert metrics.fraction_to_count(7 / 25, 25) == 7
    assert metrics.fraction_to_count(11 / 20, 100) == 55


def test_fraction_to_count_equals_exact_integer_ceil_on_rational_grid():
    for b in (17, 24, 25, 96, 100, 120):
        for a in range(1, b + 1):
            for n in range(1, 131):
                exact = min(max(-(-a * n // b), 1), n)  # ceil(a * n / b) in integer arithmetic
                assert metrics.fraction_to_count(a / b, n) == exact, (a, b, n)


def test_fraction_to_count_clipped_and_never_weaker_than_the_fraction():
    for n in (1, 2, 24, 96, 120):
        for q in (1e-9, 1 / 120, 0.25, 0.5, 34 / 120, 0.999, 1.0):
            k = metrics.fraction_to_count(q, n)
            assert 1 <= k <= n
            assert k >= q * n - 1e-9


def test_fraction_to_count_non_decreasing_in_q_and_in_n_lists():
    qs = sorted({a / b for b in (24, 96, 120) for a in range(1, b + 1)})
    for n in (1, 7, 24, 96, 120):
        ks = [metrics.fraction_to_count(q, n) for q in qs]
        assert ks == sorted(ks)
    for q in qs[::17] + [1.0]:
        ks = [metrics.fraction_to_count(q, n) for n in range(1, 131)]
        assert ks == sorted(ks)


@pytest.mark.parametrize("q", [0, 0.0, -0.1, 1.0000001, 2, float("nan")])
def test_fraction_to_count_rejects_q_outside_unit_interval(q):
    with pytest.raises(ValueError, match="outside"):
        metrics.fraction_to_count(q, 24)


@pytest.mark.parametrize("n", [0, -3])
def test_fraction_to_count_rejects_empty_list_count(n):
    with pytest.raises(ValueError, match="n_lists"):
        metrics.fraction_to_count(0.5, n)


# --- pairwise IoU and factor contrasts -----------------------------------------

def _toy_sets():
    return {"1": {1, 2, 3, 4}, "2": {1, 2, 3}, "3": {1, 2}, "4": {5}, "5": {1, 2, 3}, "6": {1, 2}}


def test_pairwise_jaccard_rows_order_and_values():
    got = metrics.pairwise_jaccard({"10": {"a", "b"}, "2": {"b", "c"}, "3": {"a"}})
    assert got[["run_1", "run_2"]].values.tolist() == [[2, 2], [2, 3], [2, 10], [3, 3], [3, 10], [10, 10]]
    assert got["iou"].tolist() == [1.0, 0.0, 1 / 3, 1.0, 0.5, 1.0]
    assert got["iou_2dp"].tolist() == [1.0, 0.0, 0.33, 1.0, 0.5, 1.0]
    sets = _toy_sets()
    assert len(metrics.pairwise_jaccard(sets)) == len(sets) * (len(sets) + 1) // 2


def test_pairwise_jaccard_agrees_with_the_scalar_functions():
    sets = _toy_sets()
    for r in metrics.pairwise_jaccard(sets).itertuples():
        a, b = sets[str(r.run_1)], sets[str(r.run_2)]
        assert r.iou == metrics.jaccard(a, b) and r.iou_2dp == metrics.jaccard_2dp(a, b)


def _toy_runs():
    return pd.DataFrame({
        "Sample": ["A", "A", "A", "A", "B", "B"],
        "F1": ["x", "y", "x", "y", "x", "y"],
        "F2": ["p", "p", "q", "q", "p", "p"],
    }, index=[1, 2, 3, 4, 5, 6])


def test_factor_contrast_means_counts_ordered_pairs_that_differ_in_one_factor_only():
    got = metrics.factor_contrast_means(metrics.pairwise_jaccard(_toy_sets()), _toy_runs(), {"one": "F1", "two": "F2"})
    assert got.equals(metrics.factor_contrast_means(
        metrics.pairwise_jaccard(_toy_sets()), _toy_runs(), {"one": "F1", "two": "F2"}, factors=["F1", "F2"]))
    assert got[["Sample", "Parameter", "n_pairs"]].values.tolist() == [["A", "one", 4], ["B", "one", 2], ["A", "two", 4]]
    a_one, b_one, a_two = got.iloc[0], got.iloc[1], got.iloc[2]
    # sample A, F1 differs: runs 1-2 (3/4) and 3-4 (0), each in both directions
    assert a_one["mean_iou"] == pytest.approx(0.375) and a_one["mean_iou_2dp"] == pytest.approx(0.375)
    # sample B, F1 differs: runs 5-6 (2/3, rounded 0.67)
    assert b_one["mean_iou"] == pytest.approx(2 / 3) and b_one["mean_iou_2dp"] == pytest.approx(0.67)
    # sample A, F2 differs: runs 1-3 (1/2) and 2-4 (0)
    assert a_two["mean_iou"] == pytest.approx(0.25)


def test_factor_contrast_means_ignores_pairs_from_different_samples():
    runs = _toy_runs()
    runs.loc[[5, 6], "F2"] = "q"  # B now shares F2 with A's q runs; still never paired with A
    got = metrics.factor_contrast_means(metrics.pairwise_jaccard(_toy_sets()), runs, {"one": "F1"}, factors=["F1", "F2"])
    assert got["Sample"].tolist() == ["A", "B"] and got["n_pairs"].tolist() == [4, 2]


def test_factor_contrast_means_needs_the_uncontrasted_factors_to_hold_them_equal():
    pairs = metrics.pairwise_jaccard(_toy_sets())
    only_f1 = metrics.factor_contrast_means(pairs, _toy_runs(), {"one": "F1"})  # F2 not held equal: pairs 1-4, 2-3 count
    assert only_f1["n_pairs"].tolist() == [8, 2]
    held = metrics.factor_contrast_means(pairs, _toy_runs(), {"one": "F1"}, factors=["F1", "F2"])
    assert held["n_pairs"].tolist() == [4, 2]


def test_ordered_pairs_mirrors_every_pair_but_keeps_self_pairs_once():
    pairs = metrics.pairwise_jaccard(_toy_sets())
    both = metrics.ordered_pairs(pairs)
    n = len(_toy_sets())
    assert len(both) == n * n
    lookup = {(r.run_1, r.run_2): r.iou for r in both.itertuples()}
    assert all(lookup[(b, a)] == v for (a, b), v in lookup.items())
    assert both[both["run_1"] == both["run_2"]].shape[0] == n


def test_legacy_pair_labels_lowercase_everything_but_the_environment_initial():
    runs = pd.DataFrame({
        "Mapper": ["BWA", "Bowtie"], "VariantCaller": ["Mutect", "Strelka"], "isTrimmed": ["YES", "no"],
        "baseRecalibration": ["NO", "yes"], "Duplicates": ["MARK", "DELETE"], "TestCaseNo": ["296", "61"],
        "Sample": ["EA", "LL"], "Environment": ["Altay + COSAP", "Uhem + COSAP"]}, index=[296, 61])
    got = metrics.legacy_pair_labels(runs)
    assert got.tolist() == ["bw/mu/yes/no/m/296/ea/A", "bo/st/no/yes/d/61/ll/U"]
    assert got.index.tolist() == [296, 61]


def test_venn_regions_counts_exclusive_regions():
    sets = {"A": {1, 2, 3, 4}, "B": {3, 4, 5}, "C": {4, 6}}
    got = metrics.venn_regions(sets)
    assert got.values.tolist() == [
        ["A", 1, 2], ["A + B", 2, 1], ["A + B + C", 3, 1], ["B", 1, 1], ["C", 1, 1]]
    assert got["size"].sum() == len(set().union(*sets.values()))
    assert metrics.venn_regions({"A": {1}, "B": {1}}).values.tolist() == [["A + B", 2, 1]]


def test_consensus_curve_equals_ensemble_call_scored_with_the_scalar_functions():
    rng = np.random.RandomState(5)
    universe = ["v{}".format(i) for i in range(60)]
    truth = set(universe[:25])
    sets = {str(k): set(rng.choice(universe, size=rng.randint(10, 40), replace=False)) for k in range(7)}
    curve = metrics.consensus_curve(sets, truth)
    assert curve["n"].tolist() == list(range(1, 8))
    for row in curve.itertuples():
        for mode, suffix in (("ge", "ge"), ("eq", "eq")):
            called = metrics.ensemble_call(list(sets.values()), row.n, mode)
            assert getattr(row, "size_" + suffix) == len(called)
            assert getattr(row, "tp_" + suffix) == len(called & truth)
            if called:
                p, r, f = metrics.prf1(called, truth)
                assert (getattr(row, "precision_" + suffix), getattr(row, "recall_" + suffix),
                        getattr(row, "f1_" + suffix)) == (p, r, f)


def test_consensus_curve_is_nan_where_the_consensus_is_empty_and_ge_is_monotone():
    sets = {"a": {1, 2, 3}, "b": {1, 2}, "c": {1}}
    curve = metrics.consensus_curve(sets, {1, 2, 9})
    assert curve["size_ge"].tolist() == [3, 2, 1] and curve["size_eq"].tolist() == [1, 1, 1]
    empty = metrics.consensus_curve({"a": {1}, "b": {1}, "c": {1}}, {1})
    assert empty["size_eq"].tolist() == [0, 0, 1] and empty["precision_eq"].isna().tolist() == [True, True, False]
    assert empty["f1_eq"].isna().tolist() == [True, True, False]
    assert (curve["size_ge"].diff().dropna() <= 0).all()
    zero = metrics.consensus_curve({"a": {5}}, {1})   # non-empty consensus with no true positive: F1 is 0, as f1()
    assert zero["f1_ge"].tolist() == [0.0] and zero["precision_ge"].tolist() == [0.0]


def test_region_scores_match_the_scalar_functions_and_handle_empty_sets():
    pred, real = {"a", "b", "c", "d"}, {"c", "d", "e"}
    n_pred, n_real, tp, p, r, f = metrics.region_scores(pred, real)
    assert (n_pred, n_real, tp) == (4, 3, 2) and (p, r, f) == metrics.prf1(pred, real)
    empty = metrics.region_scores(set(), real)
    assert empty[:3] == (0, 3, 0) and np.isnan(empty[3]) and empty[4] == 0.0 and np.isnan(empty[5])
    assert np.isnan(metrics.region_scores(pred, set())[4])
    none = metrics.region_scores({"x"}, real)   # no true positive: precision, recall and F1 are 0
    assert none[3:] == (0.0, 0.0, 0.0)


def test_subset_vote_search_equals_ensemble_call_scored_with_the_scalar_functions():
    rng = np.random.RandomState(11)
    universe = ["v{}".format(i) for i in range(50)]
    truth = set(universe[:20])
    pipelines = [set(rng.choice(universe, size=rng.randint(8, 30), replace=False)) for _ in range(7)]
    matrix, truth_columns = metrics.vote_matrix(pipelines, truth)
    assert matrix.shape[0] == 7 and set(matrix.flat) <= {0.0, 1.0}
    for size in (1, 2, 3, 4):
        subsets, n_called, tp = metrics.subset_vote_search(matrix, truth_columns, size, batch=5)
        assert subsets.tolist() == [list(c) for c in itertools.combinations(range(7), size)]
        precision, recall, f1 = metrics.scores_from_counts(n_called, tp, len(truth))
        for i, subset in enumerate(subsets):
            chosen = [pipelines[j] for j in subset]
            for k in range(1, size + 1):
                called = metrics.ensemble_call(chosen, k, "ge")
                assert n_called[i, k - 1] == len(called) and tp[i, k - 1] == len(called & truth)
                if called:
                    assert (precision[i, k - 1], recall[i, k - 1], f1[i, k - 1]) == metrics.prf1(called, truth)
                else:
                    assert np.isnan(precision[i, k - 1]) and np.isnan(f1[i, k - 1])


def test_scores_from_counts_handles_empty_and_zero_true_positive_consensus():
    precision, recall, f1_value = metrics.scores_from_counts(np.array([0, 4, 5]), np.array([0, 2, 0]), 10)
    assert np.isnan(precision[0]) and np.isnan(f1_value[0]) and recall[0] == 0.0
    assert precision[1] == 0.5 and recall[1] == 0.2 and f1_value[1] == pytest.approx(2 * 0.5 * 0.2 / 0.7)
    assert (precision[2], recall[2], f1_value[2]) == (0.0, 0.0, 0.0)


def test_best_configuration_takes_the_first_maximum_and_ignores_nan():
    scores = {2: np.array([[0.5, 0.7], [0.7, np.nan]]), 3: np.array([[0.1, 0.7, 0.2], [0.9, 0.9, 0.3]])}
    assert metrics.best_configuration(scores) == (3, 1, 1, 0.9)
    assert metrics.best_configuration({2: scores[2]}) == (2, 2, 0, 0.7)   # ties: lowest row, then lowest k
    assert metrics.best_configuration({3: np.array([[0.7, 0.7, 0.7]]), 2: np.array([[0.7, 0.1]])}) == (2, 1, 0, 0.7)
