"""Toy-data checks for swb.metrics. Numbers are worked out by hand in the comments."""
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
