"""Precision, recall, F1, Jaccard, ensemble calls and TMB, ported from legacy/StabilityAnalysis.ipynb.

Inputs are sets of CHROM_POS_REF_ALT1 keys. Zero-division behaviour is the legacy one:
precision and recall raise ZeroDivisionError on an empty prediction / truth set.
"""
import operator
from collections import Counter

import numpy as np


def true_pos(pred, real):
    return len(pred & real)


def false_pos(pred, real):
    return len(pred - real)


def false_neg(pred, real):
    return len(real - pred)


def precision(pred, real):
    """TP / (TP + FP). Raises ZeroDivisionError if pred is empty."""
    tp = true_pos(pred, real)
    fp = false_pos(pred, real)
    return (1.0 * tp) / (tp + fp)


def recall(pred, real):
    """TP / (TP + FN). Raises ZeroDivisionError if real is empty."""
    tp = true_pos(pred, real)
    fn = false_neg(pred, real)
    return (1.0 * tp) / (tp + fn)


def f1(pred, real):
    """Harmonic mean of precision and recall; 0 when both are 0."""
    p = precision(pred, real)
    r = recall(pred, real)
    if p + r == 0:
        return 0
    return 2.0 * ((p * r) / (p + r))


def prf1(pred, real):
    """(precision, recall, f1) of pred against real."""
    return precision(pred, real), recall(pred, real), f1(pred, real)


def jaccard(a, b):
    """|A n B| / |A u B|, unrounded; 0 when both are empty (legacy: intersection_over_union)."""
    a, b = set(a), set(b)
    union = len(a | b)
    if union == 0:
        return 0
    return len(a & b) / float(union)


def jaccard_2dp(a, b):
    """Jaccard rounded to 2 decimals; needs sets, raises ZeroDivisionError when both are empty (legacy: inter_union)."""
    intersection = len(a.intersection(b))
    union = len(a.union(b))
    return round(intersection / union, 2)


def variant_counts(sets):
    """{variant: number of the given sets that contain it}."""
    counts = Counter()
    for s in sets:
        counts.update(s)
    return dict(counts)


def ensemble_call(sets, min_count, mode="ge"):
    """Variants present in >= (mode 'ge'), == ('eq') or <= ('le') min_count of the given sets.

    Legacy: variant_dict_greater / variant_dict / variant_dict_lower, cell 208.
    """
    ops = {"ge": operator.ge, "eq": operator.eq, "le": operator.le}
    compare = ops[mode]
    return {v for v, c in variant_counts(sets).items() if compare(c, min_count)}


def fraction_to_count(q: float, n_lists: int) -> int:
    """Support count k = ceil(q * n_lists), clipped to [1, n_lists].

    Transfers a consensus support fraction q to n_lists call sets (for ensemble_call). Ceil, so the
    transferred rule is never weaker than the rule it was chosen under. The 1e-9 epsilon stops float
    error from lifting an exact ratio just above an integer, e.g. 7/25 * 25 = 7.000000000000001, which
    would ceil to 8 instead of 7.
    """
    if not 0.0 < q <= 1.0:
        raise ValueError(f"q={q} outside (0, 1]")
    if n_lists < 1:
        raise ValueError(f"n_lists={n_lists} must be at least 1")
    return int(min(max(int(np.ceil(q * n_lists - 1e-9)), 1), n_lists))


def tmb(n_variants, region_size_bp):
    """Somatic variants per megabase."""
    return n_variants / (region_size_bp / 1_000_000)
