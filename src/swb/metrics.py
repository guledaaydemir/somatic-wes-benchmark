"""Precision, recall, F1, Jaccard, ensemble calls and TMB, ported from legacy/StabilityAnalysis.ipynb.

Inputs are sets of CHROM_POS_REF_ALT1 keys. Zero-division behaviour is the legacy one:
precision and recall raise ZeroDivisionError on an empty prediction / truth set.
"""
import itertools
import operator
from collections import Counter

import numpy as np
import pandas as pd


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


def pairwise_jaccard(sets):
    """DataFrame (run_1, run_2, iou, iou_2dp): one row per unordered pair of runs, self-pairs included.

    sets is {run key: variant set}; keys must be integers or integer strings and are ordered numerically, so
    run_1 <= run_2. iou is jaccard(); iou_2dp is jaccard_2dp(), the rounded value of the legacy pair tables.
    """
    runs = sorted(sets, key=int)
    rows = [(int(a), int(b), jaccard(sets[a], sets[b]), jaccard_2dp(sets[a], sets[b]))
            for a, b in itertools.combinations_with_replacement(runs, 2)]
    return pd.DataFrame(rows, columns=["run_1", "run_2", "iou", "iou_2dp"])


def ordered_pairs(pairs):
    """Both directions of every pair: pairwise_jaccard output plus its mirror image (self-pairs once), as in the legacy tables."""
    swapped = pairs[pairs["run_1"] != pairs["run_2"]].rename(columns={"run_1": "run_2", "run_2": "run_1"})
    return pd.concat([pairs, swapped], ignore_index=True)


def legacy_pair_labels(runs):
    """The text key of each run in the legacy pair tables, as a Series indexed like runs.

    mapper[:2]/caller[:2]/isTrimmed/baseRecalibration/duplicates[:1]/TestCaseNo/sample, lower case, then "/" and
    the first letter of the environment (not lower-cased): "bw/mu/yes/no/m/296/ea/A". runs needs the columns
    Mapper, VariantCaller, isTrimmed, baseRecalibration, Duplicates, TestCaseNo, Sample and Environment.
    """
    text = (runs["Mapper"].str[:2] + "/" + runs["VariantCaller"].str[:2] + "/" + runs["isTrimmed"] + "/"
            + runs["baseRecalibration"] + "/" + runs["Duplicates"].str[:1] + "/" + runs["TestCaseNo"].astype(str)
            + "/" + runs["Sample"]).str.lower()
    return text + "/" + runs["Environment"].str[0]


def factor_contrast_means(pairs, runs, contrasts, factors=None, sample="Sample"):
    """Mean IoU over ordered run pairs from one sample that differ in exactly one factor.

    pairs is pairwise_jaccard output; runs has one row per run, indexed by the run number, with the `sample`
    column and every factor column. contrasts is {parameter name: the factor that must differ}; every other
    factor in `factors` (default: the contrasted factors only, so pass all design factors when some are not
    contrasted) and the sample must be equal. Each unordered pair counts in both directions, as in the legacy
    tables. Returns (Sample, Parameter, n_pairs, mean_iou_2dp, mean_iou), one row per sample and parameter that
    has pairs, in the order of contrasts.
    """
    factors = list(dict.fromkeys(list(factors or []) + list(contrasts.values())))
    ordered = ordered_pairs(pairs)
    first = runs.loc[ordered["run_1"]].reset_index(drop=True)
    second = runs.loc[ordered["run_2"]].reset_index(drop=True)
    same = {f: (first[f] == second[f]).to_numpy() for f in factors}
    same_sample = (first[sample] == second[sample]).to_numpy()
    out = []
    for name, factor in contrasts.items():
        mask = same_sample & ~same[factor]
        for other in factors:
            if other != factor:
                mask = mask & same[other]
        picked = ordered[mask].assign(**{sample: first.loc[mask, sample].to_numpy()})
        grouped = picked.groupby(sample).agg(
            n_pairs=("iou", "size"), mean_iou_2dp=("iou_2dp", "mean"), mean_iou=("iou", "mean")).reset_index()
        grouped.insert(1, "Parameter", name)
        out.append(grouped)
    return pd.concat(out, ignore_index=True)


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


def venn_regions(sets):
    """DataFrame (region, n_sets, size): the non-empty exclusive regions of a Venn diagram of {name: set}.

    region joins the names of the sets that hold a variant with " + ", in the order of `sets`; size counts the
    variants in exactly those sets. Largest first, ties by region name.
    """
    names = list(sets)
    counts = Counter(tuple(n for n in names if v in sets[n]) for v in set().union(*sets.values()))
    rows = [(" + ".join(members), len(members), size) for members, size in counts.items()]
    return pd.DataFrame(sorted(rows, key=lambda r: (-r[2], r[0])), columns=["region", "n_sets", "size"])


def consensus_curve(sets, truth):
    """DataFrame of the consensus of `sets` scored against `truth`, for every minimum number of votes n = 1 .. len(sets).

    Columns: n; size_ge, tp_ge, precision_ge, recall_ge, f1_ge for the variants called by at least n of the sets
    (legacy: variant_dict_greater); size_eq, tp_eq, precision_eq, recall_eq, f1_eq for exactly n
    (variant_dict). tp is the number of consensus variants in truth. The scores use the arithmetic of precision(),
    recall() and f1(), so they equal those functions on ensemble_call() sets. Precision and F1 are NaN where the
    consensus is empty (the scalar functions raise there).
    """
    sets = list(sets.values()) if isinstance(sets, dict) else list(sets)
    votes = variant_counts(sets)
    support = np.fromiter(votes.values(), dtype=np.int64, count=len(votes))
    in_truth = np.fromiter((v in truth for v in votes), dtype=bool, count=len(votes))
    n_max = len(sets)
    size_eq = np.bincount(support, minlength=n_max + 1)[1:n_max + 1]
    tp_eq = np.bincount(support[in_truth], minlength=n_max + 1)[1:n_max + 1]
    size_ge, tp_ge = size_eq[::-1].cumsum()[::-1], tp_eq[::-1].cumsum()[::-1]

    def scores(size, tp):
        precision = np.full(len(size), np.nan)
        np.divide(1.0 * tp, size, out=precision, where=size > 0)
        recall = tp / float(len(truth))
        total = precision + recall
        f1_value = np.full(len(size), np.nan)
        np.divide(precision * recall, total, out=f1_value, where=(size > 0) & (total > 0))
        f1_value = np.where(size > 0, 2.0 * np.where(total > 0, f1_value, 0.0), np.nan)
        return precision, recall, f1_value

    out = {"n": np.arange(1, n_max + 1)}
    for suffix, size, tp in (("ge", size_ge, tp_ge), ("eq", size_eq, tp_eq)):
        precision, recall, f1_value = scores(size, tp)
        out.update({"size_" + suffix: size, "tp_" + suffix: tp, "precision_" + suffix: precision,
                    "recall_" + suffix: recall, "f1_" + suffix: f1_value})
    return pd.DataFrame(out)


def vote_matrix(pipeline_sets, truth):
    """(matrix, truth_columns) for subset_vote_search: pipelines x variants, and the columns of the truth variants.

    pipeline_sets is a list of variant sets (its order is the row order). The columns are the union of the sets;
    matrix[i, j] is 1.0 if pipeline i calls variant j. Truth variants that no pipeline calls have no column.
    """
    variants = sorted(set().union(*pipeline_sets))
    column = {v: j for j, v in enumerate(variants)}
    matrix = np.zeros((len(pipeline_sets), len(variants)), dtype=np.float32)
    for i, called in enumerate(pipeline_sets):
        matrix[i, [column[v] for v in called]] = 1.0
    return matrix, np.array([column[v] for v in variants if v in truth], dtype=np.int64)


def subset_vote_search(matrix, truth_columns, size, batch=2000):
    """(subsets, n_called, tp) over every combination of `size` pipelines, for every voting threshold k = 1 .. size.

    subsets is (n_subsets, size) row indices in itertools.combinations order. n_called[i, k - 1] is the number of
    variants called by at least k pipelines of subset i (legacy: ensemble_call(..., k, "ge")) and tp[i, k - 1] how many
    of them are truth variants. Votes are counted with a matrix product, so a few million subsets take seconds.
    """
    subsets = np.array(list(itertools.combinations(range(matrix.shape[0]), size)), dtype=np.int16)
    n_called = np.zeros((len(subsets), size), dtype=np.int32)
    tp = np.zeros((len(subsets), size), dtype=np.int32)
    for start in range(0, len(subsets), batch):
        rows = subsets[start:start + batch]
        selector = np.zeros((len(rows), matrix.shape[0]), dtype=np.float32)
        selector[np.arange(len(rows))[:, None], rows] = 1.0
        votes = selector @ matrix
        in_truth = votes[:, truth_columns]
        for k in range(1, size + 1):
            n_called[start:start + batch, k - 1] = (votes >= k).sum(axis=1)
            tp[start:start + batch, k - 1] = (in_truth >= k).sum(axis=1)
    return subsets, n_called, tp


def scores_from_counts(n_called, tp, n_truth):
    """(precision, recall, f1) arrays from consensus sizes and true positives, with the arithmetic of precision(), recall(), f1().

    NaN where nothing is called (the scalar functions raise there).
    """
    n_called, tp = np.asarray(n_called), np.asarray(tp)
    precision = np.full(n_called.shape, np.nan)
    np.divide(1.0 * tp, n_called, out=precision, where=n_called > 0)
    recall = tp / float(n_truth)
    total = precision + recall
    f1_value = np.full(n_called.shape, np.nan)
    np.divide(precision * recall, total, out=f1_value, where=(n_called > 0) & (total > 0))
    f1_value = np.where(n_called > 0, 2.0 * np.where(total > 0, f1_value, 0.0), np.nan)
    return precision, recall, f1_value


def best_configuration(scores):
    """(size, k, row, value) of the highest score over {size: array (n_subsets, size)}, column k - 1 being the threshold k.

    The first maximum wins: the smallest size, then the lowest row, then the lowest k. NaN scores are ignored.
    """
    best = None
    for size in sorted(scores):
        row, column = np.unravel_index(np.nanargmax(scores[size]), scores[size].shape)
        value = scores[size][row, column]
        if best is None or value > best[3]:
            best = (size, int(column) + 1, int(row), float(value))
    return best


def region_scores(predicted, truth):
    """(n_predicted, n_truth, tp, precision, recall, f1) of a call set and a truth set that were filtered alike.

    TP = |predicted & truth|; precision = TP / |predicted|; recall = TP / |truth|; F1 their harmonic mean, with the
    arithmetic of precision(), recall() and f1(). Unlike those, an empty predicted set gives NaN precision and F1
    instead of an exception, and F1 is 0 when precision and recall are both 0. An empty truth set gives NaN recall and F1.
    """
    tp = len(predicted & truth)
    precision = tp / len(predicted) if predicted else float("nan")
    recall = tp / len(truth) if truth else float("nan")
    if precision != precision or recall != recall:
        f1_value = float("nan")
    elif precision + recall == 0:
        f1_value = 0.0
    else:
        f1_value = 2.0 * ((precision * recall) / (precision + recall))
    return len(predicted), len(truth), tp, precision, recall, f1_value


def tmb(n_variants, region_size_bp):
    """Somatic variants per megabase."""
    return n_variants / (region_size_bp / 1_000_000)
