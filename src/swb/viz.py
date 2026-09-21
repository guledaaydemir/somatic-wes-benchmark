"""Shared plotting helpers, ported from legacy/StabilityAnalysis.ipynb."""
import os
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch
from scipy.cluster.hierarchy import leaves_list, linkage


def save_fig(fig, stem, dpi=300):
    """Write <stem>.png and <stem>.pdf (dpi 300, tight bbox), as the legacy cells did.

    The PDF creation date is dropped so reruns give byte-identical files.
    """
    os.makedirs(os.path.dirname(stem) or ".", exist_ok=True)
    fig.savefig(stem + ".png", dpi=dpi, bbox_inches="tight")
    fig.savefig(stem + ".pdf", dpi=dpi, bbox_inches="tight", metadata={"CreationDate": None})


def adjust_box_widths(g, fac):
    """Adjust the widths of a seaborn-generated boxplot."""
    for ax in g.axes:
        for c in ax.get_children():
            if isinstance(c, PathPatch):
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5 * (xmin + xmax)
                xhalf = 0.5 * (xmax - xmin)

                xmin_new = xmid - fac * xhalf
                xmax_new = xmid + fac * xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])


def clean_bed_name(name):
    """Figure label for a GIAB stratification BED name (legacy: last definition, cell 238).

    Earlier cells (224, 229, 230) use different labels, see docs/legacy-deviations.md D8.
    """
    key = name.lower().replace(".bed", "").replace("grch38_", "")
    special = {
        "alllowmapandsegdupregions": "All Low Map\nAnd Segdup Regions",
        "notinalllowmapandsegdupregions": "Not In All Low Map\nAnd Segdup Regions",
        "alldifficultregions": "All Difficult\nRegions",
        "notinalldifficultregions": "Not In All\nDifficult Regions",
        "findings from this study": "Findings from\nThis Study",
    }
    return special.get(key, key.replace("_", " ").title())


def cluster_pivot(df):
    """Combined_Label x Gene Symbol count matrix with rows ordered by center clustering (legacy: process_dataframe).

    Returns (pivot_df, row_colors, center_linkage, lut). Center colours follow sorted center names;
    the legacy lookup iterated a set, so colours changed between runs (docs/legacy-deviations.md D7).
    """
    grouped_data = df.groupby(["Combined_Column", "Gene Symbol", "Sample"])["Count"].sum().reset_index()
    grouped_data["Combined_Label"] = grouped_data["Combined_Column"] + " | " + grouped_data["Sample"].astype(str)

    pivot_df = grouped_data.pivot(index="Combined_Label", columns="Gene Symbol", values="Count")
    pivot_df = pivot_df.fillna(0)

    sample = grouped_data["Sample"].astype(str)
    lut = dict(zip(sorted(set(sample)), sns.color_palette("muted")))
    row_colors = pd.Series(pivot_df.index, index=pivot_df.index).map(lambda x: lut[x.split(" | ")[1]])

    center_order = pd.Series(row_colors.index).map(lambda x: x.split(" | ")[1])
    center_linkage = linkage(center_order.factorize()[0].reshape(-1, 1), method="average", metric="euclidean")
    center_ordered = leaves_list(center_linkage)

    pivot_df = pivot_df.iloc[center_ordered, :]
    row_colors = row_colors.iloc[center_ordered]
    return pivot_df, row_colors, center_linkage, lut


# Categorical slots 1-5 of the reference palette on a light surface (validated: adjacent CVD and normal-vision
# floors pass; aqua, yellow and magenta are under 3:1 contrast, so the values are also given as a table).
# Fixed by sample, never by rank: a sample keeps its colour in every figure.
SAMPLE_COLOURS = {"EA": "#2a78d6", "FD": "#eb6834", "IL": "#1baf7a", "LL": "#eda100", "NC": "#e87ba4"}
_SURFACE, _INK, _INK_2, _MUTED, _GRID, _BASELINE = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


def iou_by_factor_bars(table, value="mean_iou_2dp", title=None, subtitle=None):
    """Grouped bars of a mean IoU (0-1): one group per Parameter, one bar per Sample. Returns (fig, ax).

    table has Sample, Parameter and the `value` column (metrics.factor_contrast_means output). Groups keep the
    order of Parameter in the table; bars follow SAMPLE_COLOURS order. A Sample without a colour is an error.
    title and subtitle are set left-aligned above the plot.
    """
    unknown = sorted(set(table["Sample"]) - set(SAMPLE_COLOURS))
    if unknown:
        raise ValueError("no colour assigned to samples {}".format(unknown))
    params = list(dict.fromkeys(table["Parameter"]))
    samples = [s for s in SAMPLE_COLOURS if s in set(table["Sample"])]
    width = 0.8 / len(samples)
    fig, ax = plt.subplots(figsize=(8, 3.8), facecolor=_SURFACE)
    ax.set_facecolor(_SURFACE)
    for j, sample in enumerate(samples):
        rows = table[table["Sample"] == sample].set_index("Parameter")[value]
        heights = [rows.get(p, np.nan) for p in params]
        x = np.arange(len(params)) + (j - (len(samples) - 1) / 2) * width
        ax.bar(x, heights, width=width, color=SAMPLE_COLOURS[sample], edgecolor=_SURFACE, linewidth=1.5,
               label=sample, zorder=3)
    ax.set_xticks(range(len(params)))
    ax.set_xticklabels([textwrap.fill(p, 12, break_long_words=False) for p in params], color=_INK_2, fontsize=9)
    ax.set_ylim(0, 1.04)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"], color=_MUTED, fontsize=8)
    ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
    ax.tick_params(axis="both", length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(_BASELINE)
    ax.legend(title="Sample", ncol=len(samples), frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17),
              fontsize=9, title_fontsize=9, labelcolor=_INK_2)
    ax.get_legend().get_title().set_color(_INK_2)
    if title:
        ax.set_title(title, loc="left", color=_INK, fontsize=11, pad=22 if subtitle else 10)
    if subtitle:
        ax.text(0, 1.03, subtitle, transform=ax.transAxes, color=_INK_2, fontsize=8.5)
    return fig, ax


# Mapper colours: categorical slots 1 and 2 (the first three slots validate for scatter forms, all pairs).
MAPPER_COLOURS = {"BWA": "#2a78d6", "BOWTIE": "#eb6834"}


def run_strips(table, value, ylabel, scale=1.0, samples=None, panel_notes=None, caption=None, title=None, subtitle=None,
               seed=0):
    """One panel per Sample (shared y axis): a point per run, callers along x. Returns (fig, axes).

    table has Sample, VariantCaller, Mapper, baseRecalibration and isTrimmed ("YES" or "NO") and the column
    `value`, plotted as value / scale (ylabel names it). The y axis starts at 0 and has a tick step of 1, 2, 5, 10,
    ... chosen to give at most 8 ticks.
    Colour is the mapper (MAPPER_COLOURS). The marker gives the combination "Trimming, BaseRecalibration": circle
    for trimming YES and triangle for NO, filled for base recalibration YES and hollow for NO; the legend lists
    YES, YES / YES, NO / NO, YES / NO, NO. Points are offset by mapper and jittered horizontally with
    RandomState(seed), so reruns draw the same figure.
    samples is the panel order, left to right (default: SAMPLE_COLOURS order) and must be the samples of the
    table. panel_notes ({sample: text}) is written below each panel; caption below the legend.
    """
    unknown = sorted(set(table["Sample"]) - set(SAMPLE_COLOURS))
    if unknown:
        raise ValueError("unknown samples {}".format(unknown))
    bad = sorted(set(table["Mapper"]) - set(MAPPER_COLOURS))
    if bad:
        raise ValueError("no colour assigned to mappers {}".format(bad))
    samples = list(samples) if samples is not None else [s for s in SAMPLE_COLOURS if s in set(table["Sample"])]
    if sorted(samples) != sorted(set(table["Sample"])):
        raise ValueError("samples {} are not the samples of the table {}".format(samples, sorted(set(table["Sample"]))))
    callers = sorted(set(table["VariantCaller"]))
    rng = np.random.RandomState(seed)
    fig, axes = plt.subplots(1, len(samples), sharey=True, figsize=(11.5, 5.0), facecolor=_SURFACE, squeeze=False)
    axes = axes[0]
    for ax, sample in zip(axes, samples):
        ax.set_facecolor(_SURFACE)
        rows = table[table["Sample"] == sample]
        for mapper, colour in MAPPER_COLOURS.items():
            for trim, marker in (("YES", "o"), ("NO", "^")):
                for recal, filled in (("YES", True), ("NO", False)):
                    pick = rows[(rows["Mapper"] == mapper) & (rows["isTrimmed"] == trim)
                                & (rows["baseRecalibration"] == recal)]
                    x = (np.array([callers.index(c) for c in pick["VariantCaller"]], dtype=float)
                         + (-0.18 if mapper == "BWA" else 0.18) + rng.uniform(-0.07, 0.07, len(pick)))
                    ax.scatter(x, pick[value] / scale, s=22 if marker == "o" else 26, marker=marker,
                               linewidths=1.2, zorder=3, facecolors=colour if filled else _SURFACE, edgecolors=colour)
        ax.set_xlim(-0.6, len(callers) - 0.4)
        ax.set_xticks(range(len(callers)))
        ax.set_xticklabels([textwrap.fill(c, 8, break_long_words=False) if c != "SomaticSniper" else "Somatic\nSniper"
                            for c in callers], color=_INK_2, fontsize=7.5)
        ax.set_title(sample, loc="left", color=_INK, fontsize=10)
        ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
        ax.tick_params(axis="both", length=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(_BASELINE)
        if panel_notes and sample in panel_notes:
            ax.text(0.5, -0.2, panel_notes[sample], transform=ax.transAxes, ha="center", va="top",
                    color=_INK_2, fontsize=8.5, linespacing=1.4)
    top = float(np.ceil(table[value].max() / scale))
    step = next(st for st in (1, 2, 5, 10, 20, 50, 100) if top / st <= 8)
    axes[0].set_ylim(0, top)
    axes[0].set_yticks(range(0, int(top) + 1, step))
    axes[0].set_ylabel(ylabel, color=_INK_2, fontsize=9)
    for label in axes[0].get_yticklabels():
        label.set_color(_MUTED)
        label.set_fontsize(8)
    mapper_handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=c, markeredgecolor=c,
                                 label=m) for m, c in MAPPER_COLOURS.items()]
    combo_handles = [
        plt.Line2D([], [], marker="o" if trim == "YES" else "^", linestyle="", markersize=6,
                   markerfacecolor=_INK_2 if recal == "YES" else _SURFACE, markeredgecolor=_INK_2,
                   markeredgewidth=1.2, label="{}, {}".format(trim, recal))
        for trim in ("YES", "NO") for recal in ("YES", "NO")]
    for handles, legend_title, ncol, x in (
            (mapper_handles, "Mapper", 2, 0.27), (combo_handles, "Trimming, BaseRecalibration", 4, 0.66)):
        legend = fig.legend(handles=handles, title=legend_title, ncol=ncol, frameon=False, loc="lower center",
                            bbox_to_anchor=(x, 0.06), fontsize=9, title_fontsize=9, labelcolor=_INK_2)
        legend.get_title().set_color(_INK_2)
    fig.subplots_adjust(wspace=0.08, bottom=0.38, top=0.78 if title else 0.9)
    if caption:
        fig.text(0.075, 0.015, caption, color=_INK_2, fontsize=8.5)
    if title:
        fig.text(0.075, 0.93, title, color=_INK, fontsize=11)
    if subtitle:
        fig.text(0.075, 0.87, subtitle, color=_INK_2, fontsize=8.5)
    return fig, axes


def elapsed_time_strips(table, samples=None, panel_notes=None, caption=None, title=None, subtitle=None, seed=0):
    """run_strips of elapsed_seconds, in hours. Returns (fig, axes)."""
    return run_strips(table, "elapsed_seconds", "Elapsed time (hours)", scale=3600.0, samples=samples,
                      panel_notes=panel_notes, caption=caption, title=title, subtitle=subtitle, seed=seed)


YIELD_CAPTION = ("Centers are ordered by total sequencing yield (left to right); total bases and tumour/normal coverage "
                 "from Table 1 are given below each panel.")


def yield_panel_notes(yield_table):
    """{sample: "18.6 Gb\\n43x/58x"}: total bases and tumour/normal coverage, for run_strips panel_notes."""
    return {r.Sample: "{:.1f} Gb\n{}\u00d7/{}\u00d7".format(r.total_gb, r.tumour_coverage, r.normal_coverage)
            for r in yield_table.itertuples()}


# Centre colours of the TMB box plots: blue, orange, green, red, purple for EA, LL, NC, FD, IL (documented palette slots
# 1, 2, 3, 8, 7; the emerald of slot 3 is the green, as the dark green fails the colour-blindness check beside orange).
# The order is also the order of the boxes within a group and of the legend.
CENTER_BOX_COLOURS = {"EA": "#2a78d6", "LL": "#eb6834", "NC": "#1baf7a", "FD": "#e34948", "IL": "#4a3aa7"}
_CALLER_LABELS = {"Mutect": "MuTect2", "Strelka": "Strelka2", "SomaticSniper": "SomaticSniper"}
_MAPPER_CALLER = [(m, c) for c in ("Mutect", "Strelka", "SomaticSniper") for m in ("BWA", "BOWTIE")]
_RECAL_TRIM = [("YES", "YES"), ("NO", "NO"), ("NO", "YES"), ("YES", "NO")]


def _tint(colour, share=0.25):
    """The colour mixed with white, `share` of it: the fill of a box whose outline and median use the full colour."""
    return tuple(1 - share * (1 - c) for c in to_rgb(colour))


def tmb_box_panels(table, value="tmb", footnote=None):
    """Two box-plot panels of TMB by centre (colour, CENTER_BOX_COLOURS order). Returns (fig, axes).

    table has Sample, Mapper, VariantCaller, baseRecalibration, isTrimmed and `value`. Left panel: mapper + variant
    caller (BWA and BOWTIE with MuTect2, Strelka2, SomaticSniper); right panel: base recalibration + trimming
    (YES + YES, NO + NO, NO + YES, YES + NO). Boxes are a light tint of the centre colour with the outline, median
    line, whiskers (1.5 IQR) and the outliers beyond them in the full colour; a category and centre without runs is a gap. Raises ValueError for a centre
    without a colour, or for a run outside the categories.
    """
    unknown = sorted(set(table["Sample"]) - set(CENTER_BOX_COLOURS))
    if unknown:
        raise ValueError("no colour assigned to centres {}".format(unknown))
    pairs = set(zip(table["Mapper"], table["VariantCaller"]))
    if pairs - set(_MAPPER_CALLER):
        raise ValueError("mapper and caller outside the categories: {}".format(sorted(pairs - set(_MAPPER_CALLER))))
    combos = set(zip(table["baseRecalibration"], table["isTrimmed"]))
    if combos - set(_RECAL_TRIM):
        raise ValueError("base recalibration and trimming outside the categories: {}".format(sorted(combos - set(_RECAL_TRIM))))
    panels = [
        ("Box plot of TMB by Mapper and Variant Caller", "Mapper + Variant Caller",
         ["{} + {}".format(m, _CALLER_LABELS[c]) for m, c in _MAPPER_CALLER],
         lambda r: (r["Mapper"], r["VariantCaller"]), _MAPPER_CALLER),
        ("Box plot of TMB by Base Recalibration and Trimming", "Base Recalibration + Trimming",
         ["{} + {}".format(r, t) for r, t in _RECAL_TRIM],
         lambda r: (r["baseRecalibration"], r["isTrimmed"]), _RECAL_TRIM),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8), facecolor="#ffffff")
    top = float(np.ceil(table[value].max()))
    step = next(st for st in (1, 2, 5, 10, 20, 50, 100) if top / st <= 8)
    centers = list(CENTER_BOX_COLOURS)
    slot, half = 0.15, 0.06   # box slot width and half the box width: a gap is left between neighbouring boxes
    for ax, (title, xlabel, labels, key, keys) in zip(axes, panels):
        ax.set_facecolor("#ffffff")
        keyed = table.assign(_key=[key(r) for _, r in table.iterrows()])
        for i, k in enumerate(keys):
            for j, center in enumerate(centers):
                values = keyed.loc[(keyed["_key"] == k) & (keyed["Sample"] == center), value].to_numpy()
                if not len(values):
                    continue
                colour = CENTER_BOX_COLOURS[center]
                ax.boxplot(values, positions=[i + (j - (len(centers) - 1) / 2) * slot], widths=2 * half,
                           patch_artist=True, manage_ticks=False, zorder=3,
                           boxprops=dict(facecolor=_tint(colour), edgecolor=colour, linewidth=1.3),
                           medianprops=dict(color=colour, linewidth=2.4),
                           whiskerprops=dict(color=colour, linewidth=1.2), capprops=dict(color=colour, linewidth=1.2),
                           flierprops=dict(marker="o", markersize=3.5, markerfacecolor="none", markeredgecolor=colour))
        ax.set_title(title, loc="left", color=_INK, fontsize=11)
        ax.set_xlabel(xlabel, color=_INK_2, fontsize=10)
        ax.set_ylabel("TMB (mutations/Mb)", color=_INK_2, fontsize=10)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", color=_INK_2, fontsize=9)
        ax.set_xlim(-0.5, len(labels) - 0.5)
        ax.set_ylim(0, top)
        ax.set_yticks(range(0, int(top) + 1, step))
        ax.tick_params(axis="y", labelcolor=_MUTED, labelsize=8, length=0)
        ax.tick_params(axis="x", length=0)
        ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(_BASELINE)
    handles = [Patch(facecolor=_tint(colour), edgecolor=colour, linewidth=1.3, label=center)
               for center, colour in CENTER_BOX_COLOURS.items()]
    legend = fig.legend(handles=handles, title="Centers", loc="center left", bbox_to_anchor=(0.905, 0.55), frameon=False,
                        fontsize=10, title_fontsize=10, labelcolor=_INK_2)
    legend.get_title().set_color(_INK_2)
    fig.subplots_adjust(left=0.06, right=0.89, bottom=0.30, top=0.91, wspace=0.16)
    if footnote:
        fig.text(0.06, 0.02, footnote, color=_INK_2, fontsize=8.5)
    return fig, axes


_CALLER_MARKERS = {"Mutect": "o", "SomaticSniper": "s", "Strelka": "^"}


def precision_recall_panels(table, samples=None, panel_notes=None, caption=None, title=None, subtitle=None):
    """One panel per Sample: precision (y) against recall (x), a point per run. Returns (fig, axes).

    table has Sample, Mapper, VariantCaller, Precision and Recall. Colour is the mapper (MAPPER_COLOURS), the marker
    the variant caller (circle Mutect, square SomaticSniper, triangle Strelka). samples, panel_notes and caption
    are as in run_strips.
    """
    bad = sorted(set(table["Mapper"]) - set(MAPPER_COLOURS))
    if bad:
        raise ValueError("no colour assigned to mappers {}".format(bad))
    callers = sorted(set(table["VariantCaller"]) - set(_CALLER_MARKERS))
    if callers:
        raise ValueError("no marker assigned to callers {}".format(callers))
    samples = list(samples) if samples is not None else [s for s in SAMPLE_COLOURS if s in set(table["Sample"])]
    if sorted(samples) != sorted(set(table["Sample"])):
        raise ValueError("samples {} are not the samples of the table".format(samples))
    fig, axes = plt.subplots(1, len(samples), sharex=True, sharey=True, figsize=(11.5, 4.9), facecolor=_SURFACE,
                             squeeze=False)
    axes = axes[0]
    for ax, sample in zip(axes, samples):
        ax.set_facecolor(_SURFACE)
        rows = table[table["Sample"] == sample]
        for mapper, colour in MAPPER_COLOURS.items():
            for caller, marker in _CALLER_MARKERS.items():
                pick = rows[(rows["Mapper"] == mapper) & (rows["VariantCaller"] == caller)]
                ax.scatter(pick["Recall"], pick["Precision"], s=26, marker=marker, color=colour, edgecolors=_SURFACE,
                           linewidths=0.6, zorder=3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "0.5", "1"])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.set_xlabel("Recall", color=_INK_2, fontsize=9)
        ax.set_title(sample, loc="left", color=_INK, fontsize=10)
        ax.grid(True, color=_GRID, linewidth=0.8, zorder=0)
        ax.tick_params(axis="both", length=0, labelcolor=_MUTED, labelsize=8)
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
        if panel_notes and sample in panel_notes:
            ax.text(0.5, -0.3, panel_notes[sample], transform=ax.transAxes, ha="center", va="top", color=_INK_2,
                    fontsize=8.5, linespacing=1.4)
    axes[0].set_ylabel("Precision", color=_INK_2, fontsize=9)
    mapper_handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=6, color=c, label=m)
                      for m, c in MAPPER_COLOURS.items()]
    caller_handles = [plt.Line2D([], [], marker=m, linestyle="", markersize=6, color=_INK_2, label=c)
                      for c, m in _CALLER_MARKERS.items()]
    for handles, legend_title, ncol, x in ((mapper_handles, "Mapper", 2, 0.27), (caller_handles, "Variant caller", 3, 0.66)):
        legend = fig.legend(handles=handles, title=legend_title, ncol=ncol, frameon=False, loc="lower center",
                            bbox_to_anchor=(x, 0.06), fontsize=9, title_fontsize=9, labelcolor=_INK_2)
        legend.get_title().set_color(_INK_2)
    fig.subplots_adjust(wspace=0.12, bottom=0.36, top=0.78 if title else 0.9)
    if caption:
        fig.text(0.075, 0.015, caption, color=_INK_2, fontsize=8.5)
    if title:
        fig.text(0.075, 0.93, title, color=_INK, fontsize=11)
    if subtitle:
        fig.text(0.075, 0.87, subtitle, color=_INK_2, fontsize=8.5)
    return fig, axes


_FACTOR_LEVELS = [("Mapper", ["BWA", "BOWTIE"]), ("VariantCaller", ["Mutect", "SomaticSniper", "Strelka"]),
                  ("baseRecalibration", ["YES", "NO"]), ("isTrimmed", ["YES", "NO"])]


def p_value_text(p):
    """'< 0.001' below 0.001, else '= ' and four decimals (the legacy plot's format)."""
    return "< 0.001" if p < 0.001 else "= {:.4f}".format(p)


def factor_boxplots(table, value, tests, ylabel=None, title=None, subtitle=None):
    """Four box-plot panels of `value` by the levels of Mapper, VariantCaller, baseRecalibration and isTrimmed.

    tests is {factor: (label, p-value)} written under each panel title, e.g. ("paired t-test", 1.4e-31). Returns
    (fig, axes). The boxes carry no colour code: the levels are named on the x axis.
    """
    missing = [f for f, _ in _FACTOR_LEVELS if f not in tests]
    if missing:
        raise ValueError("no test for {}".format(missing))
    fig, axes = plt.subplots(1, 4, sharey=True, figsize=(11, 4.2), facecolor=_SURFACE,
                             gridspec_kw={"width_ratios": [2, 3, 2, 2]})
    lo, hi = float(table[value].min()), float(table[value].max())
    for ax, (factor, levels) in zip(axes, _FACTOR_LEVELS):
        ax.set_facecolor(_SURFACE)
        data = [table.loc[table[factor] == lv, value].to_numpy() for lv in levels]
        ax.boxplot(data, positions=range(len(levels)), widths=0.5, patch_artist=True, manage_ticks=False, zorder=3,
                   boxprops=dict(facecolor=_tint(SAMPLE_COLOURS["EA"]), edgecolor=SAMPLE_COLOURS["EA"], linewidth=1.3),
                   medianprops=dict(color=SAMPLE_COLOURS["EA"], linewidth=2.4),
                   whiskerprops=dict(color=SAMPLE_COLOURS["EA"], linewidth=1.2),
                   capprops=dict(color=SAMPLE_COLOURS["EA"], linewidth=1.2),
                   flierprops=dict(marker="o", markersize=3.5, markerfacecolor="none",
                                   markeredgecolor=SAMPLE_COLOURS["EA"]))
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels(["Somatic\nSniper" if lv == "SomaticSniper" else lv for lv in levels], color=_INK_2, fontsize=9)
        ax.set_xlim(-0.6, len(levels) - 0.4)
        label, p = tests[factor]
        ax.set_title(factor, loc="left", color=_INK, fontsize=10, pad=18)
        ax.text(0, 1.02, "{}, p {}".format(label, p_value_text(p)), transform=ax.transAxes, color=_INK_2, fontsize=8.5)
        ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
        ax.tick_params(axis="both", length=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(_BASELINE)
    pad = 0.05 * (hi - lo)
    axes[0].set_ylim(max(0.0, lo - pad), min(1.0, hi + pad) if hi <= 1 else hi + pad)
    axes[0].set_ylabel(ylabel or value, color=_INK_2, fontsize=9)
    for label in axes[0].get_yticklabels():
        label.set_color(_MUTED)
        label.set_fontsize(8)
    fig.subplots_adjust(wspace=0.08, bottom=0.14, top=0.72 if title else 0.8)
    if title:
        fig.text(0.075, 0.94, title, color=_INK, fontsize=11)
    if subtitle:
        fig.text(0.075, 0.88, subtitle, color=_INK_2, fontsize=8.5)
    return fig, axes


def variance_bars(percent, title=None, subtitle=None):
    """Horizontal bars of the share of the total sum of squares of each source, largest first. Returns (fig, ax).

    percent is a Series indexed by statsmodels term names ("C(isTrimmed):C(Mapper)"); the labels read
    "isTrimmed x Mapper".
    """
    percent = percent.sort_values(ascending=False)
    labels = [name.replace("C(", "").replace(")", "").replace(":", " \u00d7 ") for name in percent.index]
    fig, ax = plt.subplots(figsize=(8, 0.42 * len(percent) + 1.6), facecolor=_SURFACE)
    ax.set_facecolor(_SURFACE)
    y = np.arange(len(percent))[::-1]
    ax.barh(y, percent.to_numpy(), height=0.62, color=SAMPLE_COLOURS["EA"], zorder=3)
    for yi, v in zip(y, percent.to_numpy()):
        ax.text(v + 0.006 * percent.max(), yi, "{:.1f}".format(v), va="center", color=_INK_2, fontsize=8.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, color=_INK_2, fontsize=9)
    ax.set_xlim(0, float(percent.max()) * 1.1)
    ax.set_xlabel("Share of the total sum of squares (%)", color=_INK_2, fontsize=9)
    ax.xaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
    ax.tick_params(axis="both", length=0)
    ax.tick_params(axis="x", labelcolor=_MUTED, labelsize=8)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(_BASELINE)
    fig.subplots_adjust(left=0.42, right=0.96, top=0.86 if title else 0.96, bottom=0.14)
    if title:
        fig.text(0.02, 0.95, title, color=_INK, fontsize=11)
    if subtitle:
        fig.text(0.02, 0.905, subtitle, color=_INK_2, fontsize=8.5)
    return fig, ax


# Colours of the levels in the precision/recall box plots: blue, orange (mapper); green, red, purple (caller); brown, pink
# (YES, NO). Documented palette slots 1, 2, 3, 8, 7 and 5; brown (#a0632a) is the only value outside it, chosen because it
# passes the lightness, chroma and normal-vision checks. Every box also carries its level name below it.
PARAMETER_COLOURS = {"BWA": "#2a78d6", "BOWTIE": "#eb6834", "Mutect": "#1baf7a", "Strelka": "#e34948",
                     "SomaticSniper": "#4a3aa7", "YES": "#a0632a", "NO": "#e87ba4"}
_PARAMETERS = [
    ("Mapper", "Mapper", [("BWA", "BWA"), ("BOWTIE", "Bowtie")]),
    ("VariantCaller", "Variant Caller", [("Mutect", "MuTect2"), ("Strelka", "Strelka2"), ("SomaticSniper", "SomaticSniper")]),
    ("baseRecalibration", "Base Recalibration", [("YES", "YES"), ("NO", "NO")]),
    ("isTrimmed", "Trimming", [("YES", "YES"), ("NO", "NO")]),
]


def precision_recall_boxplots(table, p_values, title=None, footnote=None):
    """Box plots of Precision (upper panel) and Recall (lower panel) by the levels of four parameters. Returns (fig, axes).

    table has Mapper, VariantCaller, baseRecalibration, isTrimmed, Precision and Recall. p_values is
    {"Precision": {factor: p}, "Recall": {factor: p}} for the factors Mapper, VariantCaller, baseRecalibration and
    isTrimmed; each is written below its parameter, after the parameter's name. The boxes of a level cover every run
    with that level and are coloured by PARAMETER_COLOURS (YES and NO share their colours between base recalibration
    and trimming); the level names stand under the boxes and one legend below the panels names the colours.
    """
    for metric in ("Precision", "Recall"):
        missing = [f for f, _, _ in _PARAMETERS if f not in p_values.get(metric, {})]
        if missing:
            raise ValueError("no p-value for {} of {}".format(missing, metric))
    for factor, _, levels in _PARAMETERS:
        outside = sorted(set(table[factor]) - {lv for lv, _ in levels})
        if outside:
            raise ValueError("levels of {} outside the plot: {}".format(factor, outside))
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 8.4), facecolor="#ffffff")
    position, groups = 0.0, []
    for factor, name, levels in _PARAMETERS:
        xs = [position + i for i in range(len(levels))]
        groups.append((factor, name, levels, xs))
        position = xs[-1] + 1.6
    for ax, metric in zip(axes, ("Precision", "Recall")):
        ax.set_facecolor("#ffffff")
        ticks, labels = [], []
        for factor, name, levels, xs in groups:
            for (level, label), x in zip(levels, xs):
                values = table.loc[table[factor] == level, metric].to_numpy()
                colour = PARAMETER_COLOURS[level]
                ax.boxplot(values, positions=[x], widths=0.62, patch_artist=True, manage_ticks=False, zorder=3,
                           boxprops=dict(facecolor=_tint(colour), edgecolor=colour, linewidth=1.3),
                           medianprops=dict(color=colour, linewidth=2.4),
                           whiskerprops=dict(color=colour, linewidth=1.2), capprops=dict(color=colour, linewidth=1.2),
                           flierprops=dict(marker="o", markersize=3.5, markerfacecolor="none", markeredgecolor=colour))
                ticks.append(x)
                labels.append("Somatic\nSniper" if label == "SomaticSniper" else label)
            ax.text(sum(xs) / len(xs), -0.2, "{}\np {}".format(name, p_value_text(p_values[metric][factor])),
                    transform=ax.get_xaxis_transform(), ha="center", va="top", color=_INK_2, fontsize=9.5, linespacing=1.5)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, color=_INK_2, fontsize=8.5)
        ax.set_xlim(-0.7, groups[-1][3][-1] + 0.7)
        ax.set_ylim(0, 1.02)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.set_ylabel(metric, color=_INK_2, fontsize=10)
        ax.tick_params(axis="y", labelcolor=_MUTED, labelsize=8, length=0)
        ax.tick_params(axis="x", length=0)
        ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(_BASELINE)
    for handles_spec, legend_title, x in (
            ([("BWA", "BWA"), ("BOWTIE", "Bowtie")], "Mapper", 0.2),
            ([("Mutect", "MuTect2"), ("Strelka", "Strelka2"), ("SomaticSniper", "SomaticSniper")], "Variant caller", 0.5),
            ([("YES", "YES"), ("NO", "NO")], "Base recalibration, trimming", 0.8)):
        handles = [Patch(facecolor=_tint(PARAMETER_COLOURS[k]), edgecolor=PARAMETER_COLOURS[k], linewidth=1.3, label=label)
                   for k, label in handles_spec]
        legend = fig.legend(handles=handles, title=legend_title, ncol=len(handles), frameon=False, loc="lower center",
                            bbox_to_anchor=(x, 0.03), fontsize=9, title_fontsize=9, labelcolor=_INK_2)
        legend.get_title().set_color(_INK_2)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.93 if title else 0.97, bottom=0.2, hspace=0.62)
    if title:
        fig.text(0.08, 0.965, title, color=_INK, fontsize=11)
    if footnote:
        fig.text(0.08, 0.005, footnote, color=_INK_2, fontsize=8)
    return fig, axes


_SIGNIFICANCE_RED = "#e34948"


def anova_contribution_significance(table, threshold_p=0.05, title=None, subtitle=None, footnote=None):
    """Bars of each ANOVA term's contribution (%, left axis) and open circles of -log10(p) (right axis). Returns (fig, ax, ax2).

    table has label, percent and minus_log10_p, one row per term, drawn in the order of the rows. A red dashed line marks
    p = threshold_p on the right axis. The two axes have different scales: each is named after its marks in words.
    """
    for column in ("label", "percent", "minus_log10_p"):
        if column not in table.columns:
            raise ValueError("no column " + column)
    n = len(table)
    x = np.arange(n)
    fig, ax = plt.subplots(figsize=(max(9.0, 0.36 * n + 3.0), 5.8), facecolor=_SURFACE)
    ax.set_facecolor(_SURFACE)
    ax.bar(x, table["percent"].to_numpy(), width=0.62, color=SAMPLE_COLOURS["EA"], zorder=3)
    ax.set_ylabel("Contribution (%), bars", color=_INK_2, fontsize=9)
    ax.set_ylim(0, float(table["percent"].max()) * 1.08)
    ax.set_xticks(x)
    ax.set_xticklabels(table["label"], rotation=90, color=_INK_2, fontsize=8)
    ax.set_xlim(-0.7, n - 0.3)
    ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
    ax.tick_params(axis="both", length=0)
    ax.tick_params(axis="y", labelcolor=_MUTED, labelsize=8)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(_BASELINE)
    ax2 = ax.twinx()
    ax2.set_facecolor("none")
    threshold = -np.log10(threshold_p)
    top = max(float(table["minus_log10_p"].max()), threshold) * 1.08
    ax2.set_ylim(0, top)
    ax2.axhline(threshold, color=_SIGNIFICANCE_RED, linestyle="--", linewidth=1.3, zorder=4)
    ax2.scatter(x, table["minus_log10_p"].to_numpy(), s=46, facecolors="none", edgecolors=_SIGNIFICANCE_RED,
                linewidths=1.7, zorder=5)
    ax2.set_ylabel("\u2212log10(p-value), circles", color=_INK_2, fontsize=9)
    ax2.tick_params(axis="y", length=0, labelcolor=_MUTED, labelsize=8)
    for side in ("top", "left", "right", "bottom"):
        ax2.spines[side].set_visible(False)
    handles = [Patch(facecolor=SAMPLE_COLOURS["EA"], label="Contribution (%)"),
               plt.Line2D([], [], marker="o", linestyle="", markersize=7, markerfacecolor="none",
                          markeredgecolor=_SIGNIFICANCE_RED, markeredgewidth=1.7, label="\u2212log10(p-value)"),
               plt.Line2D([], [], color=_SIGNIFICANCE_RED, linestyle="--", linewidth=1.3,
                          label="p = {}".format(threshold_p))]
    fig.legend(handles=handles, ncol=3, frameon=False, loc="upper right", bbox_to_anchor=(0.95, 0.97), fontsize=9,
               labelcolor=_INK_2)
    fig.subplots_adjust(left=0.09, right=0.91, top=0.84 if title else 0.94, bottom=0.36 if footnote else 0.25)
    if title:
        fig.text(0.09, 0.95, title, color=_INK, fontsize=11)
    if subtitle:
        fig.text(0.09, 0.895, subtitle, color=_INK_2, fontsize=8.5)
    if footnote:
        fig.text(0.09, 0.02, footnote, color=_INK_2, fontsize=8, linespacing=1.5)
    return fig, ax, ax2


TRUTH_COLOUR = "#898781"  # the high-confidence truth set in the Venn figure: neutral, so the centres carry the colours


def venn_panels(data, colours, titles=("Pseudovenn", "Venn")):
    """Left a pseudovenn, right a venn of six sets, as in the legacy figure. Returns (fig, axes).

    data is {name: set} with exactly six sets and colours their fill colours in the same order (the `venn` package
    draws regions with the sizes of the exclusive parts and the names in a legend).
    """
    from venn import pseudovenn, venn

    if len(data) != 6 or len(colours) != 6:
        raise ValueError("six sets and six colours are needed, got {} and {}".format(len(data), len(colours)))
    fig, axes = plt.subplots(1, 2, figsize=(17, 8.5), facecolor="#ffffff")
    for ax, title in zip(axes, titles):
        ax.set_facecolor("#ffffff")
        ax.set_title(title, loc="left", color=_INK, fontsize=11)
    pseudovenn(data, cmap=list(colours), fontsize=9, legend_loc="upper left", ax=axes[0])
    venn(data, cmap=list(colours), fontsize=9, legend_loc="upper left", ax=axes[1])
    return fig, axes


def consensus_performance(curve, title=None, subtitle=None, footnote=None):
    """F1-score, recall and precision of the "at least n" consensus against n, with the maximum F1 marked. Returns (fig, ax).

    curve is metrics.consensus_curve output. The first maximum of F1 (the legacy idxmax) gets a ring, a dashed vertical
    line, an x tick and a label. Colours are the first three palette slots: F1 blue, recall orange, precision green.
    """
    n = curve["n"].to_numpy()
    best = curve.loc[curve["f1_ge"].idxmax()]
    fig, ax = plt.subplots(figsize=(10, 5.8), facecolor=_SURFACE)
    ax.set_facecolor(_SURFACE)
    for label, column, colour in (("F1-score", "f1_ge", "#2a78d6"), ("Recall", "recall_ge", "#eb6834"),
                                  ("Precision", "precision_ge", "#1baf7a")):
        ax.plot(n, curve[column].to_numpy(), color=colour, linewidth=1.4, marker="o", markersize=3.4,
                markeredgewidth=0, label=label, zorder=3)
    ax.axvline(best["n"], color=_INK, linestyle="--", linewidth=1.0, zorder=2)
    ax.scatter([best["n"]], [best["f1_ge"]], s=90, facecolors="none", edgecolors=_INK, linewidths=1.6, zorder=6)
    ax.annotate("Max F1-score {:.2f}\nat n = {}".format(best["f1_ge"], int(best["n"])), xy=(best["n"], best["f1_ge"]),
                xytext=(best["n"] + 14, 0.34), color=_INK, fontsize=9.5, linespacing=1.4,
                arrowprops=dict(arrowstyle="-", color=_INK, linewidth=0.9), zorder=7)
    ax.set_xlim(0, len(curve) + 1)
    ax.set_ylim(0, 1.02)
    ax.set_xticks(sorted(set(range(0, len(curve) + 1, 20)) | {int(best["n"])}))
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Minimum number of pipelines that call a variant (n)", color=_INK_2, fontsize=10)
    ax.set_ylabel("Performance metrics (F1-score, recall, precision)", color=_INK_2, fontsize=10)
    ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
    ax.tick_params(axis="both", length=0, labelcolor=_MUTED, labelsize=8.5)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(_BASELINE)
    fig.legend(*ax.get_legend_handles_labels(), ncol=3, frameon=False, loc="upper right", bbox_to_anchor=(0.97, 0.985),
               fontsize=9.5, labelcolor=_INK_2)
    fig.subplots_adjust(left=0.09, right=0.97, top=0.8 if title else 0.93, bottom=0.2 if footnote else 0.13)
    if title:
        fig.text(0.09, 0.945, title, color=_INK, fontsize=11)
    if subtitle:
        fig.text(0.09, 0.895, subtitle, color=_INK_2, fontsize=8.5)
    if footnote:
        fig.text(0.09, 0.02, footnote, color=_INK_2, fontsize=8, linespacing=1.5)
    return fig, ax


# The six GIAB stratifications in the order of the manuscript, each region followed by its complement.
REGIONS = [
    ("GRCh38_segdups", "Segmental duplications", True),
    ("GRCh38_notinsegdups", "Outside segmental duplications", False),
    ("GRCh38_alllowmapandsegdupregions", "Low mappability and segmental duplications", True),
    ("GRCh38_notinalllowmapandsegdupregions", "Outside low mappability and segmental duplications", False),
    ("GRCh38_alldifficultregions", "All difficult regions", True),
    ("GRCh38_notinalldifficultregions", "Outside all difficult regions", False),
]
REGION_SHORT = {"GRCh38_segdups": "SD", "GRCh38_notinsegdups": "Outside SD",
                "GRCh38_alllowmapandsegdupregions": "LM+SD", "GRCh38_notinalllowmapandsegdupregions": "Outside LM+SD",
                "GRCh38_alldifficultregions": "Difficult", "GRCh38_notinalldifficultregions": "Outside difficult"}
_INSIDE_COLOUR, _OUTSIDE_COLOUR = "#eb6834", "#2a78d6"   # inside a complex region, and the region's complement


def region_boxplots(table, n_truth, title=None, subtitle=None, footnote=None):
    """Box plots of precision, recall and F1 per stratification region, one panel per metric. Returns (fig, axes).

    table has region (the REGIONS names), precision, recall and f1, one row per run and region; NaN values are left out.
    n_truth is {region: number of truth variants in it}, written in the short label (REGION_SHORT: SD is segmental
    duplications, LM+SD low mappability and segmental duplications; a footnote should say so). A region and its
    complement are neighbours; the region is orange and its complement blue.
    """
    unknown = sorted(set(table["region"]) - {r for r, _, _ in REGIONS})
    if unknown:
        raise ValueError("unknown regions {}".format(unknown))
    fig, axes = plt.subplots(1, 3, sharey=True, figsize=(15, 6.3), facecolor="#ffffff")
    positions = [0, 1, 2.5, 3.5, 5, 6]
    labels = ["{} (n = {})".format(REGION_SHORT[region], n_truth[region]) for region, _, _ in REGIONS]
    for ax, (metric, name) in zip(axes, (("precision", "Precision"), ("recall", "Recall"), ("f1", "F1-score"))):
        ax.set_facecolor("#ffffff")
        for x, (region, _, inside) in zip(positions, REGIONS):
            values = table.loc[table["region"] == region, metric].dropna().to_numpy()
            if not len(values):
                continue
            colour = _INSIDE_COLOUR if inside else _OUTSIDE_COLOUR
            ax.boxplot(values, positions=[x], widths=0.7, patch_artist=True, manage_ticks=False, zorder=3,
                       boxprops=dict(facecolor=_tint(colour), edgecolor=colour, linewidth=1.3),
                       medianprops=dict(color=colour, linewidth=2.4),
                       whiskerprops=dict(color=colour, linewidth=1.2), capprops=dict(color=colour, linewidth=1.2),
                       flierprops=dict(marker="o", markersize=3.5, markerfacecolor="none", markeredgecolor=colour))
        ax.set_title(name, loc="left", color=_INK, fontsize=11)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", color=_INK_2, fontsize=8.5)
        ax.set_xlim(-0.7, 6.7)
        ax.set_ylim(0, 1.02)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.tick_params(axis="both", length=0)
        ax.tick_params(axis="y", labelcolor=_MUTED, labelsize=8.5)
        ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(_BASELINE)
    handles = [Patch(facecolor=_tint(_INSIDE_COLOUR), edgecolor=_INSIDE_COLOUR, linewidth=1.3, label="Inside the region"),
               Patch(facecolor=_tint(_OUTSIDE_COLOUR), edgecolor=_OUTSIDE_COLOUR, linewidth=1.3, label="Outside it")]
    fig.legend(handles=handles, ncol=2, frameon=False, loc="upper right", bbox_to_anchor=(0.99, 0.985), fontsize=9.5,
               labelcolor=_INK_2)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.82 if title else 0.92, bottom=0.32 if footnote else 0.24, wspace=0.06)
    if title:
        fig.text(0.05, 0.955, title, color=_INK, fontsize=11)
    if subtitle:
        fig.text(0.05, 0.905, subtitle, color=_INK_2, fontsize=8.5)
    if footnote:
        fig.text(0.05, 0.015, footnote, color=_INK_2, fontsize=8, linespacing=1.5)
    return fig, axes


# Vote thresholds 1..6 take the first six palette slots in fixed order (brown is the validated sixth).
THRESHOLD_COLOURS = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#e34948", 5: "#4a3aa7", 6: "#a0632a"}
# The six mapper-caller combinations of the individual pipelines: colour and marker shape both carry the identity.
COMBO_STYLE = [
    (("BOWTIE", "Mutect"), "Bowtie + MuTect2", "#2a78d6", "o"),
    (("BOWTIE", "SomaticSniper"), "Bowtie + SomaticSniper", "#eb6834", "s"),
    (("BOWTIE", "Strelka"), "Bowtie + Strelka2", "#1baf7a", "^"),
    (("BWA", "Mutect"), "BWA + MuTect2", "#e34948", "D"),
    (("BWA", "SomaticSniper"), "BWA + SomaticSniper", "#4a3aa7", "v"),
    (("BWA", "Strelka"), "BWA + Strelka2", "#a0632a", "P"),
]
_REFERENCE_RED = "#e34948"


def ensemble_f1_grid(individual, ensembles, reference, centers, sizes=(2, 3, 4, 5, 6), max_lines=2000, seed=0,
                     title=None, subtitle=None, footnote=None):
    """F1 of every pipeline (first column) and of every pipeline combination at every voting threshold, one row per centre.

    individual: DataFrame with Sample, Mapper, VariantCaller, isTrimmed, baseRecalibration and f1, one row per pipeline.
    ensembles: {(centre, size): array (n_combinations, size)} of F1, column k - 1 being the threshold "at least k".
    reference: {centre: F1} for the red dashed row line. Each ensemble panel marks its maximum with a red-filled marker and
    the value. Every point is drawn; a seeded sample of at most `max_lines` combinations gets the grey line across
    thresholds, and each combination has one small horizontal offset (same at every threshold). Returns (fig, axes).
    """
    rng = np.random.RandomState(seed)
    fig, axes = plt.subplots(len(centers), 1 + len(sizes), figsize=(14, 2.3 * len(centers) + 1.9), facecolor=_SURFACE,
                             sharey=True, gridspec_kw=dict(width_ratios=[1.15] + [1] * len(sizes)))
    for row, centre in enumerate(centers):
        for col in range(1 + len(sizes)):
            ax = axes[row, col]
            ax.set_facecolor(_SURFACE)
            ax.axhline(reference[centre], color=_REFERENCE_RED, linestyle="--", linewidth=1.0, zorder=1)
            if col == 0:
                rows = individual[individual["Sample"] == centre]
                for x, ((mapper, caller), _, colour, marker) in enumerate(COMBO_STYLE):
                    sub = rows[(rows["Mapper"] == mapper) & (rows["VariantCaller"] == caller)]
                    offsets = (sub["isTrimmed"].eq("YES") * 2 + sub["baseRecalibration"].eq("YES")).to_numpy() * 0.14 - 0.21
                    ax.scatter(x + offsets, sub["f1"], s=26, color=colour, marker=marker, linewidths=0, zorder=3)
                ax.set_xlim(-0.7, 5.7)
                ax.set_xticks([])
                continue
            size = sizes[col - 1]
            scores = ensembles[(centre, size)]
            offset = rng.uniform(-0.28, 0.28, size=len(scores))
            if len(scores) > max_lines:
                drawn = rng.choice(len(scores), size=max_lines, replace=False)
            else:
                drawn = np.arange(len(scores))
            xs = np.arange(1, size + 1)
            segments = [np.column_stack([xs + offset[i], scores[i]]) for i in drawn]
            ax.add_collection(LineCollection(segments, colors=_MUTED, linewidths=0.35, alpha=0.3, zorder=2, rasterized=True))
            for k in range(1, size + 1):
                ax.scatter(k + offset, scores[:, k - 1], s=2.5, color=THRESHOLD_COLOURS[k], alpha=0.35, linewidths=0,
                           zorder=3, rasterized=True)
            best = np.unravel_index(np.nanargmax(scores), scores.shape)
            ax.scatter([best[1] + 1 + offset[best[0]]], [scores[best]], s=42, color=_REFERENCE_RED, edgecolors=_INK,
                       linewidths=0.8, zorder=6)
            ax.annotate("{:.4f}".format(scores[best]), xy=(best[1] + 1 + offset[best[0]], scores[best]), xytext=(0, 7),
                        textcoords="offset points", ha="center", color=_REFERENCE_RED, fontsize=8.5, zorder=7)
            ax.set_xlim(0.4, size + 0.6)
            ax.set_xticks(xs)
        for col in range(1 + len(sizes)):
            ax = axes[row, col]
            ax.set_ylim(0, 1.04)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.tick_params(axis="both", length=0, labelcolor=_MUTED, labelsize=8.5)
            ax.yaxis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
            for side in ("top", "right", "left"):
                ax.spines[side].set_visible(False)
            ax.spines["bottom"].set_color(_BASELINE)
        axes[row, 0].set_ylabel("{}\nF1-score".format(centre), color=_INK_2, fontsize=10)
    axes[0, 0].set_title("Individual", color=_INK, fontsize=10.5)
    for col, size in enumerate(sizes, start=1):
        axes[0, col].set_title("{} pipelines".format(size), color=_INK, fontsize=10.5)
        axes[-1, col].set_xlabel("Min voting required", color=_INK_2, fontsize=9)
    axes[-1, 0].set_xlabel("Pipeline", color=_INK_2, fontsize=9)
    handles = [Line2D([], [], linestyle="", marker=marker, color=colour, markersize=6.5, label=label)
               for _, label, colour, marker in COMBO_STYLE]
    fig.legend(handles=handles, ncol=6, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.05 if footnote else 0.005),
               fontsize=9, labelcolor=_INK_2, handletextpad=0.3, columnspacing=1.6)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.915 if title else 0.94, bottom=0.125 if footnote else 0.1, wspace=0.08,
                        hspace=0.22)
    if title:
        fig.text(0.07, 0.975, title, color=_INK, fontsize=11)
    if subtitle:
        fig.text(0.07, 0.955, subtitle, color=_INK_2, fontsize=8.5)
    if footnote:
        fig.text(0.07, 0.008, footnote, color=_INK_2, fontsize=8, linespacing=1.5)
    return fig, axes
