"""Toy-data checks for swb.viz."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

from swb import viz  # noqa: E402


def make_fig():
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.plot([0, 1, 2], [0, 1, 4])
    return fig


def test_save_fig_writes_png_and_pdf_byte_identically(tmp_path):
    fig = make_fig()
    viz.save_fig(fig, str(tmp_path / "one" / "fig"))
    viz.save_fig(fig, str(tmp_path / "two" / "fig"))
    plt.close(fig)
    for ext in [".png", ".pdf"]:
        first = (tmp_path / "one" / ("fig" + ext)).read_bytes()
        assert first and first == (tmp_path / "two" / ("fig" + ext)).read_bytes()


def test_clean_bed_name():
    assert viz.clean_bed_name("GRCh38_alldifficultregions.bed") == "All Difficult\nRegions"
    assert viz.clean_bed_name("GRCh38_segdups.bed") == "Segdups"
    assert viz.clean_bed_name("Findings from this study") == "Findings from\nThis Study"


def test_cluster_pivot_shape_and_colours():
    df = pd.DataFrame({
        "Combined_Column": ["m1", "m1", "m2", "m2", "m1", "m2"],
        "Gene Symbol": ["G1", "G2", "G1", "G2", "G1", "G1"],
        "Sample": ["A", "A", "A", "A", "B", "B"],
        "Count": [1, 2, 3, 4, 5, 6],
    })
    pivot, row_colors, linkage_matrix, lut = viz.cluster_pivot(df)
    assert sorted(pivot.index) == ["m1 | A", "m1 | B", "m2 | A", "m2 | B"]
    assert list(pivot.columns) == ["G1", "G2"]
    assert pivot.loc["m1 | A"].tolist() == [1, 2]
    assert pivot.loc["m2 | B"].tolist() == [6, 0]  # missing gene filled with 0
    assert list(lut) == ["A", "B"]  # sorted center names
    assert len(row_colors) == 4 and linkage_matrix.shape == (3, 4)


def test_iou_by_factor_bars_colours_follow_the_sample_and_a_missing_value_is_a_gap():
    table = pd.DataFrame({
        "Sample": ["EA", "EA", "NC"], "Parameter": ["a", "b", "b"], "mean_iou_2dp": [0.25, 0.75, 0.5]})
    fig, ax = viz.iou_by_factor_bars(table)  # NC has no value in group "a"
    try:
        heights = [b.get_height() for b in ax.patches]
        assert heights[:2] == [0.25, 0.75] and np.isnan(heights[2]) and heights[3] == 0.5
        assert to_rgb(ax.patches[0].get_facecolor()) == to_rgb(viz.SAMPLE_COLOURS["EA"])
        assert to_rgb(ax.patches[3].get_facecolor()) == to_rgb(viz.SAMPLE_COLOURS["NC"])
        assert [t.get_text() for t in ax.get_xticklabels()] == ["a", "b"]
        assert [t.get_text() for t in ax.get_legend().get_texts()] == ["EA", "NC"]
    finally:
        plt.close(fig)


def test_iou_by_factor_bars_rejects_a_sample_without_a_colour():
    table = pd.DataFrame({"Sample": ["XX"], "Parameter": ["a"], "mean_iou_2dp": [0.5]})
    with pytest.raises(ValueError, match="no colour"):
        viz.iou_by_factor_bars(table)


def _time_table():
    return pd.DataFrame({
        "Sample": ["EA", "EA", "NC"], "VariantCaller": ["Mutect", "Strelka", "Mutect"], "Mapper": ["BWA", "BOWTIE", "BWA"],
        "baseRecalibration": ["YES", "NO", "YES"], "isTrimmed": ["YES", "NO", "NO"],
        "elapsed_seconds": [7200, 3600, 10800]})


def test_elapsed_time_strips_has_a_panel_per_sample_and_is_deterministic():
    fig1, axes1 = viz.elapsed_time_strips(_time_table())
    fig2, axes2 = viz.elapsed_time_strips(_time_table())
    try:
        assert [a.get_title(loc="left") for a in axes1] == ["EA", "NC"]
        assert axes1[0].get_ylim() == (0, 3.0)
        points = lambda axes: [c.get_offsets().tolist() for a in axes for c in a.collections]
        assert points(axes1) == points(axes2)
        heights = sorted(y for c in axes1[0].collections for _, y in c.get_offsets())
        assert heights == [1.0, 2.0]
        # trimming is the marker shape: circle (YES) and triangle (NO) points exist in panel EA
        shapes = {c.get_paths()[0].vertices.shape[0] for c in axes1[0].collections if len(c.get_offsets())}
        assert len(shapes) == 2
    finally:
        plt.close(fig1)
        plt.close(fig2)


def test_elapsed_time_strips_panel_order_and_notes():
    fig, axes = viz.elapsed_time_strips(_time_table(), samples=["NC", "EA"], panel_notes={"NC": "10.0 Gb", "EA": "20.0 Gb"},
                                        caption="a caption")
    try:
        assert [a.get_title(loc="left") for a in axes] == ["NC", "EA"]
        assert [[t.get_text() for t in a.texts] for a in axes] == [["10.0 Gb"], ["20.0 Gb"]]
        assert "a caption" in [t.get_text() for t in fig.texts]
        legends = {l.get_title().get_text(): [t.get_text() for t in l.get_texts()] for l in fig.legends}
        assert legends == {"Mapper": ["BWA", "BOWTIE"],
                           "Trimming, BaseRecalibration": ["YES, YES", "YES, NO", "NO, YES", "NO, NO"]}
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="are not the samples"):
        viz.elapsed_time_strips(_time_table(), samples=["EA"])


def test_elapsed_time_strips_keeps_its_title_and_subtitle():
    fig, _ = viz.elapsed_time_strips(_time_table(), title="The title", subtitle="The subtitle")
    try:
        assert {"The title", "The subtitle"} <= {t.get_text() for t in fig.texts}
    finally:
        plt.close(fig)


def test_run_strips_plots_any_column_on_its_own_axis():
    table = _time_table().assign(tmb=[12.0, 31.0, 18.0])
    fig, axes = viz.run_strips(table, "tmb", "TMB (variants per Mb)")
    try:
        assert axes[0].get_ylabel() == "TMB (variants per Mb)"
        assert axes[0].get_ylim() == (0, 31.0)
        assert list(axes[0].get_yticks()) == [0, 5, 10, 15, 20, 25, 30]
        assert sorted(y for c in axes[0].collections for _, y in c.get_offsets()) == [12.0, 31.0]
    finally:
        plt.close(fig)


def test_elapsed_time_strips_rejects_unknown_sample_or_mapper():
    with pytest.raises(ValueError, match="unknown samples"):
        viz.elapsed_time_strips(_time_table().assign(Sample="XX"))
    with pytest.raises(ValueError, match="no colour assigned to mappers"):
        viz.elapsed_time_strips(_time_table().assign(Mapper="HISAT"))


def test_yield_panel_notes():
    table = pd.DataFrame({"Sample": ["LL"], "total_gb": [18.6], "tumour_coverage": [43], "normal_coverage": [58]})
    assert viz.yield_panel_notes(table) == {"LL": "18.6 Gb\n43\u00d7/58\u00d7"}
    assert viz.YIELD_CAPTION.startswith("Centers are ordered by total sequencing yield")


def _tmb_table():
    rows = []
    for centre, base in (("EA", 10.0), ("LL", 20.0)):
        for i, (mapper, caller) in enumerate([("BWA", "Mutect"), ("BOWTIE", "Strelka")]):
            for k, (recal, trim) in enumerate([("YES", "YES"), ("NO", "NO"), ("NO", "YES"), ("YES", "NO")]):
                rows.append((centre, mapper, caller, recal, trim, base + i + k))
    return pd.DataFrame(rows, columns=["Sample", "Mapper", "VariantCaller", "baseRecalibration", "isTrimmed", "tmb"])


def test_tmb_box_panels_follow_the_specified_layout():
    fig, axes = viz.tmb_box_panels(_tmb_table(), footnote="a note")
    try:
        assert [a.get_title(loc="left") for a in axes] == [
            "Box plot of TMB by Mapper and Variant Caller", "Box plot of TMB by Base Recalibration and Trimming"]
        assert [a.get_xlabel() for a in axes] == ["Mapper + Variant Caller", "Base Recalibration + Trimming"]
        assert [a.get_ylabel() for a in axes] == ["TMB (mutations/Mb)"] * 2
        assert [t.get_text() for t in axes[0].get_xticklabels()] == [
            "BWA + MuTect2", "BOWTIE + MuTect2", "BWA + Strelka2", "BOWTIE + Strelka2",
            "BWA + SomaticSniper", "BOWTIE + SomaticSniper"]
        assert [t.get_text() for t in axes[1].get_xticklabels()] == ["YES + YES", "NO + NO", "NO + YES", "YES + NO"]
        assert all(t.get_rotation() == 45 for a in axes for t in a.get_xticklabels())
        legend = fig.legends[0]
        assert legend.get_title().get_text() == "Centers"
        assert [t.get_text() for t in legend.get_texts()] == ["EA", "LL", "NC", "FD", "IL"]
        # two centres x 2 mapper-caller categories drawn on the left, 2 x 4 on the right
        assert [len(a.patches) for a in axes] == [4, 8]
        assert to_rgb(axes[0].patches[0].get_edgecolor()) == to_rgb(viz.CENTER_BOX_COLOURS["EA"])  # outline in the centre colour
        assert to_rgb(axes[0].patches[0].get_facecolor()) == pytest.approx(viz._tint(viz.CENTER_BOX_COLOURS["EA"]), abs=1e-3)
        assert "a note" in [t.get_text() for t in fig.texts]
    finally:
        plt.close(fig)


def test_tmb_box_panels_reject_unknown_centres_and_runs_outside_the_categories():
    with pytest.raises(ValueError, match="no colour assigned to centres"):
        viz.tmb_box_panels(_tmb_table().assign(Sample="XX"))
    with pytest.raises(ValueError, match="mapper and caller"):
        viz.tmb_box_panels(_tmb_table().assign(Mapper="HISAT"))
    with pytest.raises(ValueError, match="base recalibration and trimming"):
        viz.tmb_box_panels(_tmb_table().assign(isTrimmed="MAYBE"))


def _metric_table():
    rows = []
    for sample in ("EA", "LL"):
        for i, mapper in enumerate(("BWA", "BOWTIE")):
            for j, caller in enumerate(("Mutect", "SomaticSniper", "Strelka")):
                for k, (recal, trim) in enumerate((("YES", "YES"), ("NO", "NO"))):
                    rows.append((sample, mapper, caller, recal, trim, 0.5 + 0.05 * (i + j), 0.4 + 0.05 * (j + k)))
    return pd.DataFrame(rows, columns=["Sample", "Mapper", "VariantCaller", "baseRecalibration", "isTrimmed",
                                       "Precision", "Recall"])


def test_precision_recall_panels_have_a_panel_per_sample_and_a_point_per_run():
    table = _metric_table()
    fig, axes = viz.precision_recall_panels(table, samples=["LL", "EA"], panel_notes={"LL": "n1", "EA": "n2"},
                                            caption="cap", title="T")
    try:
        assert [a.get_title(loc="left") for a in axes] == ["LL", "EA"]
        assert [a.get_ylabel() for a in axes] == ["Precision", ""] and all(a.get_xlabel() == "Recall" for a in axes)
        assert [sum(len(c.get_offsets()) for c in a.collections) for a in axes] == [12, 12]
        assert {l.get_title().get_text() for l in fig.legends} == {"Mapper", "Variant caller"}
        assert {"cap", "T"} <= {t.get_text() for t in fig.texts}
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="no colour assigned to mappers"):
        viz.precision_recall_panels(table.assign(Mapper="HISAT"))
    with pytest.raises(ValueError, match="no marker assigned to callers"):
        viz.precision_recall_panels(table.assign(VariantCaller="X"))


def test_p_value_text():
    assert viz.p_value_text(1.4e-31) == "< 0.001" and viz.p_value_text(0.44548) == "= 0.4455"


def test_factor_boxplots_show_the_levels_and_the_test_of_each_factor():
    tests = {"Mapper": ("paired t-test", 1e-9), "VariantCaller": ("one-way ANOVA", 0.5), "baseRecalibration": ("paired t-test", 0.02),
             "isTrimmed": ("paired t-test", 0.4455)}
    fig, axes = viz.factor_boxplots(_metric_table(), "Precision", tests, ylabel="Precision", title="T")
    try:
        assert [a.get_title(loc="left") for a in axes] == ["Mapper", "VariantCaller", "baseRecalibration", "isTrimmed"]
        assert [[t.get_text() for t in a.get_xticklabels()] for a in axes][1] == ["Mutect", "Somatic\nSniper", "Strelka"]
        assert [[t.get_text() for t in a.texts] for a in axes] == [
            ["paired t-test, p < 0.001"], ["one-way ANOVA, p = 0.5000"], ["paired t-test, p = 0.0200"],
            ["paired t-test, p = 0.4455"]]
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="no test for"):
        viz.factor_boxplots(_metric_table(), "Precision", {"Mapper": ("t", 1.0)})


def test_variance_bars_label_sources_readably_and_sort_them():
    percent = pd.Series({"C(isTrimmed):C(Mapper)": 1.8, "C(VariantCaller)": 53.9, "C(Center)": 19.4})
    fig, ax = viz.variance_bars(percent, title="T")
    try:
        labels = [t.get_text() for t in ax.get_yticklabels()]
        assert labels == ["VariantCaller", "Center", "isTrimmed \u00d7 Mapper"]  # largest first
        assert list(ax.get_yticks()) == [2, 1, 0]  # and drawn on top
        assert [round(p.get_width(), 1) for p in ax.patches] == [53.9, 19.4, 1.8]
        assert [t.get_text() for t in ax.texts] == ["53.9", "19.4", "1.8"]
    finally:
        plt.close(fig)


def _pr_p_values():
    row = {"Mapper": 1e-9, "VariantCaller": 0.5, "baseRecalibration": 0.02, "isTrimmed": 0.4455}
    return {"Precision": dict(row), "Recall": dict(row, isTrimmed=0.1)}


def test_precision_recall_boxplots_layout_colours_and_p_values():
    table = _metric_table()
    fig, axes = viz.precision_recall_boxplots(table, _pr_p_values(), title="T", footnote="a note")
    try:
        assert [a.get_ylabel() for a in axes] == ["Precision", "Recall"]
        levels = ["BWA", "Bowtie", "MuTect2", "Strelka2", "Somatic\nSniper", "YES", "NO", "YES", "NO"]
        assert [[t.get_text() for t in a.get_xticklabels()] for a in axes] == [levels, levels]
        assert [len(a.patches) for a in axes] == [9, 9]   # one box per level and parameter
        edge = [to_rgb(p.get_edgecolor()) for p in axes[0].patches]
        wanted = [viz.PARAMETER_COLOURS[k] for k in ("BWA", "BOWTIE", "Mutect", "Strelka", "SomaticSniper", "YES", "NO", "YES", "NO")]
        assert edge == [to_rgb(c) for c in wanted]
        texts = [[t.get_text() for t in a.texts] for a in axes]
        assert texts[0] == ["Mapper\np < 0.001", "Variant Caller\np = 0.5000", "Base Recalibration\np = 0.0200", "Trimming\np = 0.4455"]
        assert texts[1][3] == "Trimming\np = 0.1000"
        assert {l.get_title().get_text(): [t.get_text() for t in l.get_texts()] for l in fig.legends} == {
            "Mapper": ["BWA", "Bowtie"], "Variant caller": ["MuTect2", "Strelka2", "SomaticSniper"],
            "Base recalibration, trimming": ["YES", "NO"]}
        assert {"T", "a note"} <= {t.get_text() for t in fig.texts}
    finally:
        plt.close(fig)


def test_precision_recall_boxplots_reject_a_missing_p_value_and_unknown_levels():
    p_values = _pr_p_values()
    del p_values["Recall"]["Mapper"]
    with pytest.raises(ValueError, match="no p-value"):
        viz.precision_recall_boxplots(_metric_table(), p_values)
    with pytest.raises(ValueError, match="outside the plot"):
        viz.precision_recall_boxplots(_metric_table().assign(Mapper="HISAT"), _pr_p_values())


def test_anova_contribution_significance_draws_bars_circles_and_the_threshold():
    table = pd.DataFrame({"label": ["VC", "CEN", "MP + VC"], "percent": [54.0, 19.4, 8.0], "minus_log10_p": [900.0, 120.0, 0.5]})
    fig, ax, ax2 = viz.anova_contribution_significance(table, title="T", footnote="key")
    try:
        assert [round(p.get_height(), 1) for p in ax.patches] == [54.0, 19.4, 8.0]
        assert [t.get_text() for t in ax.get_xticklabels()] == ["VC", "CEN", "MP + VC"]
        assert ax.get_ylabel().endswith("bars") and ax2.get_ylabel().endswith("circles")
        circles = ax2.collections[0]
        assert circles.get_offsets()[:, 1].tolist() == [900.0, 120.0, 0.5]
        assert to_rgb(circles.get_edgecolor()[0]) == to_rgb("#e34948") and circles.get_facecolor().shape[0] == 0
        line = ax2.lines[0]
        assert line.get_linestyle() == "--" and list(line.get_ydata()) == pytest.approx([-np.log10(0.05)] * 2)
        assert {"T", "key"} <= {t.get_text() for t in fig.texts}
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="no column"):
        viz.anova_contribution_significance(table.drop(columns="percent"))


def test_venn_panels_need_six_sets_and_draw_two_panels():
    data = {n: {i, i + 1, 100} for i, n in enumerate("ABCDEF")}
    colours = ["#2a78d6", "#eb6834", "#1baf7a", "#e34948", "#4a3aa7", "#898781"]
    fig, axes = viz.venn_panels(data, colours)
    try:
        assert [a.get_title(loc="left") for a in axes] == ["Pseudovenn", "Venn"]
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="six sets"):
        viz.venn_panels({"A": {1}, "B": {2}}, colours)


def test_consensus_performance_marks_the_first_maximum_f1():
    curve = pd.DataFrame({
        "n": [1, 2, 3, 4], "f1_ge": [0.2, 0.9, 0.9, 0.5], "recall_ge": [1.0, 0.9, 0.8, 0.4], "precision_ge": [0.1, 0.9, 0.95, 1.0]})
    fig, ax = viz.consensus_performance(curve, title="T", footnote="note")
    try:
        assert [l.get_label() for l in ax.lines[:3]] == ["F1-score", "Recall", "Precision"]
        assert ax.lines[3].get_xdata()[0] == 2                      # the vertical line: first maximum, n = 2
        assert 2 in ax.get_xticks() and ax.get_ylim()[0] == 0
        assert [t.get_text() for t in ax.texts] == ["Max F1-score 0.90\nat n = 2"]
        assert {"T", "note"} <= {t.get_text() for t in fig.texts}
    finally:
        plt.close(fig)


def test_region_boxplots_pair_each_region_with_its_complement():
    regions = [r for r, _, _ in viz.REGIONS]
    table = pd.DataFrame([(r, 0.5 + 0.01 * i + 0.02 * j, 0.6, 0.55) for i, r in enumerate(regions) for j in range(4)],
                         columns=["region", "precision", "recall", "f1"])
    table.loc[0, "precision"] = float("nan")
    n_truth = {r: 10 * (i + 1) for i, r in enumerate(regions)}
    fig, axes = viz.region_boxplots(table, n_truth, title="T", footnote="note")
    try:
        assert [a.get_title(loc="left") for a in axes] == ["Precision", "Recall", "F1-score"]
        assert [len(a.patches) for a in axes] == [6, 6, 6]
        labels = [t.get_text() for t in axes[0].get_xticklabels()]
        assert labels[0] == "SD (n = 10)" and labels[1] == "Outside SD (n = 20)" and labels[2] == "LM+SD (n = 30)"
        edge = [to_rgb(p.get_edgecolor()) for p in axes[1].patches]
        assert edge == [to_rgb(c) for c in ["#eb6834", "#2a78d6"] * 3]
        assert {"T", "note"} <= {t.get_text() for t in fig.texts}
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="unknown regions"):
        viz.region_boxplots(table.assign(region="GRCh38_x"), n_truth)


def test_ensemble_f1_grid_has_a_row_per_centre_and_marks_each_panel_maximum():
    combos = [(m, c) for m, c, _, _ in [(x[0][0], x[0][1], 0, 0) for x in viz.COMBO_STYLE]]
    individual = pd.DataFrame([(cen, m, c, t, r, 0.5 + 0.01 * i) for cen in ("A", "B") for i, (m, c) in enumerate(combos)
                               for t in ("YES", "NO") for r in ("YES", "NO")],
                              columns=["Sample", "Mapper", "VariantCaller", "isTrimmed", "baseRecalibration", "f1"])
    rng = np.random.RandomState(0)
    ensembles = {(cen, size): rng.uniform(0.1, 0.9, size=(7, size)) for cen in ("A", "B") for size in (2, 3)}
    ensembles[("B", 3)][4, 1] = 0.95
    fig, axes = viz.ensemble_f1_grid(individual, ensembles, {"A": 0.6, "B": 0.7}, ["A", "B"], sizes=(2, 3), max_lines=3,
                                     title="T", footnote="note")
    try:
        assert axes.shape == (2, 3)
        assert [a.get_title() for a in axes[0]] == ["Individual", "2 pipelines", "3 pipelines"]
        assert [t.get_text() for t in axes[1, 2].texts] == ["0.9500"]              # the maximum of that panel, at k = 2
        assert [t.get_text() for t in axes[0, 1].texts] == ["{:.4f}".format(ensembles[("A", 2)].max())]
        assert len(axes[1, 2].collections[0].get_segments()) == 3                   # the line sample
        assert all(a.get_ylim() == (0, 1.04) for a in axes.flat)                    # one shared F1 scale
        assert {"T", "note"} <= {t.get_text() for t in fig.texts}
    finally:
        plt.close(fig)
