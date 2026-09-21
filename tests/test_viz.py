"""Toy-data checks for swb.viz."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

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
