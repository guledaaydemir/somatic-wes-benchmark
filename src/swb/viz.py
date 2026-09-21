"""Shared plotting helpers, ported from legacy/StabilityAnalysis.ipynb."""
import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import PathPatch
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
