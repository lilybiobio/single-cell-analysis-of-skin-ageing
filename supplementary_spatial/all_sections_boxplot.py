"""Plot paired ATF3 correlations for all sections, independently of map selection."""

from pathlib import Path
import argparse
import json
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

BASAL = "Basal keratinocytes"
FIBRO = "Combined fibroblasts"


def analyze(path, manifest):
    table = pd.read_csv(path)
    pairs = table[table.celltype.isin([BASAL, FIBRO])].pivot(
        index="dataset", columns="celltype", values="rho_absolute"
    )
    expected = json.loads(manifest.read_text())["analysis_sections"]
    if set(pairs.index) != set(expected) or pairs.isna().any().any():
        raise ValueError("Every manifest section must have both correlation values.")
    pairs = pairs.loc[expected, [BASAL, FIBRO]]
    test = wilcoxon(pairs[BASAL], pairs[FIBRO], alternative="two-sided", method="auto")
    return pairs, {
        "n_sections": len(pairs),
        "wilcoxon_statistic": float(test.statistic),
        "p_value": float(test.pvalue),
        "unit": "section",
        "paired": True,
        "mean_basal": pairs[BASAL].mean(),
        "mean_fibroblast": pairs[FIBRO].mean(),
    }


def draw(pairs, result):
    fig, ax = plt.subplots(figsize=(4.4, 4.8))
    colors = ["#2F75B5", "#E07A3F"]
    for row in pairs.to_numpy():
        ax.plot([0, 1], row, color="#C6C6C6", linewidth=0.65, alpha=0.8, zorder=1)
    boxes = ax.boxplot(
        [pairs[BASAL], pairs[FIBRO]],
        positions=[0, 1],
        widths=0.38,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
    )
    for patch, color in zip(boxes["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
    for x, name, color in zip([0, 1], [BASAL, FIBRO], colors):
        ax.scatter(np.repeat(x, len(pairs)), pairs[name], color=color, s=18, zorder=3)
    ax.set_xticks([0, 1], ["Basal\nkeratinocytes", "Fibroblasts"])
    ax.set_ylabel("Spearman correlation with normalized ATF3")
    ax.set_title(
        f"All {result['n_sections']} sections\nPaired Wilcoxon p = {result['p_value']:.4f}"
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=here / "saved_results/atf3_celltype_correlations_by_section.csv",
    )
    parser.add_argument("--sections", type=Path, default=here / "sections.json")
    parser.add_argument("--output", type=Path, default=here / "generated")
    parser.add_argument("--format", choices=["png", "pdf", "both"], default="both")
    args = parser.parse_args()
    pairs, result = analyze(args.input, args.sections)
    args.output.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.output / "all_sections_paired_values.csv")
    (args.output / "all_sections_statistics.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    fig = draw(pairs, result)
    for ext in (["png", "pdf"] if args.format == "both" else [args.format]):
        fig.savefig(
            args.output / f"atf3_all_sections_correlation_boxplot.{ext}", dpi=300
        )
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
