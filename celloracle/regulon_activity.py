"""Summarize full-atlas ATF3 AUCell and plot saved basal-subcluster activity z scores."""

from pathlib import Path
import argparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def summarize_loom(path):
    import h5py

    with h5py.File(path, "r") as f:
        auc = np.asarray(f["col_attrs/RegulonsAUC"][:]["ATF3(+)"], float)
        celltype = f["col_attrs/celltype"][:].astype(str)
    table = (
        pd.DataFrame({"celltype": celltype, "auc": auc})
        .groupby("celltype")
        .auc.agg(["mean", "std", "count"])
    )
    # Exact formula recovered in 18_scenic_full_downstream.py. This is an
    # atlas-standardized mean activity, not a formal regulon specificity score.
    table["z"] = (table["mean"] - auc.mean()) / auc.std(ddof=0)
    return table.reset_index()


def draw(table):
    basal = (
        table.set_index("celltype")
        .loc[[f"Basal cell-{i}" for i in range(1, 8)]]
        .sort_values("z", ascending=False)
    )
    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    x = np.arange(len(basal))
    ax.bar(x, basal.z, color="#4C72A8", width=0.72)
    ax.set_xticks(x, [name.replace("Basal cell-", "BC-") for name in basal.index])
    for pos, value in enumerate(basal.z):
        ax.text(pos, value + 0.015, f"{value:+.2f}", ha="center", fontsize=9)
    ax.set_ylim(-0.05, max(basal.z) * 1.22)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Mean ATF3 AUCell relative to atlas\n(units of per-cell atlas SD)")
    ax.set_title("ATF3 regulon activity across basal subclusters")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig, basal


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--loom",
        type=Path,
        help="Optional: recompute summaries from full-atlas AUCell loom",
    )
    parser.add_argument(
        "--table",
        type=Path,
        default=here / "saved_results/scenic_full/ATF3_RSS_full.csv",
    )
    parser.add_argument("--output", type=Path, default=here / "generated")
    parser.add_argument("--format", choices=["png", "pdf", "both"], default="both")
    args = parser.parse_args()
    table = summarize_loom(args.loom) if args.loom else pd.read_csv(args.table)
    fig, basal = draw(table)
    args.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output / "ATF3_atlas_standardized_activity.csv", index=False)
    basal.to_csv(args.output / "ATF3_basal_activity.csv")
    for ext in (["png", "pdf"] if args.format == "both" else [args.format]):
        fig.savefig(
            args.output / f"Fig_ATF3_regulon_specificity_basal_subclusters.{ext}",
            dpi=300,
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
