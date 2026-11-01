"""Recompute the gene-recovery table and plot both attached percentage summaries."""

from pathlib import Path
import argparse
import json
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

EXPRESSION_MIN = 0.10
AGE_GAP_MIN = 0.05
RECOVERY_EPSILON = 0.005


def recovery_table(source):
    d = pd.read_csv(source)
    required = ["gene", "young", "old", "sim_young_mean"]
    if not set(required).issubset(d.columns) or d.gene.duplicated().any():
        raise ValueError("Expected unique genes and young/old/sim_young_mean columns")
    valid = d[required[1:]].notna().all(axis=1)
    excluded_missing = int((~valid).sum())
    d = d.loc[valid, required].copy()
    d["age_gap"] = d.old - d.young
    d["gap_before"] = abs(d.age_gap)
    d["gap_after"] = abs(d.sim_young_mean - d.young)
    d["recovery"] = d.gap_before - d.gap_after
    d["expr_max"] = d[["young", "old"]].max(axis=1)
    retained = (
        (d.expr_max >= EXPRESSION_MIN)
        & (d.gap_before >= AGE_GAP_MIN)
        & (d.gene != "ATF3")
    )
    d["included"] = retained
    audit = d.copy()
    d = d.loc[retained].copy()
    if d.empty:
        raise ValueError("No genes passed the documented expression/age-gap filters")
    d["rec_frac"] = d.recovery / d.gap_before
    d["category"] = np.select(
        [d.recovery > RECOVERY_EPSILON, d.recovery < -RECOVERY_EPSILON],
        ["restored", "disrupted"],
        default="unchanged",
    )
    d = d.sort_values("recovery", ascending=False)
    counts = d.category.value_counts()
    restored, disrupted, unchanged = [
        int(counts.get(x, 0)) for x in ["restored", "disrupted", "unchanged"]
    ]
    summary = dict(
        n_genes=len(d),
        restored=restored,
        disrupted=disrupted,
        unchanged=unchanged,
        impacted=restored + disrupted,
        not_restored=disrupted + unchanged,
        missing_simulation_rows=excluded_missing,
        expression_min=EXPRESSION_MIN,
        age_gap_min=AGE_GAP_MIN,
        recovery_epsilon=RECOVERY_EPSILON,
    )
    return d, audit, summary


def draw(summary, kind):
    total = summary["n_genes"]
    if kind == "impacted":
        values = [summary["impacted"], summary["unchanged"]]
        labels = ["Shifted (restored or disrupted)", "Unchanged"]
        color = "#4575B4"
        title = "Age-associated modeled genes affected by ATF3 simulation"
    else:
        values = [summary["restored"], summary["not_restored"]]
        labels = ["Restored toward young", "Not restored"]
        color = "#1A9641"
        title = "Modeled gene-expression recovery toward young"
    percentages = np.array(values) * 100 / total
    fig, ax = plt.subplots(figsize=(7.2, 3.1))
    left = 0
    for i, (value, label, percent) in enumerate(zip(values, labels, percentages)):
        ax.barh(
            0,
            percent,
            left=left,
            color=[color, "#BBBBBB"][i],
            height=0.62,
            label=f"{label} ({value}, {percent:.0f}%)",
        )
        ax.text(
            left + percent / 2,
            0,
            f"{percent:.0f}%",
            ha="center",
            va="center",
            fontsize=16,
            color="white" if i == 0 else "#333333",
            weight="bold",
        )
        left += percent
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_xlabel(f"Percentage of {total} eligible modeled genes")
    ax.set_title(title, fontsize=11)
    ax.legend(
        frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28), fontsize=9
    )
    fig.tight_layout()
    return fig


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=here / "saved_results/celloracle/full_genome_dose_response.csv",
    )
    parser.add_argument("--output", type=Path, default=here / "generated")
    parser.add_argument("--format", choices=["png", "pdf", "both"], default="both")
    args = parser.parse_args()
    table, audit, summary = recovery_table(args.input)
    args.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output / "full_genome_recovery_ranked_youngmean.csv", index=False)
    audit.to_csv(args.output / "gene_inclusion_audit.csv", index=False)
    (args.output / "recovery_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    for kind, name in [
        ("impacted", "Fig_bar_age_genes_impacted"),
        ("restored", "Fig_bar_network_restoration"),
    ]:
        fig = draw(summary, kind)
        for ext in (["png", "pdf"] if args.format == "both" else [args.format]):
            fig.savefig(args.output / f"{name}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
