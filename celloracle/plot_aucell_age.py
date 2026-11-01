"""Plot ATF3 AUCell by age using the recovered per-cell calculations.

Rectangles show cell-level means and normal 95% CIs, not quartiles.
Points are orig.ident sample means. Sites from the same donor are not
independent donors. These historical tests are retained for traceability;
see README before interpreting cell/sample p-values as age evidence.
Requires ATF3_AUCELL_LOOM; writes ATF3_FIGURE_OUT (default generated/scenic_full)."""

import os
import numpy as np
import pandas as pd
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

OUT = os.environ.get("ATF3_FIGURE_OUT", "generated/scenic_full")
os.makedirs(OUT, exist_ok=True)
AUCELL = os.environ["ATF3_AUCELL_LOOM"]
BASAL_STEM = ["Basal cell-2", "Basal cell-3", "Basal cell-4", "Basal cell-6"]
Y_COL, O_COL = ("#c6dbef", "#08519c")
Y_EDGE, O_EDGE = ("#08306b", "#08306b")


def load():
    with h5py.File(AUCELL, "r") as f:
        rec = f["col_attrs/RegulonsAUC"][:]
        names = list(rec.dtype.names)
        atf3 = (
            rec["ATF3(+)"] if "ATF3(+)" in names else rec[names[names.index("ATF3(+)")]]
        )
        celltype = f["col_attrs/celltype"][:].astype(str)
        age = f["col_attrs/age_group"][:].astype(str)
        sample = f["col_attrs/orig.ident"][:].astype(str)
    return pd.DataFrame(
        {
            "atf3": np.asarray(atf3, float),
            "celltype": celltype,
            "age": age,
            "sample": sample,
        }
    )


def is_kera(ct):
    return (
        ct.str.startswith("Basal cell")
        | ct.str.startswith("Spinous cell")
        | (ct == "Granular cell")
    )


def stats_for(sub):
    """Return dict of per-cell and per-sample statistics for a Y/O subset."""
    y = sub.loc[sub.age == "Y", "atf3"].values
    o = sub.loc[sub.age == "O", "atf3"].values
    mw_cell = stats.mannwhitneyu(y, o, alternative="two-sided")
    dm = sub.groupby(["sample", "age"])["atf3"].mean().reset_index()
    dy = dm.loc[dm.age == "Y", "atf3"].values
    do = dm.loc[dm.age == "O", "atf3"].values
    mw_don = stats.mannwhitneyu(dy, do, alternative="two-sided")
    ylo, yhi = mean_ci(y)
    olo, ohi = mean_ci(o)
    return {
        "n_cells_young": len(y),
        "n_cells_old": len(o),
        "n_samples_young": len(dy),
        "n_samples_old": len(do),
        "cell_mean_young": y.mean(),
        "cell_mean_old": o.mean(),
        "cell_ci_lo_young": ylo,
        "cell_ci_hi_young": yhi,
        "cell_ci_lo_old": olo,
        "cell_ci_hi_old": ohi,
        "cell_median_young": np.median(y),
        "cell_median_old": np.median(o),
        "sample_mean_young": dy.mean(),
        "sample_mean_old": do.mean(),
        "pct_change_old_vs_young": 100 * (o.mean() - y.mean()) / y.mean(),
        "p_percell_mw": mw_cell.pvalue,
        "p_persample_mw": mw_don.pvalue,
        "_y": y,
        "_o": o,
        "_dy": dy,
        "_do": do,
        "_ci": [(ylo, yhi), (olo, ohi)],
    }


def mean_ci(x, conf=0.95):
    """95% CI of the cell-level mean (normal approx; n is in the thousands so CLT holds)."""
    m = x.mean()
    sem = x.std(ddof=1) / np.sqrt(len(x))
    z = stats.norm.ppf(0.5 + conf / 2)
    return (m - z * sem, m + z * sem)


def draw_box(ax, s, title):
    """Box = 95% CI of the cell-level mean; center line = mean; dots = per-sample means."""
    means = [s["cell_mean_young"], s["cell_mean_old"]]
    for pos, (lo, hi), m, face in zip([0, 1], s["_ci"], means, [Y_COL, O_COL]):
        ax.add_patch(
            plt.Rectangle(
                (pos - 0.275, lo),
                0.55,
                hi - lo,
                facecolor=face,
                edgecolor="#08306b",
                linewidth=1.0,
                alpha=0.95,
                zorder=3,
            )
        )
        ax.plot([pos - 0.275, pos + 0.275], [m, m], color="black", lw=1.4, zorder=4)
    rng = np.random.RandomState(0)
    for pos, vals in [(0, s["_dy"]), (1, s["_do"])]:
        jit = pos + (rng.rand(len(vals)) - 0.5) * 0.18
        ax.scatter(
            jit,
            vals,
            s=22,
            color="#222",
            edgecolor="white",
            linewidth=0.5,
            zorder=5,
            alpha=0.9,
        )
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [
            f"Young\n(n={s['n_cells_young']} cells,\n{s['n_samples_young']} samples)",
            f"Old\n(n={s['n_cells_old']} cells,\n{s['n_samples_old']} samples)",
        ],
        fontsize=8,
    )
    ax.set_ylabel("ATF3(+) regulon activity (AUCell)", fontsize=9)
    ax.set_title(
        f"{title}\nMann-Whitney  p(cell)={fmt_p(s['p_percell_mw'])}   p(sample)={s['p_persample_mw']:.2g}",
        fontsize=9,
    )
    allv = np.concatenate([s["_dy"], s["_do"], np.ravel(s["_ci"])])
    lo, hi = (allv.min(), allv.max())
    pad = 0.12 * (hi - lo)
    top = hi + pad
    ax.plot(
        [0, 0, 1, 1],
        [top, top + pad * 0.35, top + pad * 0.35, top],
        color="#333",
        lw=1.1,
    )
    ax.text(
        0.5,
        top + pad * 0.4,
        p_to_star(s["p_persample_mw"]),
        ha="center",
        va="bottom",
        fontsize=11,
    )
    ax.set_ylim(lo - pad, top + pad * 1.5)


def fmt_p(p):
    return "<1e-300" if p == 0 else f"{p:.1e}"


def p_to_star(p):
    return (
        "****"
        if p < 0.0001
        else "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    )


def main():
    df = load()
    print(f"[box] loaded {len(df)} cells")
    pops = {
        "basal": ("Basal stem cells (BC2/3/4/6)", df[df.celltype.isin(BASAL_STEM)]),
        "allkera": (
            "All keratinocytes (Basal+Spinous+Granular)",
            df[is_kera(df.celltype)],
        ),
        "allcells": ("All cells (full atlas)", df),
    }
    rows = []
    computed = {}
    for key, (title, sub) in pops.items():
        s = stats_for(sub)
        computed[key] = (title, s)
        rows.append(
            {
                "population": title,
                **{k: v for k, v in s.items() if not k.startswith("_")},
            }
        )
        print(
            f"[box] {title}: young={s['cell_mean_young']:.4f} old={s['cell_mean_old']:.4f} ({s['pct_change_old_vs_young']:+.1f}%)  p_cell={s['p_percell_mw']:.2e}  p_sample={s['p_persample_mw']:.3g}"
        )
    pd.DataFrame(rows).to_csv(
        os.path.join(OUT, "ATF3_AUCell_age_box_summary.csv"), index=False
    )
    for key, (title, s) in computed.items():
        fig, ax = plt.subplots(figsize=(3.8, 3.6))
        draw_box(ax, s, title)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"fig_ATF3_AUCell_age_box_{key}.pdf"))
        fig.savefig(os.path.join(OUT, f"fig_ATF3_AUCell_age_box_{key}.png"), dpi=200)
        plt.close(fig)
    for tag, keys, w in [
        ("3panel", ["basal", "allkera", "allcells"], 11),
        ("2panel", ["allkera", "allcells"], 7.6),
    ]:
        fig, axes = plt.subplots(1, len(keys), figsize=(w, 3.9))
        for ax, key in zip(np.atleast_1d(axes), keys):
            draw_box(ax, computed[key][1], computed[key][0])
        fig.suptitle(
            "ATF3(+) regulon activity (AUCell): Young vs Old", fontsize=11, y=1.02
        )
        fig.tight_layout()
        fig.savefig(
            os.path.join(OUT, f"fig_ATF3_AUCell_age_box_{tag}.pdf"), bbox_inches="tight"
        )
        fig.savefig(
            os.path.join(OUT, f"fig_ATF3_AUCell_age_box_{tag}.png"),
            dpi=200,
            bbox_inches="tight",
        )
        plt.close(fig)
    print("[box] wrote 3-panel + 2-panel + individual figures + summary CSV")


if __name__ == "__main__":
    main()
