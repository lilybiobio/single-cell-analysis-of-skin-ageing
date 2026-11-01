"""Compute section-level ATF3 correlations for every section in sections.json."""

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

FIBROBLAST_COLUMNS = [
    "c2l_APOD+ fibroblasts",
    "c2l_SFRP2+ fibroblasts",
    "c2l_POSTN+ fibroblasts",
    "c2l_PTGDS+ fibroblasts",
]
BASAL = "Basal keratinocytes"
FIBRO = "Combined fibroblasts"


def vector(values):
    return (
        values.toarray().ravel()
        if sparse.issparse(values)
        else np.asarray(values).ravel()
    )


def read_section(path):
    import scanpy as sc

    adata = sc.read_h5ad(path)
    if adata.raw is None or "ATF3" not in adata.raw.var_names:
        raise ValueError(f"Missing raw ATF3 counts: {path.name}")
    counts = vector(adata.raw[:, "ATF3"].X).astype(float)
    library = vector(adata.raw.X.sum(axis=1)).astype(float)
    if np.any(library <= 0) or not np.isfinite(library).all():
        raise ValueError(f"Invalid library sizes: {path.name}")
    columns = ["c2l_Basal keratinocytes", *FIBROBLAST_COLUMNS]
    missing = set(columns) - set(adata.obs.columns)
    if missing:
        raise ValueError(f"{path.name}: missing abundance columns {sorted(missing)}")
    values = {
        "ATF3": np.log1p(counts / library * 10_000),
        BASAL: adata.obs["c2l_Basal keratinocytes"].to_numpy(float),
        FIBRO: adata.obs[FIBROBLAST_COLUMNS].sum(axis=1).to_numpy(float),
    }
    if any(not np.isfinite(v).all() for v in values.values()):
        raise ValueError(f"Nonfinite expression/abundance values: {path.name}")
    return adata, values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", type=Path, required=True)
    parser.add_argument(
        "--sections", type=Path, default=Path(__file__).with_name("sections.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("generated/section_correlations.csv")
    )
    args = parser.parse_args()
    expected = json.loads(args.sections.read_text())["analysis_sections"]
    observed = {p.stem for p in args.datasets.glob("*.h5ad")}
    if observed != set(expected):
        raise ValueError(
            f"Dataset manifest mismatch: missing={set(expected)-observed}; extra={observed-set(expected)}"
        )
    rows = []
    for section in expected:
        adata, values = read_section(args.datasets / f"{section}.h5ad")
        for celltype in (BASAL, FIBRO):
            rho = spearmanr(values["ATF3"], values[celltype]).statistic
            if not np.isfinite(rho):
                raise ValueError(f"Undefined correlation: {section}, {celltype}")
            rows.append(
                {
                    "dataset": section,
                    "region": "Body" if section.startswith("body") else "Face",
                    "celltype": celltype,
                    "n_spots": adata.n_obs,
                    "rho_absolute": rho,
                }
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
