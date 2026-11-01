"""Draw representative tissue maps; this file performs no hypothesis tests."""

from pathlib import Path
import argparse
import importlib.util
import json
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

source = Path(__file__).with_name("section_correlations.py")
spec = importlib.util.spec_from_file_location("section_correlations", source)
sections = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sections)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, default=Path("generated/representative_spatial_maps")
    )
    args = parser.parse_args()
    manifest = json.loads(Path(__file__).with_name("sections.json").read_text())
    representatives = manifest["representative_sections"]
    fig, axes = plt.subplots(
        len(representatives), 3, figsize=(10, 3.1 * len(representatives))
    )
    for row, (identifier, label) in enumerate(representatives.items()):
        adata, values = sections.read_section(args.datasets / f"{identifier}.h5ad")
        if "spatial" in adata.obsm:
            coords = np.asarray(adata.obsm["spatial"])
        elif "X_spatial" in adata.obsm:
            coords = np.asarray(adata.obsm["X_spatial"])
        else:
            coords = adata.obs[["array_col", "array_row"]].to_numpy()
        for column, (feature, palette) in enumerate(
            [(sections.FIBRO, "Reds"), (sections.BASAL, "Blues"), ("ATF3", "Greens")]
        ):
            ax = axes[row, column]
            points = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values[feature],
                cmap=palette,
                s=4,
                linewidths=0,
                rasterized=True,
            )
            ax.invert_yaxis()
            ax.set_aspect("equal")
            ax.set_title(f"{label}\n{feature}", fontsize=9)
            ax.set_xlabel("Spatial 1")
            ax.set_ylabel("Spatial 2")
            fig.colorbar(points, ax=ax, shrink=0.6)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(args.output.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
