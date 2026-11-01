"""Convert exported Matrix Market counts and metadata to a raw-count pySCENIC loom.
Set SCENIC_OUT; default prefix is full. Existing loom files are not overwritten."""

import sys
import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import loompy

OUT_DIR = os.environ.get("SCENIC_OUT", "generated/scenic")


def main(out_name="full"):
    mtx_path = os.path.join(OUT_DIR, f"{out_name}_counts.mtx")
    genes_path = os.path.join(OUT_DIR, f"{out_name}_genes.tsv")
    cells_path = os.path.join(OUT_DIR, f"{out_name}_barcodes.tsv")
    meta_path = os.path.join(OUT_DIR, f"{out_name}_meta.tsv")
    loom_path = os.path.join(OUT_DIR, f"{out_name}.loom")
    print(f"[loom] reading {mtx_path}")
    M = scipy.io.mmread(mtx_path).tocsc().astype(np.int32)
    genes = [ln.strip() for ln in open(genes_path)]
    cells = [ln.strip() for ln in open(cells_path)]
    assert M.shape[0] == len(genes), f"rows {M.shape[0]} != genes {len(genes)}"
    assert M.shape[1] == len(cells), f"cols {M.shape[1]} != cells {len(cells)}"
    print(
        f"[loom] matrix {M.shape[0]} genes x {M.shape[1]} cells, nnz={M.nnz}, dtype={M.dtype}"
    )
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    meta = meta.set_index("barcode").loc[cells].reset_index()
    assert (meta["barcode"].values == np.array(cells)).all()
    print(f"[loom] meta cols: {list(meta.columns)}")
    print(f"[loom] celltype counts:")
    print(meta["celltype"].value_counts())
    print(f"[loom] age_group x site:")
    print(pd.crosstab(meta["age_group"], meta["site"]))
    row_attrs = {"Gene": np.array(genes, dtype=object)}
    col_attrs = {
        "CellID": np.array(cells, dtype=object),
        "celltype": meta["celltype"].astype(str).values,
        "age_group": meta["age_group"].astype(str).values,
        "site": meta["site"].astype(str).values,
        "orig.ident": meta["orig.ident"].astype(str).values,
        "donor_age": meta["donor_age"].astype(np.int32).values,
        "nGene": np.array((M > 0).sum(axis=0)).flatten().astype(np.int32),
        "nUMI": np.array(M.sum(axis=0)).flatten().astype(np.int32),
    }
    if os.path.exists(loom_path):
        raise FileExistsError(
            f"Choose a new output directory; refusing to overwrite {loom_path}"
        )
    print(f"[loom] writing {loom_path}")
    loompy.create(loom_path, M.tocsc(), row_attrs, col_attrs)
    with loompy.connect(loom_path, "r") as ds:
        print(f"[loom] round-trip shape: {ds.shape}")
        print(f"[loom] row attrs: {list(ds.ra.keys())}")
        print(f"[loom] col attrs: {list(ds.ca.keys())}")
        assert "ATF3" in ds.ra["Gene"], "ATF3 missing from loom!"
        atf3_idx = int(np.where(ds.ra["Gene"] == "ATF3")[0][0])
        nnz = int((ds[atf3_idx, :] > 0).sum())
        print(f"[loom] ATF3 row idx={atf3_idx} nonzero in {nnz} cells")
    print(f"[loom] OK")


if __name__ == "__main__":
    out_name = sys.argv[1] if len(sys.argv) > 1 else "full"
    main(out_name)
