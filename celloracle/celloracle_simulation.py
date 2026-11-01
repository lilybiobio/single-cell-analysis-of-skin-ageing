"""Fit the recovered CellOracle model and simulate ATF3 at the young basal mean.

Inputs: CORA_LOOM and CORA_ADJ. Outputs: CORA_OUT (default: out).
BC2/3/4/6 define basal stem cells. The model retains 3,000 HVGs plus
ATF3 and the explicitly listed additional genes.
Historical fitted outputs are preserved separately in saved_results/."""

import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
import scipy.sparse

OUT = os.environ.get("CORA_OUT", "out")
LOOM = os.environ.get("CORA_LOOM", "full.loom")
ADJ = os.environ.get("CORA_ADJ", "adjacencies.tsv")
os.makedirs(OUT, exist_ok=True)
ADDITIONAL_GENES = [
    "CDKN1A",
    "CDKN2A",
    "TP53",
    "IL6",
    "CXCL8",
    "IL1B",
    "IFNB1",
    "CXCL2",
    "ATF4",
    "DDIT3",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "MMP9",
    "MMP2",
]
BASAL_STEM = ["Basal cell-2", "Basal cell-3", "Basal cell-4", "Basal cell-6"]


def load_loom_adata(path):
    with h5py.File(path, "r") as f:
        genes = f["row_attrs/Gene"][:].astype(str)
        cells = f["col_attrs/CellID"][:].astype(str)
        mat = f["matrix"][:]
        meta = {}
        for key in ("celltype", "age_group", "site", "orig.ident"):
            if key in f["col_attrs"]:
                meta[key] = f["col_attrs"][key][:].astype(str)
    X = mat.T.astype(np.float32)
    adata = sc.AnnData(
        X=X, obs=pd.DataFrame(meta, index=cells), var=pd.DataFrame(index=genes)
    )
    adata.var_names_make_unique()
    return adata


def build_base_grn(adj_path, var_names):
    print(f"[co] loading adjacencies", flush=True)
    adj = pd.read_csv(adj_path, sep="\t")
    adj = adj[adj["target"].isin(var_names) & adj["TF"].isin(var_names)].copy()
    adj = (
        adj.sort_values("importance", ascending=False)
        .groupby("TF", group_keys=False)
        .head(20)
        .reset_index(drop=True)
    )
    tf_scores = adj.groupby("TF")["importance"].sum().sort_values(ascending=False)
    keep_tfs = tf_scores.head(300).index.tolist()
    if "ATF3" not in keep_tfs:
        keep_tfs.append("ATF3")
    adj = adj[adj["TF"].isin(keep_tfs)].copy()
    print(f"[co] base GRN: {len(adj):,} edges, {adj['TF'].nunique()} TFs", flush=True)
    base_grn = adj.pivot_table(
        index="target", columns="TF", values="importance", aggfunc="max", fill_value=0.0
    ).reset_index()
    base_grn.insert(1, "peak_id", base_grn["target"])
    base_grn.rename(columns={"target": "gene_short_name"}, inplace=True)
    return base_grn


def extract_sim_old(oracle, old_m):
    if "simulated_count" not in oracle.adata.layers:
        raise RuntimeError("CellOracle did not produce simulated_count")
    sim_X = oracle.adata.layers["simulated_count"]
    if scipy.sparse.issparse(sim_X):
        sim_X = sim_X.toarray()
    return np.asarray(sim_X, dtype=np.float64)[old_m].mean(axis=0)


def main():
    import celloracle as co

    adata = load_loom_adata(LOOM)
    print(f"[co] {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)
    tmp = adata.copy()
    sc.pp.normalize_total(tmp, target_sum=10000.0)
    sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(tmp, n_top_genes=3000, flavor="seurat")
    hv = set(tmp.var_names[tmp.var["highly_variable"].values])
    hv.add("ATF3")
    for g in ADDITIONAL_GENES:
        if g in adata.var_names:
            hv.add(g)
    adata = adata[:, [g for g in adata.var_names if g in hv]].copy()
    print(f"[co] kept {adata.n_vars} genes", flush=True)
    emb = adata.copy()
    sc.pp.normalize_total(emb, target_sum=10000.0)
    sc.pp.log1p(emb)
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=50)
    sc.pp.neighbors(emb, n_neighbors=15, n_pcs=30)
    sc.tl.umap(emb)
    adata.obsm["X_pca"] = emb.obsm["X_pca"]
    adata.obsm["X_umap"] = emb.obsm["X_umap"]
    adata.obs["cluster"] = adata.obs["celltype"]
    base_grn = build_base_grn(ADJ, list(adata.var_names))
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata, cluster_column_name="cluster", embedding_name="X_umap"
    )
    oracle.import_TF_data(TF_info_matrix=base_grn)
    print("[co] perform_PCA", flush=True)
    oracle.perform_PCA()
    print("[co] knn_imputation", flush=True)
    oracle.knn_imputation(
        n_pca_dims=30, k=24, balanced=True, b_sight=100, b_maxl=15, n_jobs=8
    )
    print("[co] get_links", flush=True)
    links = oracle.get_links(cluster_name_for_GRN_unit="cluster", alpha=10, n_jobs=16)
    links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    print("[co] fit_GRN_for_simulation", flush=True)
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
    var_names = list(oracle.adata.var_names)
    obs_df = oracle.adata.obs
    stem_mask = obs_df["celltype"].isin(BASAL_STEM).values
    young_m = stem_mask & (obs_df["age_group"].values == "Y")
    old_m = stem_mask & (obs_df["age_group"].values == "O")
    orig_X = oracle.adata.X
    if scipy.sparse.issparse(orig_X):
        orig_X = orig_X.toarray()
    orig_X = np.asarray(orig_X, dtype=np.float64)
    young_mean_expr = orig_X[young_m].mean(axis=0)
    old_mean_expr = orig_X[old_m].mean(axis=0)
    oX = oracle.adata[:, "ATF3"].X
    oX = (
        oX.toarray().flatten()
        if scipy.sparse.issparse(oX)
        else np.asarray(oX).flatten()
    )
    atf3_max = float(oX.max())
    atf3_young = float(oX[young_m].mean())
    atf3_old = float(oX[old_m].mean())
    print(
        f"[co] ATF3 imputed: max={atf3_max:.4f} young_stem_mean={atf3_young:.4f} old_stem_mean={atf3_old:.4f}",
        flush=True,
    )
    doses = {"young_mean": atf3_young}
    pd.DataFrame([{"dose": k, "atf3_target": v} for k, v in doses.items()]).to_csv(
        os.path.join(OUT, "dose_levels.csv"), index=False
    )
    print("[co] dose targets:", {k: round(v, 4) for k, v in doses.items()}, flush=True)
    res = pd.DataFrame(
        {"gene": var_names, "young": young_mean_expr, "old": old_mean_expr}
    )
    csv_path = os.path.join(OUT, "full_genome_dose_response.csv")
    for name, val in doses.items():
        print(f"[co] simulate_shift dose={name} ATF3={val:.4f}", flush=True)
        try:
            oracle.simulate_shift(
                perturb_condition={"ATF3": val},
                n_propagation=3,
                clip_delta_X=True,
                ignore_warning=name != "young_mean",
            )
            res[f"sim_{name}"] = extract_sim_old(oracle, old_m)
        except Exception as e:
            print(f"   !! dose {name} FAILED: {e}", flush=True)
            res[f"sim_{name}"] = np.nan
            res.to_csv(csv_path, index=False)
            raise RuntimeError(f"Simulation failed for {name}; checkpoint saved") from e
        gb = np.abs(old_mean_expr - young_mean_expr)
        ga = np.abs(res[f"sim_{name}"].values - young_mean_expr)
        rec = gb - ga
        flt = (
            (gb >= 0.05)
            & (np.maximum(young_mean_expr, old_mean_expr) >= 0.1)
            & (np.array(var_names) != "ATF3")
        )
        r = rec[flt]
        nr = (r > 0.005).sum()
        nd = (r < -0.005).sum()
        nf = len(r) - nr - nd
        print(
            f"   -> n={len(r)} restored={nr} ({100 * nr / len(r):.0f}%) disrupted={nd} ({100 * nd / len(r):.0f}%) flat={nf} ({100 * nf / len(r):.0f}%)",
            flush=True,
        )
        res.to_csv(csv_path, index=False)
        print(f"[co] checkpoint saved ({name})", flush=True)
    res.to_csv(csv_path, index=False)
    print(f"[co] wrote full_genome_dose_response.csv ({len(res)} genes)", flush=True)
    print("[co] DONE", flush=True)


if __name__ == "__main__":
    main()
