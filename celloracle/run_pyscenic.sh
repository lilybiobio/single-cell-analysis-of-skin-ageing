#!/usr/bin/env bash
# pySCENIC pipeline (GRN + ctx + AUCell).
# Set WORKDIR to the exported loom directory and RESOURCES to reference files.
# Usage:  ./run_pyscenic.sh <loom_name_without_ext>   [e.g. kera or full]
set -euo pipefail

LOOM_NAME="${1:-full}"
WORKDIR="${WORKDIR:?Set WORKDIR to the exported loom directory}"
RESOURCES="${RESOURCES:-$WORKDIR/resources}"
N_WORKERS="${N_WORKERS:-8}"

LOOM="$WORKDIR/${LOOM_NAME}.loom"
TFS="$RESOURCES/hs_hgnc_tfs.txt"
RANK="$RESOURCES/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
MOTIF="$RESOURCES/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

for f in "$LOOM" "$TFS" "$RANK" "$MOTIF"; do
  [ -f "$f" ] || { echo "[run] MISSING: $f" >&2; exit 2; }
done

OUTDIR="$WORKDIR/out_${LOOM_NAME}"
mkdir -p "$OUTDIR"
ADJ="$OUTDIR/adjacencies.tsv"
REG="$OUTDIR/regulons.csv"
AUC="$OUTDIR/aucell.loom"

echo "[run] loom=$LOOM  workers=$N_WORKERS  outdir=$OUTDIR"
echo "[run] stage 1: pyscenic grn"
time pyscenic grn \
  "$LOOM" "$TFS" \
  -o "$ADJ" \
  --num_workers "$N_WORKERS" \
  --method grnboost2 \
  --seed 777 \
  2>&1 | tee "$OUTDIR/grn.log"

echo "[run] stage 2: pyscenic ctx (motif pruning)"
time pyscenic ctx \
  "$ADJ" "$RANK" \
  --annotations_fname "$MOTIF" \
  --expression_mtx_fname "$LOOM" \
  --mode "dask_multiprocessing" \
  --num_workers "$N_WORKERS" \
  --output "$REG" \
  --mask_dropouts \
  2>&1 | tee "$OUTDIR/ctx.log"

echo "[run] stage 3: pyscenic aucell"
time pyscenic aucell \
  "$LOOM" "$REG" \
  --output "$AUC" \
  --num_workers "$N_WORKERS" \
  2>&1 | tee "$OUTDIR/aucell.log"

echo "[run] DONE — outputs in $OUTDIR"
ls -lh "$OUTDIR"
