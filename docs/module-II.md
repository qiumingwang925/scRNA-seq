# Module II — scNexus-Integrate

Merge several processed samples, correct batch effects, benchmark the result,
and annotate the integrated object.

**Pipeline:** Upload & Merge → Integration → Benchmarking → Annotation

| | |
|---|---|
| **Input** | Two or more processed Seurat `.rds` objects (Module I output) |
| **Output** | Annotated, integrated Seurat `.rds` |
| **Next module** | [Module III](module-III.md) / [Module IV](module-IV.md) |

> For the demo, upload all three `datasets/Demo_Module_II_*.rds` files —
> `21401X3_control`, `22713X2_disease`, and `24143X4_control`. Integration is
> only meaningful with more than one sample.

---

## Step 1 — Upload & Merge

Load the per-sample objects, label them, and combine into one.

1. Open the **Upload & Merge** tab.
2. Upload two or more Seurat `.rds` objects (Module I outputs), adding each one.
3. In the metadata table, assign each sample its identifiers — e.g. sample ID
   and condition/group — so batches and biology can be told apart later.
4. Click the merge action to combine them into a single object.

> 📷 **Screenshot:** _Upload & Merge tab with samples listed and metadata assigned_ — `img/module-II/01-upload-merge.png`

**Result:** one merged Seurat object carrying per-sample metadata, ready to
integrate.

---

## Step 2 — Integration

Normalize and remove batch effects between samples.

1. Open the **Integration** tab.
2. Choose a normalization method — **LogNorm** (default) or SCTransform. The
   available *Regress Factors* update to match the method.
3. Choose a batch-integration method: **CCA**, **RPCA**, **Harmony**, or
   **FastMNN**.
4. Run integration, then inspect the UMAP colored by sample — well-integrated
   data has samples overlapping rather than forming separate islands.

> 📷 **Screenshot:** _Integration tab with method selected and integrated UMAP_ — `img/module-II/02-integration.png`

**Result:** an integrated embedding where shared cell types from different
samples mix together.

---

## Step 3 — Benchmarking

Quantify how well the integration worked.

1. Open the **Benchmarking** tab.
2. Compute the integration-quality metrics (**ASW**, **LISI**).
3. Read the ranked comparison to confirm your method choice — if a different
   method scores better, go back to Step 2 and re-run with it.

> 📷 **Screenshot:** _Benchmarking tab with ranked metric comparison_ — `img/module-II/03-benchmark.png`

**Result:** a ranked table of integration quality backing your choice of method.

---

## Step 4 — Annotation

Cluster and label the integrated object.

1. Open the **Annotation** tab.
2. Explore the integrated UMAP; subset and cluster as needed.
3. Run differential expression to find cluster markers.
4. Assign cell-type labels using the manual annotation controls.
5. Save the annotated object as `.rds`.

> 📷 **Screenshot:** _Annotation tab with labeled clusters_ — `img/module-II/04-annotation.png`

**Result:** an annotated, integrated object — the input for
[Module III](module-III.md) and [Module IV](module-IV.md).
