# Module II — Multi-Sample Integration

Merge several processed samples, correct batch effects, benchmark the result,
and annotate the integrated object.

**Pipeline:** Upload & Merge → Integration → Benchmarking → Annotation

| | |
|---|---|
| **Input** | Two or more processed Seurat `.rds` objects (Module I output) |
| **Output** | Annotated, integrated Seurat `.rds` |
| **Next module** | [Module III](module-III.md) / [Module IV](module-IV.md) |

> For a quick demo without running Module I on multiple samples, you can upload
> `test-data/Manual_Annotated_21401X3.rds` to see the tabs populated, though
> integration is most meaningful with several samples.

---

## Step 1 — Upload & Merge

1. Open the **Upload & Merge** tab.
2. Upload the Seurat `.rds` objects, assign metadata (e.g. sample / condition).
3. Merge into a single object.

> 📷 **Screenshot:** _Upload & Merge tab with samples listed and metadata assigned_ — `img/module-II/01-upload-merge.png`

---

## Step 2 — Integration

1. Open the **Integration** tab.
2. Choose a normalization method (LogNormalize or SCTransform).
3. Choose a batch-integration method (CCA, RPCA, Harmony, or FastMNN) and run.

> 📷 **Screenshot:** _Integration tab with method selected and integrated UMAP_ — `img/module-II/02-integration.png`

---

## Step 3 — Benchmarking

1. Open the **Benchmarking** tab.
2. Compute integration-quality metrics (ASW, LISI) to compare methods.

> 📷 **Screenshot:** _Benchmarking tab with ranked metric comparison_ — `img/module-II/03-benchmark.png`

---

## Step 4 — Annotation

1. Open the **Annotation** tab.
2. Explore the integrated UMAP, cluster, run DE, and assign cell-type labels.
3. Save the annotated object as `.rds` for downstream modules.

> 📷 **Screenshot:** _Annotation tab with labeled clusters_ — `img/module-II/04-annotation.png`
