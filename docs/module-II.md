# Module II — scNexus-Integrate

Merge several processed samples, correct batch effects, benchmark the result,
and annotate the integrated object.

**Pipeline:** Upload & Merge → Integration → Benchmarking → Annotation

| | |
|---|---|
| **Input** | Two or more processed Seurat `.rds` objects (Module I output). Testing datasets`Demo_Module_II_21401X3_control.rds`, 'Demo_Module_II_24143X4_control.rds`, `Demo_Module_II_22713X2_disease.rds` |
| **Output** | Annotated, integrated Seurat `.rds` |
| **Next module** | [Module III](module-III.md) / [Module IV](module-IV.md) |

> For the demo, upload all three `datasets/Demo_Module_II_*.rds` files —
> `21401X3_control`, `22713X2_disease`, and `24143X4_control`. Integration is
> only meaningful with more than one sample.

---

## Step 1 — Upload & Merge

Load the per-sample objects, label them, and merge into one.

1. Open the **Upload & Merge** tab.
2. Click **Browse...**. to select two or more Seurat objects `.rds` files (Module I outputs) and upload them together.
3. In the **Metadata Assignemen** table, assign each sample its batch and group.
4. Click the **Merge Objects** to combine them into a single object.
5. Once merged, a summary will be showed under **Merge Status**

> 📷 **Screenshot:** _Upload & Merge tab with samples listed and metadata assigned_ — `img/module-II/01-upload-merge.png`

**Result:** one merged Seurat object carrying per-sample metadata, ready to
integrate.

---

## Step 2 — Integration

Normalize and transform data. Apply one or more integration methods to remove batch effects between samples.

1. Open the **Integration** tab.
2. Choose a normalization method — **Log Normalization** (default) or **SCTransform**. The
   available *Regress Factors* update to match the method.
   Note: For Log Normalization, *Regress Factor* **nCount_RNA** is recommended.
4. Choose a batch-integration method: **CCA**, **RPCA**, **Harmony**, or
   **FastMNN** (only available for Log Normalization).
5. Click **Run Pipeline** to run the integration workflows.
6. **Pipeline Status Log** summrizes the normalization and reductions of process integration methods.

> 📷 **Screenshot:** _Integration tab with method selected and integrated UMAP_ — `img/module-II/02-integration.png`

**Result:** multiple integrated PCA and UMAP embeddings added to the Seurat object.

---

## Step 3 — Benchmarking

Quantify how well the integration worked.

1. Open the **Benchmarking** tab.
2. Select **Batch Label** (batch information assigned in **Step 1**).
3. Select **Cell-Type Label** (manual or SingleR annotation).
4. Compute the integration-quality metrics (**ASW**, **LISI**, **GraphLISI**).
- ASW:
  -basw
  -casw
- LISI:
- GraphLISI:
6. Read the ranked comparison to confirm your method choice


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
