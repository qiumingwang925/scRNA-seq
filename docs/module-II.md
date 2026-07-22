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
4. Compute the integration-quality metrics (**ASW**, **LISI**, **GraphLISI**) for both *batch mixing* and *cell-type separation*.
- ASW
  - basw (batch): Worst = 1 or -1; Best = 0
  - casw (cell-type): Worst = -1; Best = 1
- LISI
  - ilisi score (batch mixing): worst = 1; best = number of batches
  - clisi score (cell-type separation): worst = number of cell types ; best = 1
- GraphLISI
  - gilisi score (batch mixing): worst = 1; best = number of batches
  - gclisi score (cell-type separation): worst = number of cell types ; best = 1
5.  Utilize **Rank Summary** and **Score Distribution** to compare method performance. *Median Score Summry* displays scores rounded to three decimal places, whereas *Rank Summary* calculates rankings using the full-precision scores.


> 📷 **Screenshot:** _Benchmarking tab with ranked metric comparison_ — `img/module-II/03-benchmark.png`

**Result:** *Rank Summary* and *Score Distribution* summarize the integration performance of each method to support method selection. Use these results together with the UMAP Visual Inspection in the next step to make the final decision.
 
---

## Step 4 — Annotation

Cluster and label the integrated object.

1. Open the **Annotation** tab.
2. Use **Select Integration** to explore the unintegrated or integrated UMAP.
3. Use **Color by Metadata** to visualize `orig.ident`, `batch`, `group`, or any cell-type `annotation`, and visually assess integration performance by evaluating
- batch mixing: same cell type from different samples overlap with each other
- cell-type separation: (sub-)cell types from different samples remain separated
5. Explore the data using the available tools, including interactive cell selection, re-clustering, coloring metadata variables, gene expression visualization, differential expression analysis, and manual cell-type labeling (same workflow as the **Module I Annotation** tab)
6. Save the annotated object as `.rds`
- Full object: multiple integrated PCA

> 📷 **Screenshot:** _Annotation tab with labeled clusters_ — `img/module-II/04-annotation.png`

**Result:** an annotated, integrated object — the input for
[Module III](module-III.md) and [Module IV](module-IV.md).
