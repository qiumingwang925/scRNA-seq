# Module III — Post-Annotation Exploration

Explore and visualize an annotated object. Unlike the earlier modules this is
not a linear pipeline — once data is uploaded, all visualization tabs are
available in any order.

**Tabs:** Upload → UMAP Cell-Type → Violin → Dot Plot → Heatmap → Differential Expression → Enrichment

| | |
|---|---|
| **Input** | Annotated Seurat `.rds` (Module II output; slimmed `data`-only export is fine) |
| **Output** | Figures and tables (downloadable per tab) |

---

## Step 1 — Upload

1. Open the **Upload** tab and load an annotated `.rds`
   (e.g. `test-data/Manual_Annotated_21401X3.rds`).
2. Confirm the structure summary looks correct.

> 📷 **Screenshot:** _Upload tab with object summary_ — `img/module-III/01-upload.png`

---

## Step 2 — UMAP Cell-Type

Highlight cell types, view gene (co-)expression, and lasso-select for re-UMAP.

> 📷 **Screenshot:** _UMAP colored by cell type_ — `img/module-III/02-umap.png`

---

## Step 3 — Violin Plot

Gene expression by cell type, with stacking and `split.by` options.

> 📷 **Screenshot:** _Violin plot for selected genes_ — `img/module-III/03-violin.png`

---

## Step 4 — Dot Plot

Gene expression dot plot across identities.

> 📷 **Screenshot:** _Dot plot for a marker gene set_ — `img/module-III/04-dot.png`

---

## Step 5 — Heatmap

Scaled-expression heatmap using HVGs or a custom gene list.

> 📷 **Screenshot:** _Heatmap of scaled expression_ — `img/module-III/05-heatmap.png`

---

## Step 6 — Differential Expression

Compare configurable cell populations, with metadata subsetting.

> 📷 **Screenshot:** _DE results table_ — `img/module-III/06-de.png`

---

## Step 7 — Enrichment

Run enrichR pathway analysis on the DE results and view the bar plot.

> 📷 **Screenshot:** _Enrichment bar plot_ — `img/module-III/07-enrichment.png`
