# Module IV — Cell-Cell Communication

Infer cell-cell communication with two parallel engines — **CellChat** and
**LIANA** — over a single shared upload. Run either or both.

**Tabs:** Upload → Run CellChat → Visualize CellChat → Run LIANA → Visualize LIANA

| | |
|---|---|
| **Input** | Annotated Seurat `.rds` (Module II output; Module III's slimmed export also works) |
| **Output** | CellChat or LIANA result object (`.rds`) |

---

## Step 1 — Upload

1. Open the **Upload** tab and load an annotated `.rds`
   (e.g. `test-data/Manual_Annotated_21401X3.rds`).
2. Pick the metadata column to use as the grouping (`group`) — usually cell type.

> 📷 **Screenshot:** _Upload tab with grouping column selected_ — `img/module-IV/01-upload.png`

---

## Step 2 — Run CellChat

Run the per-group CellChat pipeline (Secreted Signaling database), followed by a
merge for cross-group comparison.

> 📷 **Screenshot:** _CellChat compute tab after a run completes_ — `img/module-IV/02-cellchat-run.png`

---

## Step 3 — Visualize CellChat

Four subtabs: Global Network, Zoom-in (pathway / L-R pair), Signaling-Focused,
and Communication Patterns. Each panel renders on **Generate Plot**.

> 📷 **Screenshot:** _CellChat global network circle plot_ — `img/module-IV/03-cellchat-vis.png`

---

## Step 4 — Run LIANA

Run the per-group multi-method LIANA pipeline. Mouse ↔ human ortholog mapping is
applied automatically when the chosen resource requires it.

> 📷 **Screenshot:** _LIANA compute tab after a run completes_ — `img/module-IV/04-liana-run.png`

---

## Step 5 — Visualize LIANA

Three subtabs: CCC Dot Plot (single group), CCC Freq Heatmap, and CCC Freq Chord
Diagram. Each panel renders on **Generate Plot**.

> 📷 **Screenshot:** _LIANA CCC dot plot_ — `img/module-IV/05-liana-vis.png`
