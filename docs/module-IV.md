# Module IV — scNexus-Interact

Infer cell-cell communication with two parallel engines — **CellChat** and
**LIANA** — over a single shared upload. Run either or both.

**Tabs:** Upload → Run CellChat → Visualize CellChat → Run LIANA → Visualize LIANA

| | |
|---|---|
| **Input** | Annotated Seurat `.rds` (Module II output). For testing, use `Demo_Module_III_IV_(control_vs_disease).rds` |
| **Output** | CellChat or LIANA result object (`.rds`) and figures and tables |

> The single upload feeds all four downstream tabs. Every visualization panel is
> gated behind a **Generate Plot** button, so tweaking parameters queues up
> without triggering an expensive re-render until you ask for it.

---

## Step 1 — Upload

Load the object and choose the grouping used for communication.

1. Open the **Upload** tab.
2. Click **Brower...** to load an annotated `.rds` — the object [Module II](module-II.md) produced.
3. Pick the metadata column to use as **Condition / Group column** (`group`).
4. Click **Confirm & Prepare Object** and alidate the information in the *Validation* tex box.

> 📷 **Screenshot:** _Upload tab with grouping column selected_ — `img/module-IV/01-upload.png`

**Result:** a validated, grouped object shared with all four downstream tabs.

---

## Step 2 — Run CellChat

Compute communication with the CellChat engine.

1. Open the **Run CellChat** tab.
2. Run the pipeline — CellChat is computed per group on the Secreted Signaling
   database, then merged across groups for cross-group comparison.

> 📷 **Screenshot:** _CellChat compute tab after a run completes_ — `img/module-IV/02-cellchat-run.png`

**Result:** a CellChat result object, downloadable as `.rds`. (This step is
memory-hungry — give the app plenty of memory if it dies mid-run.)

---

## Step 3 — Visualize CellChat

Explore the CellChat result across four subtabs. Each panel renders on
**Generate Plot**.

1. **Global Network** — overall interaction counts / strengths between groups.
2. **Zoom-in** — a chosen signaling pathway or ligand-receptor pair.
3. **Signaling-Focused** — outgoing/incoming signaling roles per group.
4. **Communication Patterns** — including the manifold-learning embedding.

> 📷 **Screenshot:** _CellChat global network circle plot_ — `img/module-IV/03-cellchat-vis.png`

**Result:** downloadable figures for each selected view.

---

## Step 4 — Run LIANA

Compute communication with the LIANA engine.

1. Open the **Run LIANA** tab.
2. Run the per-group multi-method pipeline. Mouse ↔ human ortholog mapping is
   applied automatically when the chosen resource is not `MouseConsensus`.

> 📷 **Screenshot:** _LIANA compute tab after a run completes_ — `img/module-IV/04-liana-run.png`

**Result:** a LIANA result object, downloadable as `.rds`.

---

## Step 5 — Visualize LIANA

Explore the LIANA result across three subtabs. Each panel renders on
**Generate Plot**.

1. **CCC Dot Plot** — single group.
2. **CCC Freq Heatmap** — N-group grid.
3. **CCC Freq Chord Diagram** — N-group grid.

> 📷 **Screenshot:** _LIANA CCC dot plot_ — `img/module-IV/05-liana-vis.png`

**Result:** downloadable figures for each selected view.
