# Module I — scNexus-Process

Process a single scRNA-seq sample from a raw Cell Ranger matrix to an annotated
Seurat object, one pipeline step per tab.

**Pipeline:** Import → QC Filtering → PCA → Doublet Removal → Cell Cycle Scoring → Annotation

| | |
|---|---|
| **Input** | `datasets/Demo_Module_I_21401X3.rds` (demo), or a Cell Ranger MEX folder |
| **Output** | Processed Seurat `.rds`, downloadable at every step |
| **Next module** | [Module II](module-II.md) (integrate several of these objects) |

> Every tab is available at all times. If you jump ahead before finishing the
> previous step, the app warns with a "proceed anyway" dialog rather than
> blocking you — for the happy path, just go left to right.

---

## Step 1 — Import

Load a sample into a Seurat object and compute per-cell QC metrics. There are
two ways to load — you only need one.

**Option 2 — Processed Seurat object** (the demo path)

1. Open the **Import** tab.
2. Under *Option 2*, upload `datasets/Demo_Module_I_21401X3.rds`. Any missing QC
   metrics are filled in on load.

**Option 1 — Cell Ranger MEX folder** (starting from raw data)

1. Open the **Import** tab.
2. Under *Option 1*, click **Select** and browse to a Cell Ranger MEX folder —
   the one containing `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`.
3. Optionally adjust the sample ID and the cell/feature thresholds.
4. Click **Convert**. The app builds the Seurat object and computes
   mitochondrial, ribosomal, and hemoglobin percentages automatically.

> 📷 **Screenshot:** _Import tab showing the two load options_ — `img/module-I/01-import.png`

**Result:** a summary of the loaded object (cell and feature counts) and QC
metric columns ready for the next step.

---

## Step 2 — QC Filtering

Remove low-quality cells with adjustable thresholds and a live pass/fail preview.

1. Open the **QC** tab.
2. Adjust the sliders for feature count, UMI count, and mitochondrial % — the
   plots recolor cells as pass/fail in real time.
3. Click the filter action to apply the thresholds.

> 📷 **Screenshot:** _QC tab showing threshold sliders and pass/fail scatter_ — `img/module-I/02-qc.png`

**Result:** the object is subset to passing cells; the retained/removed counts
are reported.

---

## Step 3 — PCA

Normalize the data and run PCA.

1. Open the **PCA** tab.
2. Choose a normalization method — **LogNormalize** (default) or SCTransform.
3. Run PCA and inspect the elbow, loading, and heatmap plots to judge how many
   principal components to carry forward.

> 📷 **Screenshot:** _PCA tab with elbow plot after running PCA_ — `img/module-I/03-pca.png`

**Result:** a dimensionally-reduced object with PCA embeddings.

---

## Step 4 — Doublet Removal

Detect and remove likely doublets with DoubletFinder.

1. Open the **Doublet** tab.
2. Run doublet detection.
3. Review the UMAP with singlets/doublets highlighted, then apply removal.

> 📷 **Screenshot:** _Doublet tab UMAP colored by singlet/doublet_ — `img/module-I/04-doublet.png`

**Result:** predicted doublets removed from the object.

---

## Step 5 — Cell Cycle Scoring

Score each cell for cell-cycle phase (mouse gene sets).

1. Open the **Cell Cycle** tab.
2. Click **Run Cell Cycle** and view the S and G2M phase scores on the UMAP.

> 📷 **Screenshot:** _Cell Cycle tab UMAP colored by S and G2M phase scores_ — `img/module-I/05-cellcycle.png`

**Result:** `S.Score`, `G2M.Score`, and `Phase` added to the object.

---

## Step 6 — Annotation

Assign cell-type labels via SingleR (automatic) and/or manual annotation.

1. Open the **Annotation** tab.
2. Run SingleR auto-annotation, then refine with the manual annotation controls
   as needed.

> 📷 **Screenshot:** _Annotation tab with SingleR labels on the UMAP_ — `img/module-I/06-annotation.png`

**Result:** cell-type labels in the object metadata.

---

## Saving your object

Each step exposes a **Save** panel to download the current Seurat object as
`.rds`. Save after annotation to produce the input for [Module II](module-II.md).

> 📷 **Screenshot:** _Save panel with export options_ — `img/module-I/07-save.png`
