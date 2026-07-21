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


**Option 1 — Cell Ranger MEX folder** (starting from raw data)

1. Open the **Import** tab.
2. Under *Option 1*, click **Select** and browse to a Cell Ranger MEX folder —
   containing `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`.
3. The created sample ID will be saved in meta.data under **"orig.ident"** as converting to a Seurat Object.
4. Optionally, adjust the sample ID and the cell threshold (minimum number of cells expressing a gene for the gene to be retained) feature threshold (minimum number of detected genes for the cell to be retained).
5. Click **Convert**. The app builds the Seurat object and computes
   mitochondrial, ribosomal, and hemoglobin percentages automatically.

**Option 2 — Processed Seurat object** (starting from processed data)

1. Open the **Import** tab.
2. Under *Option 2*, click **Browse...** to select a **Seurat Object (.rds)** file from your local computer. Any missing QC metrics are filled in on load.
3. For testing, upload `Demo_Module_I_21401X3.rds`. 

> 📷 **Screenshot:** _Import tab showing the two load options_ — `img/module-I/01-import.png`

**Result:** a summary of the loaded/converted Seurat object (total cell counts) and its QC metric (nFearure_RNA, nCount_RNA, percent.mt, percent.rp, and percent.hb) ready for the next step.

---

## Step 2 — QC Filtering

Remove low-quality cells with adjustable thresholds and a live pass/fail preview.

1. Open the **QC Filtering** tab.
2. Adjust the sliders for **nfeature(nFeature_RNA) Range**: nFeature_RNA is the number of unique genes detecked in each cell. **Low nFeature_RNA (<200)** indicates empty droplets or poor quality cells. Very high nFeature_RNA may imply possible doublets, multiplets, or highly trancriptionally active cell. No upper limit recommended. 
3. Adjust the sliders for **nCount(nCount_RNA) Ranger**: nCount_RNA is the total number of RNA modulecules detected in each cell. **Low nCount_RNA (<500)** indicates low-quality cell. Very high nCount_RNA may imply potential doublets, multiplets, or highly trancriptionally active cell. No upper limit recommended. 
5. Adjust the **Mitochondrial% (percent.mt)** below **10-20%** to remove damaged or dying cells qne low quality droplets.
6. **Ribosomal% (percent.rp)** adjustment is optional. Percentage of toal RNA counts derived from ribosomal protein genes. Use percent.rp with other QC matrics for low quality cell detection.
7. **Hemoglobin% (percent.hb)** adjustment is optional. Percentage of total RNA counts derived from hemoglobin genes can help detect erythroid cells or blooad contaimation. Interpret by cell type and study design.
8. Click **Plot** to visuliz **QC Matrix 1** and **QC Matrix 2** in **Plot Type** (**Scatter,Violin,Density** plots); Click **Update Selection** to re-color cells as pass/fail in real time.
9. Click **Filter Low Quality Cells** to apply the filting thresholds for next step.

> 📷 **Screenshot:** _QC tab showing threshold sliders and pass/fail scatter_ — `img/module-I/02-qc.png`

**Result:** the Seurat object is subset to cells passed QC; the retained cell number is updated in **QC Plot**.
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
