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

**Result:** a summary of the loaded or converted Seurat object is displayed, including the total number of cells and the QC metrics saved in the object's 'meta.data': 'nFearure_RNA', 'nCount_RNA', 'percent.mt', 'percent.rp', and 'percent.hb'. The object is then ready for the next step.

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

---

## Step 3 — PCA

Normalize the data and run PCA.

1. Open the **PCA** tab.
2. Choose a normalization/tranformation method
- **Log Normalization** (default): Normalizes gene expression counts by dividing each cell's counts by its total counts, multiplying by a scale factor (10,000 by default), and applying a log(x + 1) transformation.
- **SCTransform**: Normalizes gene expression counts using a regularized negative binomial regression model that corrects for sequencing depth and stabilizes variance. The resulting Pearson residuals are used as normalized expression values.
3. **Run PCA** and inspect the elbow, loading, and heatmap plots by clicking **Plot** and **Update Plot** to determine the number of principal components (PCs) to retain for downstream analyses.
- **Elbow Plot**: Shows the variance explained by each PC. Select the number of PCs near the "elbow", where additional PCs contribute only small increases in explained variance.
- **Loading Plot**: Displays the genes that contribute most strongly to each PC. Check whether the top genes represent meaningful biological signals (e.g., identification biomarkers) rather than technical artifacts (e.g., mitochondrial, ribosomal, or cell-cycle genes, if not expected).
- **Heatmap Plot**: Visualizes the expression patterns of the top genes across cells for each PC. Retain PCs that separate biologically distinct cell populations and show interpretable gene expression patterns.
**Tips**:For most datasets, retaining approximately 20–30 PCs provides a reasonable starting point. The optimal number should be guided by the elbow plot and the biological relevance of the loading and heatmap plots.
> 📷 **Screenshot:** _PCA tab with elbow plot after running PCA_ — `img/module-I/03-pca.png`

**Result:** Seurat Object adds a dimensionally-reduced PCA embeddings.

---

## Step 4 — Doublet Removal

Detect and remove likely doublets with **DoubletFinder**.

1. Open the **Doublet Removal** tab.
2. Select **Assay Settings** used to generate your dataset. This setting determines the default expected doublet rate used by DoubletFinder when estimating the expected number of doublets.
- 10x Genomics **High Throughput v3.1**: Default expected doublet rate = **4%**
- 10x Genomics **Standard v3.1**: Default expected doublet rate = **8%**
Note: These values are recommended starting points based on the assay chemistry. If the expected doublet rate for your experiment is known (e.g., from the number of recovered cells or the sequencing provider's recommendations), you can manually adjust the value.
3. **Calculate Default parameters** including pK (Optimal) and pN(Default).
4. **Run doublet detection**
5. Review the UMAP and scatter plot (nCounts_RNA vs nFeature_RNA) with singlets/doublets highlighted, then apply **Removal Doublets**.

> 📷 **Screenshot:** _Doublet tab UMAP colored by singlet/doublet_ — `img/module-I/04-doublet.png`

**Result:** Predicted doublets removed from the Seurat object.

---

## Step 5 — Cell Cycle Scoring

Score each cell for cell-cycle phase (mouse gene sets).

1. Open the **Cell Cycle** tab.
2. Click **Run Cell Cycle** and view the S and G2M phase scores on the UMAP.
3. Optionally, **Download Seurat Object** at this stage, as cell annotation may require additional time or multiple rounds of refinement.

> 📷 **Screenshot:** _Cell Cycle tab UMAP colored by S and G2M phase scores_ — `img/module-I/05-cellcycle.png`

**Result:** `S.Score`, `G2M.Score`, and `Phase` added to the meta.date of Seurat object.

---

## Step 6 — Annotation

Assign cell-type labels via SingleR (automatic) and/or manual annotation.

1. Open the **Annotation** tab.
2. Optionally, **Upload Seurat Object** if the dataset wasn't processed from beginning of the workflow.
3. The **Plotly** toolbar on the upper-right coner of the UMAP provides interactive tools for exploring and selecting cells:
- **Download a plot as png**: save the current plot as a PNG image.
- **Zoom**: zoom into a selected region of the plot.
- **Pan**: move the plot while maintaining the current zoom level.
- **Box Select**: select cells within a rectangular region.
- **Lasso Select**: select cells by drawing a freehand boundary.
- **Zoom out**: decrease the magnification of the plot.
- **Zoom in**: increase the magnification of the plot.
- **Autoscale**: automatically adjust the axes to fit the displayed data.
- **Reset axes**: restore the original plot view.
- **Show closest data on hover**: display information for the data point nearest the cursor.
- **Compare data on hover**: display information for all nearby data points at the cursor location.
5. **Manual Annotation**
- **Display UMAP from** shows the all cells and selected cells (after subset and re-clusring).
- **PCs for Subset UMAP** defines the number of PCs used to **Run UMAP on Selection** and re-cluering.
- **Metadata**: use **Color by Metadata** to visualize variables in metadata, adjust **Cluster Resolution** and "**Run Clustering** on all or seleted cells.
- **Expression**: to visualize gene expression in UMAP, the gene name can be selected from **Upload Biomarker CSV** file **Or Type Gene Name** directly.
- **Differential Expression**: Click **Run DE Analysis** to identfy top genes/markers for the current selected cells vs the rest. View results under sub-tab DE AnAlysis
- **Manual Annotation**: **Enter New Labels** to mark the current selected cells. Then click **Apply Label to Selection** to save the labels in `manual_annotation` of metadate.
6.**SingleR Annotation**
- **Upload Reference SCE (.rds)**: upload a SingleCellExperiment object saved as rds file as a reference.
- **Label Name**: the name (e.g. SingleR_Reference) will be saved in Seurat Object metadata. 
- Click **Run SingleR** to perform automated cell type annotation. Once finished, an **Annotation Summary** will be generated. The predicted cell type annotations can be visualizated under **Manual Annotation**.
7. **Export Annotated Object**: 
— **Full object**: raw counts, normalized/tranaformed data, scaled data, metadata, graphs, PCA embedding, UMAP embeddding
— **Module II Analysis-Ready Object**: raw counts and metadata
— **Customized Object**: selection based on needs
> 📷 **Screenshot:** _Annotation tab with SingleR labels on the UMAP_ — `img/module-I/06-annotation.png`

**Result:** `manual_annotation`and/or SingleR annotation labels add in the Seurat object metadata. Multiple versions of Seurat object can be downloaded as `.rds`.
