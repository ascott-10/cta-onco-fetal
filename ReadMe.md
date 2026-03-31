# Single-Cell RNA-Seq Processing and Annotaation Pipeline

A modular pipeline for processing and annotating scRNA-seq datasets across multiple projects.

Supports:

* Flexible input formats (10x, prebuilt .h5ad, custom matrices)
* Multi-sample projects
* Re-runnable annotation and visualization
* Persistent color mapping across datasets

---

## Project Structure


        
    /data/scottaa/cta_onco_fetal/
    ├── reference_data_common/
    │   ├── markers/
    │   ├── gene_maps/
    │   └── annotations/
    │
    └── datasets/
        └── [project_name]/
            ├── original_data/        # Raw source files (optional)
            ├── raw_data/             # Intermediate .h5ad (init, filtered, qc)
            ├── working_adata/        # Annotated AnnData objects
            ├── results/
            │   ├── figures/
            │   │   ├── qc_plots/
            │   │   └── umap/
            │   └── tables/
            └── colors.json           # Persistent color mappings
---

## Usage

The pipeline is executed through `main.py` and supports three modes:

| Mode    | Command | Examples |
|---------|--------| --------| 
| Process | python python_common_scripts/main.py process <project_name> | python python_common_scripts/main.py process fetal_gonad
| Analyze | python python_common_scripts/main.py refine <project_name>| python python_common_scripts/main.py refine fetal_gonad


---

## Pipeline Overview



### 1. Process Mode (`process`)

Handles raw data → QC .h5ad

Per sample:

1. Import

- Loads raw data depending on processing_type
- Supported:
   - prebuilt_h5
   - 10x-style matrices
    - custom formats

2. Filtering (`filter_data`)

- Removes low-quality cells
- Preserves `important_genes`

3. QC + Normalization (`run_qc`)

- Normalization and log transform
- Stores full data in .raw
- Saves QC metrics

Outputs:

    raw_data/
    ├── *_init.h5ad
    ├── *_before_qc.h5ad
    └── *_after_qc.h5ad

### 2. Refine Mode (`refine`)

Runs annotation + visualization.
Can be rerun without reprocessing raw data.

Per sample:

1. Load QC data

    `*_after_qc.h5ad`

2. Cell Annotation (`cell_annotation_wrapper`)
- Marker-based scoring (ULM)
- Leiden → cell type mapping
- Metadata annotation
- Gene binary features
3. Color Assignment (`set_obs_colors`)
- Uses persistent color maps
- Adds new entries only (no overwriting)
4. UMAP Visualization (`make_umaps`)
- Cell type plots
- Gene expression overlays
- Per-sample output

Outputs:

    working_adata/
    └── *_annotated_cellmarker.h5ad

    results/figures/umap/<sample_id>/
---

## Configuration (config.py)

All behavior is controlled centrally.

Base Paths

    BASE_DIR = "/data/scottaa/cta_onco_fetal"
    REFERENCE_DIR = BASE_DIR/reference_data_common

Global Config

    GLOBAL_CONFIG = {
        "important_genes": [...],
        "umap_obs_cols": ["leiden", "predicted_cell_type"],
    }
Color Defaults

    GENE_COLOR_MAP = {...}
    CELL_TYPE_COLORS = {"None": "#d3d3d3"}

These are:

- initial defaults only
- extended dynamically during runtime
- saved to colors.json
---

### Project Configuration

Each project defines:

    PROJECT_CONFIG = {
        "project_name": {
            "processing_type": "...",
            "sample_list": [...],              # or metadata-driven
            "sample_meta_filename": "...",     # optional
            "sample_meta_cols": [...],
            "cell_markers_file_path": ...,
            "cell_annotation_file_path": ...
        }
    }
Processing Types
| Type	|Description |
|-------|------|
|prebuilt_h5	|Already processed AnnData|
|10x|	Matrix + genes + cell metadata|


### Color Persistence (Important)

Color mappings are stored per project:

    datasets/<project_name>/colors.json

Behavior:

- Existing colors are reused
- New cell types / genes get assigned once
- Colors remain consistent across:
    - samples
    - reruns
    - annotation updates

---


## Adding a New Dataset

1. Add entry to PROJECT_CONFIG
2. Define:
- processing_type
- metadata (if needed)
- marker + annotation files
3. Place data in:

    datasets/<project_name>/original_data/

Run:

`python main.py process <project_name>`

`python main.py refine <project_name>`
