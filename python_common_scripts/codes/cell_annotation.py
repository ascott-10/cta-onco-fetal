# python_common_scripts/codes/cell_annotation.py

import os
import json
import re
import numpy as np
import pandas as pd
import random

import scanpy as sc
import scipy.sparse as sp

import matplotlib
from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
import decoupler as dc
import colorcet as cc
import glasbey

sc.settings.autosave = True  # Saves plots to figdir automatically
sc.settings.autoshow = False # Prevents plots from displaying inline

from codes.differential_gene_expression import get_ranked_genes
from config import *

def get_markers_df(markers_file_path):

    """markers_df, cellmarkers_annotation_df = get_markers_df(markers_file_path)
    """
    if os.path.exists(markers_file_path):
        cellmarkers_annotation_df = pd.read_csv(markers_file_path)
        if not set(["source", "target"]).issubset(cellmarkers_annotation_df.columns):
            raise ValueError("markers_file must have columns: ['source','target']")
        markers_df = cellmarkers_annotation_df[["source", "target"]]
    
    else:
        # Use "PanglaoDB" to create a markers_df
        try:
            markers = dc.get_resource(name="PanglaoDB", organism="human", license="academic")
        except AttributeError:
            markers = dc.op.resource("PanglaoDB", organism="human", license="academic")
        markers = markers[
            markers["human"].astype(bool)
            & markers["canonical_marker"].astype(bool)
            & (markers["human_sensitivity"].astype(float) > 0.5)
        ].copy()
        markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
        markers_df = markers.rename(
            columns={"cell_type": "source", "genesymbol": "target"}
        )[["source", "target"]]
        markers_df.to_csv(os.path.join("reference_data_common", "markers", "PanglaoDB_markers.csv"))
        cellmarkers_annotation_df = markers_df.copy()
        markers_df.to_csv(os.path.join("reference_data_common", "markers", "PanglaoDB_markers_annotation.csv"))


    return markers_df, cellmarkers_annotation_df



def get_cell_type(sample_id, adata, markers_df, cell_type_label_col="predicted_cell_type"):
    """
    Performs cell type annotation using decoupler's ULM method.

    This function is corrected to run the analysis on the complete gene set stored
    in adata.raw, ensuring that cell type markers are not missed due to HVG filtering.
    adata = get_cell_type(sample_id, adata, markers_df, cell_type_label_col="predicted_cell_type")
    """
    
    # Use the var_names from the raw object for all gene-related operations.
    
    adata_var_names_upper = adata.raw.var_names.astype(str).str.upper()

    net = markers_df.drop_duplicates(subset=["source", "target"]).copy()
    net["target"] = net["target"].astype(str).str.upper()
    net = net.reset_index(drop=True)

    # Check for overlap against the FULL gene list from adata.raw.
    overlap = set(adata_var_names_upper) & set(net["target"])
    
    print(f"[{sample_id}] Overlapping genes: {len(overlap)}")

    if len(overlap) == 0:
        print(f"[{sample_id}] No overlap between adata.raw.var_names and marker targets. Skipping ULM.")
        adata.obs[cell_type_label_col] = adata.obs["leiden"].astype(str)
        return adata

    # Run ULM on a temporary AnnData object created from adata.raw.
    # This ensures decoupler sees all genes, not just the highly variable ones.
    print(f"[{sample_id}] Running ULM on full gene set ({adata.raw.n_vars} genes).")
    
    # Create a temporary AnnData object that contains the log-normalized data for all genes.
    adata_for_ulm = adata.raw.to_adata()
    
    # copy the Leiden cluster assignments from the main 'adata' object to our temporary object so that the ranking step knows which cells belong to which cluster.
    adata_for_ulm.obs['leiden'] = adata.obs['leiden'].copy()

    source_sizes = net.groupby("source")["target"].nunique()
    valid_sources = source_sizes[source_sizes >= 3].index
    net = net[net["source"].isin(valid_sources)].copy()

    print(f"[{sample_id}] Sources kept: {len(valid_sources)} / {len(source_sizes)}")  

    X = adata_for_ulm.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    noise = np.random.normal(0, 1e-6, size=X.shape)
    adata_for_ulm.X = X + noise
    # Run ULM on the temporary object that contains all genes.

    try:
        dc.mt.ulm(data=adata_for_ulm, net=net, tmin=0)
    except AssertionError as e:
        print(f"[{sample_id}] ULM failed: {e}. Using Leiden clusters as fallback.")
        adata.obs[cell_type_label_col] = adata.obs["leiden"].astype(str)
        return adata

    # Get ULM scores from the temporary object.
    score = dc.pp.get_obsm(adata_for_ulm, key="score_ulm")
    score_df = score.to_df() if hasattr(score, "to_df") else score

   # assign the results back to the original 'adata' object.

    # Per-cell predicted cell type, assigned back to the main adata object.
    adata.obs["ulm_score_cell_type"] = score_df.idxmax(axis=1)

    # Cluster-level ranking (using the scores and leiden clusters from the temp object).
    try:
        df = dc.tl.rankby_group(adata=score, groupby="leiden", reference="rest", method="t-test_overestim_var")
    except Exception as e:
        print(f"[{sample_id}] rankby_group failed: {e}")
        df = None

    if df is None or df.empty:
        print("writing predicted cell type from leiden to ", cell_type_label_col)
        if "leiden" in adata.obs:
            adata.obs[cell_type_label_col] = adata.obs["leiden"].astype(str)
        else:
            print(f"[{sample_id}] No Leiden clusters found. Assigning all cells to 'unknown'.")
            print("writing predicted cell type 'unknown' to ", cell_type_label_col)
            adata.obs[cell_type_label_col] = pd.Categorical(["unknown"] * adata.n_obs)
    else:
        df = df[df["stat"] > 0]
        dict_ann = (df.groupby("group").head(1).set_index("group")["name"].to_dict())
        s = adata.obs["leiden"].astype(str)
        adata.obs[cell_type_label_col] = pd.Categorical(s.map(dict_ann).fillna(s))
        print("writing predicted cell type from ulm-->leiden clusters to ", cell_type_label_col)

    adata.obs[cell_type_label_col] = adata.obs[cell_type_label_col].astype("category")
    
    return adata


def annotate_cells(adata, cellmarkers_annotation_df, cell_type_label_col="predicted_cell_type"):
    hierarchy_cols = ["target", "Lineage", "tissue_class", "tissue_type", "Cell_Category"]

    # keep only columns that actually exist
    available_cols = [col for col in hierarchy_cols if col in cellmarkers_annotation_df.columns]

    # must have at least the key column to proceed
    if "target" not in available_cols:
        return adata

    def most_common(x):
        return x.mode().iloc[0] if not x.mode().empty else None

    annotation_df = (
        cellmarkers_annotation_df[available_cols]
        .groupby("target", as_index=False)
        .agg(most_common)
    )

    annotation_df[cell_type_label_col] = annotation_df["target"]
    annotation_df = annotation_df.set_index(cell_type_label_col, drop=False)

    for col in annotation_df.columns:
        if col == cell_type_label_col:
            continue
        adata.obs[col] = adata.obs[cell_type_label_col].map(annotation_df[col])

    return adata

def add_gene_binary_columns(adata, genes, threshold=0.1):

    """adata = add_gene_binary_columns(adata, genes, threshold=0.1)"""

    if adata.raw is None:
        raise ValueError("adata.raw is None")

    for g in genes:

        if g not in adata.raw.var_names:
            continue

        expr = adata.raw[:, g].X

        if sp.issparse(expr):
            expr = expr.toarray().ravel()
        else:
            expr = np.asarray(expr).ravel()

        pos_mask = expr > threshold

        adata.obs[f"{g}_binary"] = pd.Categorical(
            np.where(pos_mask, f"{g}_pos", f"{g}_neg"),
            categories=[f"{g}_neg", f"{g}_pos"],
        )

    return adata

def set_obs_colors(adata, palette, cell_type_colors, gene_color_map):

    if isinstance(palette, str):
        palette = list(cc.palette[palette])

    for col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].astype("category")

        categories = sorted(adata.obs[col].cat.categories)
        assigned = []

        for i, cat in enumerate(categories):

            if col == "predicted_cell_type":
                if cat not in cell_type_colors:
                    cell_type_colors[cat] = palette[len(cell_type_colors) % len(palette)]
                assigned.append(cell_type_colors[cat])

            elif col.endswith("_binary"):
                gene = col.replace("_binary", "")

                if cat == f"{gene}_pos":
                    if gene not in gene_color_map:
                        gene_color_map[gene] = palette[len(gene_color_map) % len(palette)]
                    assigned.append(gene_color_map[gene])
                else:
                    assigned.append("#d3d3d3")

            else:
                assigned.append(palette[i % len(palette)])

        adata.uns[f"{col}_colors"] = assigned

    return adata, gene_color_map, cell_type_colors

def cell_annotation_wrapper(sample_id, adata_qc, adata_annotated_path, markers_file_path, gene_color_map, cell_type_colors, genes_list, cell_type_label_col = "predicted_cell_type", palette="glasbey"):
    adata = adata_qc.copy()
    
    markers_df, cellmarkers_annotation_df = get_markers_df(markers_file_path)
    adata = get_cell_type(sample_id, adata, markers_df, cell_type_label_col="predicted_cell_type")
    adata = annotate_cells(adata, cellmarkers_annotation_df, cell_type_label_col="predicted_cell_type")
    adata = add_gene_binary_columns(adata, genes_list, threshold=0.1)
    adata, gene_color_map, cell_type_colors = set_obs_colors(adata, palette, cell_type_colors, gene_color_map)

    adata.write(adata_annotated_path)

    return adata, gene_color_map, cell_type_colors

# New for embryos
def adata_split(adata_female_path, adata_male_path, adata_female_annotated_path,adata_male_annotated_path, adata_qc_path):
    # Split data into male and female

    if os.path.exists(adata_female_annotated_path) and os.path.exists(adata_male_annotated_path):
        adata_f = sc.read_h5ad(adata_female_annotated_path)
        adata_m = sc.read_h5ad(adata_male_annotated_path)
    
    else:

        if os.path.exists(adata_female_path) and os.path.exists(adata_male_path):
            adata_f = sc.read_h5ad(adata_female_path)
            adata_m = sc.read_h5ad(adata_male_path)
        else:

            adata_annotated = sc.read_h5ad(adata_qc_path)

            adata_f = adata_annotated[adata_annotated.obs["sex"] == "female"].copy()
            adata_m = adata_annotated[adata_annotated.obs["sex"] == "male"].copy()


            sc.pp.neighbors(adata_f, use_rep="X_pca", n_pcs=30)
            sc.tl.umap(adata_f)

            sc.pp.neighbors(adata_m, use_rep="X_pca", n_pcs=30)
            sc.tl.umap(adata_m)

            adata_f.write(adata_female_path)
            adata_m.write(adata_male_path)

    return adata_f, adata_m
        

def natural_key(x):

    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', str(x))]

def add_stage_label(adata):
    if "stage_label" not in adata.obs:
        if "developmental_stage" in adata.obs.columns:
            adata.obs["age"] = adata.obs["developmental_stage"].str.extract(r"(\d+)").astype(float)
            ordered_ages = sorted(adata.obs["age"].dropna().unique())

            adata.obs["stage_label"] = adata.obs["age"].astype("Int64").astype(str) + "W"
            ordered_labels = [f"{int(age)}W" for age in ordered_ages]

            adata.obs["stage_label"] = pd.Categorical(
                adata.obs["stage_label"],
                categories=ordered_labels,
                ordered=True
            )
    return adata

def build_palettes(adata, adata_annotated_path, cols, existing_palettes=None, palette="tab10"):

    import copy
    all_palettes = copy.deepcopy(existing_palettes) if existing_palettes else {}

    for col in cols:
        if col not in adata.obs:
            continue

        categories = sorted([str(c) for c in adata.obs[col].dropna().unique()], key=natural_key)
        adata.obs[col] = adata.obs[col].astype(str).astype("category").cat.set_categories(categories)

        if col == "sex":
            palette_dict = {"female": "#E78AC3", "male": "#4C72B0"}

        elif col == "stage_label":
            cmap = matplotlib.colormaps.get_cmap("Purples")
            n = len(categories)
            palette_list = [cmap(i / max(n - 1, 1))[:3] for i in range(n)]
            palette_dict = {cat: to_hex(color) for cat, color in zip(categories, palette_list)}

        elif col.startswith("leiden"):
            colors = glasbey.create_palette(palette_size=len(categories))
            palette_dict = {cat: to_hex(color) for cat, color in zip(categories, colors)}

        else:
            # unified persistent behavior (celltype-like)
            key = "celltype" if col.startswith("celltype_") else col

            if key not in all_palettes:
                all_palettes[key] = {}

            palette_dict = all_palettes[key]

            missing = sorted([c for c in categories if c not in palette_dict], key=natural_key)

            if missing:
                existing_n = len(palette_dict)
                full_palette = glasbey.create_palette(palette_size=existing_n + len(missing))
                new_colors = full_palette[existing_n:]
                new_colors = [to_hex(c) for c in new_colors]

                for k, c in zip(missing, new_colors):
                    palette_dict[k] = c

        if not col.startswith("celltype_") and not col.startswith("cellstate_"):
            all_palettes[col] = palette_dict

        adata.uns[f"{col}_colors"] = [palette_dict.get(c, "#d3d3d3") for c in categories]

    adata.write(adata_annotated_path)
    return adata, all_palettes

def plot_umaps(adata, adata_annotated_path, subproject_name, obs_cols, obs_titles, leiden_cols, leiden_cols_titles, figures_dir):

    umap_dir = os.path.join(figures_dir, "umap")
    leiden_umap_dir = os.path.join(umap_dir, "leiden")
    obs_umap_dir = os.path.join(umap_dir, "cat_obs")
    os.makedirs(umap_dir, exist_ok=True)
    os.makedirs(leiden_umap_dir, exist_ok=True)
    os.makedirs(obs_umap_dir, exist_ok=True)

    if obs_cols:
        sc.settings.figdir = obs_umap_dir
        for curr_col, curr_title in zip(obs_cols, obs_titles):
            if curr_col in adata.obs.columns:
                sc.pl.umap(adata, color=curr_col, title=curr_title,frameon=False,use_raw=False, save=f"_{subproject_name}_{curr_col}.png")
                if curr_col.startswith("celltype_"):
                    sc.pl.umap(adata, color=curr_col, title=curr_title,frameon=False, legend_loc = "on data",use_raw=False, save=f"_{subproject_name}_{curr_col}_no_leg.png")

    if leiden_cols:
        sc.settings.figdir = leiden_umap_dir
        if curr_col in adata.obs.columns:
            for curr_col, curr_title in zip(leiden_cols, leiden_cols_titles):
                sc.pl.umap(adata, color=curr_col, title=curr_title, frameon=False, legend_loc = "on data",use_raw=False, save=f"_{subproject_name}_{curr_col}.png")

    adata.write(adata_annotated_path)
def run_full_annotation_pipeline(subproject_name, adata_current, adata_annotated_path, leiden_res_list, json_annotations_path, palette_path, cols_to_plot, cols_to_plot_titles, tables_dir, figures_dir, palette="tab10"):
    #0.
    adata_current = add_stage_label(adata_current)
    adata_current, leiden_cols, leiden_cols_titles, celltype_cols = get_ranked_genes(adata_current, subproject_name, json_annotations_path, leiden_res_list, tables_dir, figures_dir, always_rank=True)

    # 3. Combine columns for palettes
    all_cols = cols_to_plot + leiden_cols + celltype_cols

    if os.path.exists(palette_path):
        try:
            with open(palette_path, "r") as f:
                existing_palettes = json.load(f)
        except Exception:
            existing_palettes = {}
    else:
        existing_palettes = {}
    
    print(f"Building palettes for {len(all_cols)} columns")
    adata_current, all_palettes = build_palettes(adata_current, all_cols,existing_palettes,palette=palette)

    # 5. Plot
    plot_cols = cols_to_plot + celltype_cols

    # extend titles (simple version)
    plot_titles = (cols_to_plot_titles + [f"Predicted Cell Type \n Leiden Resolution {c.replace('celltype_leiden_res_', '')}" for c in celltype_cols])
   
    plot_umaps(adata_current, subproject_name,plot_cols, plot_titles, leiden_cols,leiden_cols_titles,figures_dir)

    
    # 6. Save palettes
    print(f"Saving palettes → {palette_path}")
    with open(palette_path, "w") as f:
        json.dump(all_palettes, f, indent=2)
    # 7. Save AnnData
    adata_current.write(adata_annotated_path)

    return adata_current