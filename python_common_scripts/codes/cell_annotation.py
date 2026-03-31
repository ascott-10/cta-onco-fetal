# python_common_scripts/codes/cell_annotation.py

import os
import numpy as np
import pandas as pd
import random

import scanpy as sc
import scipy.sparse as sp

import matplotlib.pyplot as plt
import decoupler as dc
import colorcet as cc

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