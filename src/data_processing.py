import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import pickle

def process_data():
    # Load and preprocess the dataset
    adata = sc.read_h5ad("animal_id_5.h5ad")
    genes_to_keep = [gene for gene in adata.var_names if 'Blank' not in gene]
    adata = adata[:, genes_to_keep].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=3)
    
    # Get top 10 most abundant cell types
    cell_type_counts = adata.obs['Cell_class'].value_counts()
    top_10_cell_types = cell_type_counts.head(10).index
    adata = adata[adata.obs['Cell_class'].isin(top_10_cell_types)].copy()
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Compute dimensionality reductions
    sc.pp.neighbors(adata, use_rep='X')
    sc.tl.umap(adata)
    sc.tl.tsne(adata)

    # Create clustering options and perform clustering
    n_clusters_options = [4, 5, 6, 7, 8]
    for n_clusters in n_clusters_options:
        # UMAP clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        clusters = kmeans.fit_predict(adata.obsm['X_umap'])
        adata.obs[f'umap_cluster_{n_clusters}'] = clusters.astype(str)
        
        # t-SNE clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        clusters = kmeans.fit_predict(adata.obsm['X_tsne'])
        adata.obs[f'tsne_cluster_{n_clusters}'] = clusters.astype(str)

    # Create DataFrames for visualizations
    df_umap = pd.DataFrame({
        'x': adata.obsm['spatial'][:, 0],
        'y': adata.obsm['spatial'][:, 1],
        'UMAP1': adata.obsm['X_umap'][:, 0],
        'UMAP2': adata.obsm['X_umap'][:, 1],
        'Cell Class': adata.obs['Cell_class'].values
    })

    df_tsne = pd.DataFrame({
        'x': adata.obsm['spatial'][:, 0],
        'y': adata.obsm['spatial'][:, 1],
        'tSNE1': adata.obsm['X_tsne'][:, 0],
        'tSNE2': adata.obsm['X_tsne'][:, 1],
        'Cell Class': adata.obs['Cell_class'].values
    })

    # Add cluster columns to DataFrames
    for n_clusters in n_clusters_options:
        df_umap[f'umap_cluster_{n_clusters}'] = adata.obs[f'umap_cluster_{n_clusters}'].values
        df_tsne[f'tsne_cluster_{n_clusters}'] = adata.obs[f'tsne_cluster_{n_clusters}'].values

    # Create bubble plot data
    bubble_df_dict = {}
    for n_clusters in n_clusters_options:
        # UMAP bubble plot data
        expr_df = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs.index)
        expr_df[f'umap_cluster_{n_clusters}'] = adata.obs[f'umap_cluster_{n_clusters}']
        mean_expr = expr_df.groupby(f'umap_cluster_{n_clusters}').mean().transpose()
        
        bubble_data = []
        for cluster in mean_expr.columns:
            top_genes = mean_expr[cluster].sort_values(ascending=False).head(10)
            for gene, expr in top_genes.items():
                bubble_data.append({
                    'Cluster': cluster,
                    'Gene': gene,
                    'Mean Expression': expr
                })
        
        bubble_df = pd.DataFrame(bubble_data)
        bubble_df_dict[('UMAP', n_clusters)] = bubble_df
        
        # t-SNE bubble plot data
        expr_df = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs.index)
        expr_df[f'tsne_cluster_{n_clusters}'] = adata.obs[f'tsne_cluster_{n_clusters}']
        mean_expr = expr_df.groupby(f'tsne_cluster_{n_clusters}').mean().transpose()
        
        bubble_data = []
        for cluster in mean_expr.columns:
            top_genes = mean_expr[cluster].sort_values(ascending=False).head(10)
            for gene, expr in top_genes.items():
                bubble_data.append({
                    'Cluster': cluster,
                    'Gene': gene,
                    'Mean Expression': expr
                })
        
        bubble_df = pd.DataFrame(bubble_data)
        bubble_df_dict[('TSNE', n_clusters)] = bubble_df

    # Create cell type composition data
    cell_composition_dict = {}
    for n_clusters in n_clusters_options:
        # UMAP cell composition data
        cell_counts_umap = {}
        for cluster in adata.obs[f'umap_cluster_{n_clusters}'].unique():
            cluster_mask = adata.obs[f'umap_cluster_{n_clusters}'] == cluster
            cell_type_counts = adata.obs.loc[cluster_mask, 'Cell_class'].value_counts()
            cell_counts_umap[cluster] = cell_type_counts
        
        cell_composition_dict[('UMAP', n_clusters)] = cell_counts_umap

        # t-SNE cell composition data
        cell_counts_tsne = {}
        for cluster in adata.obs[f'tsne_cluster_{n_clusters}'].unique():
            cluster_mask = adata.obs[f'tsne_cluster_{n_clusters}'] == cluster
            cell_type_counts = adata.obs.loc[cluster_mask, 'Cell_class'].value_counts()
            cell_counts_tsne[cluster] = cell_type_counts
            
        cell_composition_dict[('TSNE', n_clusters)] = cell_counts_tsne

    # Create color maps for clusters
    color_discrete_maps = {}
    for n_clusters in n_clusters_options:
        # UMAP color map
        color_discrete_maps[f'umap_cluster_{n_clusters}'] = {}
        for i in range(n_clusters):
            rgba = plt.cm.tab10(i/10)
            hex_color = '#{:02x}{:02x}{:02x}'.format(
                int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
            color_discrete_maps[f'umap_cluster_{n_clusters}'][f'UMAP_{i}'] = hex_color
        
        # t-SNE color map
        color_discrete_maps[f'tsne_cluster_{n_clusters}'] = {}
        for i in range(n_clusters):
            rgba = plt.cm.tab10(i/10)
            hex_color = '#{:02x}{:02x}{:02x}'.format(
                int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
            color_discrete_maps[f'tsne_cluster_{n_clusters}'][f'tSNE_{i}'] = hex_color

    # Save all processed data
    df_umap.to_pickle('processed_data/df_umap.pkl')
    df_tsne.to_pickle('processed_data/df_tsne.pkl')
    
    with open('processed_data/bubble_df_dict.pkl', 'wb') as f:
        pickle.dump(bubble_df_dict, f)
    
    with open('processed_data/cell_composition_dict.pkl', 'wb') as f:
        pickle.dump(cell_composition_dict, f)
    
    with open('processed_data/color_discrete_maps.pkl', 'wb') as f:
        pickle.dump(color_discrete_maps, f)
    
    with open('processed_data/n_clusters_options.pkl', 'wb') as f:
        pickle.dump(n_clusters_options, f)

if __name__ == '__main__':
    process_data() 