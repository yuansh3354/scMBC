import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import squidpy as sq
import time
import warnings
from pathlib import Path
from scipy.sparse import csr_matrix
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from sklearn.cluster import MiniBatchKMeans

warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================================
# Argument Parser
# ============================================================================
def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Universal Xenium spatial transcriptomics analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Input/Output parameters
    parser.add_argument(
        '--input_dir',
        type=str,
        required=True,
        help='Directory containing h5ad files'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        default='./results',
        help='Output directory for results [default: ./results]'
    )
    
    parser.add_argument(
        '--pattern',
        type=str,
        default='*.h5ad',
        help='File pattern for h5ad files [default: *.h5ad]'
    )
    
    # Analysis parameters
    parser.add_argument(
        '--sample_col',
        type=str,
        default='sample',
        help='Column name for sample identifier [default: sample]'
    )
    
    parser.add_argument(
        '--celltype_col',
        type=str,
        required=True,
        help='Column name for cell type annotation (required)'
    )
    
    # Neighborhood analysis parameters
    parser.add_argument(
        '--neighborhood_mode',
        type=str,
        default='knn',
        choices=['knn', 'radius'],
        help='Neighborhood definition mode [default: knn]'
    )
    
    parser.add_argument(
        '--k_neighbor',
        type=int,
        default=20,
        help='Number of nearest neighbors for knn mode [default: 20]'
    )
    
    parser.add_argument(
        '--radius',
        type=float,
        default=80.0,
        help='Radius for radius mode [default: 80.0]'
    )
    
    parser.add_argument(
        '--n_clusters',
        type=int,
        default=10,
        help='Number of clusters for neighborhood classification [default: 10]'
    )
    
    # Cell communication parameters
    parser.add_argument(
        '--comm_resource',
        type=str,
        default='cellchatdb',
        choices=['cellchatdb', 'cellphonedb', 'connectome', 'consensus'],
        help='Resource for cell communication analysis [default: cellchatdb]'
    )
    
    parser.add_argument(
        '--pval_threshold',
        type=float,
        default=0.05,
        help='P-value threshold for significant interactions [default: 0.05]'
    )
    
    parser.add_argument(
        '--expr_prop',
        type=float,
        default=0.1,
        help='Expression proportion threshold [default: 0.1]'
    )
    
    parser.add_argument(
        '--n_jobs',
        type=int,
        default=-1,
        help='Number of parallel jobs (-1 for all cores) [default: -1]'
    )
    
    # Analysis modules
    parser.add_argument(
        '--skip_neighborhood',
        action='store_true',
        help='Skip neighborhood analysis'
    )
    
    parser.add_argument(
        '--skip_communication',
        action='store_true',
        help='Skip cell communication analysis'
    )
    
    parser.add_argument(
        '--skip_visualization',
        action='store_true',
        help='Skip visualization'
    )
    
    # Visualization parameters (all optional, no defaults)
    parser.add_argument(
        '--source_celltypes',
        type=str,
        nargs='+',
        default=None,
        help='Source cell types for communication dotplot (optional, will use all if not specified)'
    )
    
    parser.add_argument(
        '--target_celltypes',
        type=str,
        nargs='+',
        default=None,
        help='Target cell types for communication dotplot (optional, will use all if not specified)'
    )
    
    parser.add_argument(
        '--ligands',
        type=str,
        nargs='+',
        default=None,
        help='Ligands of interest for filtering (optional, will show top results if not specified)'
    )
    
    parser.add_argument(
        '--receptors',
        type=str,
        nargs='+',
        default=None,
        help='Receptors of interest for filtering (optional, will show top results if not specified)'
    )
    
    parser.add_argument(
        '--top_interactions',
        type=int,
        default=50,
        help='Number of top interactions to show when no specific ligands/receptors specified [default: 50]'
    )
    
    parser.add_argument(
        '--fig_width',
        type=float,
        default=12.0,
        help='Figure width in inches [default: 12.0]'
    )
    
    parser.add_argument(
        '--fig_height',
        type=float,
        default=8.0,
        help='Figure height in inches [default: 8.0]'
    )
    
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Figure DPI [default: 300]'
    )
    
    parser.add_argument(
        '--spatial_key',
        type=str,
        default='spatial',
        help='Key for spatial coordinates in obsm [default: spatial]'
    )
    
    parser.add_argument(
        '--min_cells_per_type',
        type=int,
        default=10,
        help='Minimum number of cells per cell type to include in analysis [default: 10]'
    )
    
    return parser.parse_args()


# ============================================================================
# Data Loading Functions
# ============================================================================
def load_xenium_data(input_dir, pattern='*.h5ad'):
    """
    Load Xenium h5ad files from directory
    
    Parameters:
    -----------
    input_dir : str
        Directory containing h5ad files
    pattern : str
        File pattern to match
        
    Returns:
    --------
    adata_dict : dict
        Dictionary of AnnData objects
    """
    print("\n" + "="*80)
    print("LOADING XENIUM DATA")
    print("="*80)
    
    result_dir = Path(input_dir)
    h5ad_files = list(result_dir.glob(pattern))
    
    if len(h5ad_files) == 0:
        raise FileNotFoundError(f"No files matching pattern '{pattern}' found in {input_dir}")
    
    adata_dict = {}
    for file in h5ad_files:
        region_name = file.stem
        adata_dict[region_name] = sc.read_h5ad(file)
        print(f"Loaded: {region_name} - {adata_dict[region_name].n_obs} cells")
    
    print(f"\nTotal: {len(adata_dict)} region files loaded")
    return adata_dict


def merge_samples(adata_dict, sample_col='sample'):
    """
    Merge multiple AnnData objects and add sample identifiers
    
    Parameters:
    -----------
    adata_dict : dict
        Dictionary of AnnData objects
    sample_col : str
        Column name for sample identifier
        
    Returns:
    --------
    adata_combined : AnnData
        Merged AnnData object
    """
    print("\n" + "="*80)
    print("MERGING SAMPLES")
    print("="*80)
    
    adata_list = []
    
    for sample_key, adata in adata_dict.items():
        if sample_col not in adata.obs.columns:
            adata.obs[sample_col] = sample_key
        adata_list.append(adata)
        print(f"Processing: {sample_key}, cells: {adata.n_obs}")
    
    print("\nMerging all samples...")
    adata_combined = sc.concat(adata_list, axis=0, join='outer', index_unique='_')
    
    print(f"\nMerge completed!")
    print(f"Total cells: {adata_combined.n_obs}")
    print(f"Total samples: {len(adata_combined.obs[sample_col].unique())}")
    print(f"Sample list: {list(adata_combined.obs[sample_col].unique())}")
    
    return adata_combined


# ============================================================================
# Neighborhood Analysis Functions
# ============================================================================
def cell_neighborhood_analysis(
    adata, 
    mode='knn', 
    sample_col='sample', 
    cluster_col='celltype', 
    radius=80, 
    k_neighbor=20, 
    n_clusters=10,
    spatial_key='spatial',
    min_cells_per_type=10,
    out_dir=None
):
    """
    Perform cell neighborhood analysis
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    mode : str
        Neighborhood mode ('knn' or 'radius')
    sample_col : str
        Sample column name
    cluster_col : str
        Cell type column name
    radius : float
        Radius for radius mode
    k_neighbor : int
        Number of neighbors for knn mode
    n_clusters : int
        Number of neighborhood clusters
    spatial_key : str
        Key for spatial coordinates
    min_cells_per_type : int
        Minimum cells per cell type
    out_dir : str
        Output directory
        
    Returns:
    --------
    adata : AnnData
        Updated AnnData with neighborhood annotations
    """
    print("\n" + "="*80)
    print("CELL NEIGHBORHOOD ANALYSIS")
    print("="*80)
    print(f"Mode: {mode}")
    print(f"Cells: {adata.n_obs}, Samples: {len(adata.obs[sample_col].unique())}")
    
    # Check if spatial coordinates exist
    if spatial_key not in adata.obsm.keys():
        print(f"ERROR: Spatial coordinates '{spatial_key}' not found in adata.obsm")
        print(f"Available keys: {list(adata.obsm.keys())}")
        return adata
    
    # Filter cell types by minimum cell count
    cell_type_counts = adata.obs[cluster_col].value_counts()
    valid_cell_types = cell_type_counts[cell_type_counts >= min_cells_per_type].index
    adata_filtered = adata[adata.obs[cluster_col].isin(valid_cell_types)].copy()
    
    print(f"Filtered to {len(valid_cell_types)} cell types with ≥{min_cells_per_type} cells")
    print(f"Remaining cells: {adata_filtered.n_obs}")
    
    start_time = time.time()
    
    all_cell_types = np.sort(adata_filtered.obs[cluster_col].unique())
    cell_type_to_idx = {ct: i for i, ct in enumerate(all_cell_types)}
    n_cell_types = len(all_cell_types)
    
    sample_data = []
    samples = adata_filtered.obs[sample_col].unique()
    
    for s in samples:
        sample_start = time.time()
        sample_mask = adata_filtered.obs[sample_col] == s
        tmp = adata_filtered[sample_mask]
        x_y_coordinates = tmp.obsm[spatial_key]
        n_cells = x_y_coordinates.shape[0]
        cts = tmp.obs[cluster_col].values
        
        if mode == 'radius':
            connectivity_matrix = radius_neighbors_graph(
                x_y_coordinates, 
                radius=radius, 
                mode='connectivity'
            )
        else:
            connectivity_matrix = kneighbors_graph(
                x_y_coordinates, 
                n_neighbors=min(k_neighbor, n_cells-1), 
                mode='connectivity'
            )
        
        neighborhood_matrix = np.zeros((n_cells, n_cell_types), dtype=np.float32)
        
        for i in range(n_cells):
            neighbors = connectivity_matrix[i].indices
            if len(neighbors) > 0:
                neighbor_types = cts[neighbors]
                for ct in neighbor_types:
                    neighborhood_matrix[i, cell_type_to_idx[ct]] += 1
        
        row_sums = neighborhood_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        neighborhood_matrix = neighborhood_matrix / row_sums
        
        sample_data.append({
            'sample': s,
            'matrix': neighborhood_matrix,
            'cells': tmp.obs.index.tolist()
        })
        
        print(f"Sample {s} processed in {time.time() - sample_start:.2f}s")
    
    all_matrices = np.vstack([d['matrix'] for d in sample_data])
    
    print("\nPerforming clustering...")
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters, 
        random_state=42, 
        batch_size=min(10000, len(all_matrices)),
        n_init=10
    )
    cluster_labels = kmeans.fit_predict(all_matrices)
    
    neighborhood_dict = {}
    start_idx = 0
    for d in sample_data:
        n = len(d['cells'])
        for cell, label in zip(d['cells'], cluster_labels[start_idx:start_idx + n]):
            neighborhood_dict[cell] = f'Niche_{label}'
        start_idx += n
    
    adata.obs['neighborhood'] = adata.obs.index.map(neighborhood_dict)
    
    print(f"\nNeighborhood analysis completed in {time.time() - start_time:.2f}s")
    print(f"Identified {n_clusters} neighborhood types")
    
    if out_dir:
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        
        neighborhood_df = pd.DataFrame({
            'cell_id': list(neighborhood_dict.keys()),
            'neighborhood': list(neighborhood_dict.values())
        })
        neighborhood_df.to_csv(out_path / 'neighborhood_assignments.csv', index=False)
        print(f"Results saved to {out_path}")
    
    return adata


# ============================================================================
# Cell Communication Analysis Functions
# ============================================================================
def cell_communication_analysis(
    adata,
    sample_col='sample',
    celltype_col='celltype',
    resource_name='cellchatdb',
    expr_prop=0.1,
    n_jobs=-1,
    min_cells_per_type=10,
    out_dir=None
):
    """
    Perform cell-cell communication analysis using LIANA/CellPhoneDB
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    sample_col : str
        Sample column name
    celltype_col : str
        Cell type column name
    resource_name : str
        Communication database resource
    expr_prop : float
        Expression proportion threshold
    n_jobs : int
        Number of parallel jobs
    min_cells_per_type : int
        Minimum cells per cell type
    out_dir : str
        Output directory
        
    Returns:
    --------
    adata : AnnData
        Updated AnnData with communication results
    """
    print("\n" + "="*80)
    print("CELL COMMUNICATION ANALYSIS")
    print("="*80)
    
    try:
        import liana as li
        from liana.method import cellphonedb
    except ImportError:
        print("ERROR: liana package not installed. Install with: pip install liana")
        return adata
    
    # Filter cell types by minimum cell count
    cell_type_counts = adata.obs[celltype_col].value_counts()
    valid_cell_types = cell_type_counts[cell_type_counts >= min_cells_per_type].index
    adata_filtered = adata[adata.obs[celltype_col].isin(valid_cell_types)].copy()
    
    print(f"Filtered to {len(valid_cell_types)} cell types with ≥{min_cells_per_type} cells")
    
    cpdb_res_sample = {}
    samples = adata_filtered.obs[sample_col].drop_duplicates()
    
    for sample in samples:
        print(f"\nProcessing sample: {sample}")
        tmp = adata_filtered[adata_filtered.obs[sample_col] == sample, :]
        
        # Check if sample has enough cells
        if tmp.n_obs < 10:
            print(f"Skipping sample {sample}: too few cells ({tmp.n_obs})")
            continue
            
        cellphonedb(
            tmp,
            n_jobs=n_jobs,
            groupby=celltype_col,
            resource_name=resource_name,
            expr_prop=expr_prop,
            use_raw=False,
            key_added='cpdb_res',
            verbose=True
        )
        
        cpdb_res_sample[sample] = tmp.uns['cpdb_res']
        print(f"Sample {sample}: {len(tmp.uns['cpdb_res'])} interactions detected")
    
    if len(cpdb_res_sample) == 0:
        print("WARNING: No samples had sufficient data for communication analysis")
        return adata
    
    result_list = []
    for k, v in cpdb_res_sample.items():
        v['sample'] = k
        v['source_target_celltype'] = v['source'] + "|" + v['target']
        v['L_R'] = v['ligand'] + "->" + v['receptor']
        result_list.append(v)
    
    combined_results = pd.concat(result_list)
    adata.uns['cpdb_res'] = combined_results
    
    if out_dir:
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        combined_results.to_csv(out_path / 'cpdb_results.csv', index=False)
        print(f"\nResults saved to {out_path / 'cpdb_results.csv'}")
    
    print(f"\nTotal interactions: {len(combined_results)}")
    print(f"Unique L-R pairs: {len(combined_results['L_R'].unique())}")
    print(f"Unique cell type pairs: {len(combined_results['source_target_celltype'].unique())}")
    
    return adata


# ============================================================================
# Visualization Functions
# ============================================================================
def visualize_communication(
    adata,
    source_labels=None,
    target_labels=None,
    ligands=None,
    receptors=None,
    pval_threshold=0.05,
    top_interactions=50,
    figure_size=(12, 8),
    dpi=300,
    out_dir=None
):
    """
    Create communication dotplot visualization
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with communication results
    source_labels : list, optional
        Source cell types to display
    target_labels : list, optional
        Target cell types to display
    ligands : list, optional
        Ligands of interest
    receptors : list, optional
        Receptors of interest
    pval_threshold : float
        P-value threshold
    top_interactions : int
        Number of top interactions to show
    figure_size : tuple
        Figure size (width, height)
    dpi : int
        Figure DPI
    out_dir : str
        Output directory
    """
    print("\n" + "="*80)
    print("GENERATING COMMUNICATION VISUALIZATION")
    print("="*80)
    
    try:
        import liana as li
    except ImportError:
        print("ERROR: liana package not installed")
        return
    
    if 'cpdb_res' not in adata.uns:
        print("ERROR: No communication results found. Run communication analysis first.")
        return
    
    cpdb_res = adata.uns['cpdb_res']
    
    # If no specific filters provided, use top interactions
    if ligands is None and receptors is None:
        print(f"No specific ligands/receptors specified, showing top {top_interactions} interactions")
        
        # Sort by p-value and get top interactions
        cpdb_filtered = cpdb_res[cpdb_res['cellphone_pvals'] <= pval_threshold].copy()
        cpdb_filtered = cpdb_filtered.sort_values('cellphone_pvals').head(top_interactions)
        
        ligands = cpdb_filtered['ligand'].unique().tolist()
        receptors = cpdb_filtered['receptor'].unique().tolist()
    
    # If no cell type labels specified, use all
    if source_labels is None:
        source_labels = sorted(cpdb_res['source'].unique().tolist())
        print(f"Using all {len(source_labels)} source cell types")
    
    if target_labels is None:
        target_labels = sorted(cpdb_res['target'].unique().tolist())
        print(f"Using all {len(target_labels)} target cell types")
    
    print(f"Creating dotplot for filtered interactions...")
    print(f"Source cell types: {len(source_labels)}")
    print(f"Target cell types: {len(target_labels)}")
    print(f"Ligands: {len(ligands) if ligands else 0}")
    print(f"Receptors: {len(receptors) if receptors else 0}")
    
    # Create filter function
    if ligands is not None and receptors is not None:
        filter_fun = lambda x: (x['cellphone_pvals'] <= pval_threshold) & \
                              (x['ligand'].isin(ligands)) & \
                              (x['receptor'].isin(receptors))
    elif ligands is not None:
        filter_fun = lambda x: (x['cellphone_pvals'] <= pval_threshold) & \
                              (x['ligand'].isin(ligands))
    elif receptors is not None:
        filter_fun = lambda x: (x['cellphone_pvals'] <= pval_threshold) & \
                              (x['receptor'].isin(receptors))
    else:
        filter_fun = lambda x: x['cellphone_pvals'] <= pval_threshold
    
    try:
        fig = li.pl.dotplot(
            adata=adata,
            colour='lr_means',
            size='cellphone_pvals',
            inverse_size=True,
            source_labels=source_labels,
            target_labels=target_labels,
            figure_size=figure_size,
            filter_fun=filter_fun,
            uns_key='cpdb_res',
            return_fig=True
        )
        
        if out_dir:
            out_path = Path(out_dir)
            out_path.mkdir(parents=True, exist_ok=True)
            
            fig.save(out_path / 'ligand_receptor_dotplot.png', dpi=dpi)
            fig.save(out_path / 'ligand_receptor_dotplot.pdf', dpi=dpi)
            print(f"\nFigures saved to {out_path}")
        
        print("Visualization completed")
        
    except Exception as e:
        print(f"WARNING: Could not create dotplot: {e}")
        print("Saving filtered results to CSV instead")
        
        if out_dir:
            out_path = Path(out_dir)
            out_path.mkdir(parents=True, exist_ok=True)
            
            # Filter and save results
            filtered_res = cpdb_res[cpdb_res['cellphone_pvals'] <= pval_threshold].copy()
            if ligands is not None:
                filtered_res = filtered_res[filtered_res['ligand'].isin(ligands)]
            if receptors is not None:
                filtered_res = filtered_res[filtered_res['receptor'].isin(receptors)]
            
            filtered_res.to_csv(out_path / 'filtered_interactions.csv', index=False)
            print(f"Filtered interactions saved to {out_path / 'filtered_interactions.csv'}")


def plot_spatial_distribution(adata, celltype_col='celltype', spatial_key='spatial', out_dir=None):
    """
    Plot spatial distribution of cell types
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    celltype_col : str
        Cell type column name
    spatial_key : str
        Key for spatial coordinates
    out_dir : str
        Output directory
    """
    print("\n" + "="*80)
    print("PLOTTING SPATIAL DISTRIBUTION")
    print("="*80)
    
    if spatial_key not in adata.obsm.keys():
        print(f"WARNING: Spatial coordinates '{spatial_key}' not found, skipping spatial plot")
        return
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    try:
        sc.pl.embedding(
            adata,
            basis=spatial_key,
            color=celltype_col,
            ax=ax,
            show=False
        )
        
        if out_dir:
            out_path = Path(out_dir)
            out_path.mkdir(parents=True, exist_ok=True)
            plt.savefig(out_path / 'spatial_distribution.png', dpi=300, bbox_inches='tight')
            plt.savefig(out_path / 'spatial_distribution.pdf', bbox_inches='tight')
            print(f"Spatial plot saved to {out_path}")
        
        plt.close()
    except Exception as e:
        print(f"WARNING: Could not create spatial plot: {e}")
        plt.close()


# ============================================================================
# Main Pipeline
# ============================================================================
def main():
    """Main analysis pipeline"""
    args = parse_arguments()
    
    print("\n" + "="*80)
    print("XENIUM SPATIAL TRANSCRIPTOMICS ANALYSIS PIPELINE")
    print("="*80)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"File pattern: {args.pattern}")
    print(f"Cell type column: {args.celltype_col}")
    
    start_time = time.time()
    
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load data
    adata_dict = load_xenium_data(args.input_dir, args.pattern)
    
    # Merge samples
    adata = merge_samples(adata_dict, sample_col=args.sample_col)
    
    # Check if cell type column exists
    if args.celltype_col not in adata.obs.columns:
        print(f"\nERROR: Cell type column '{args.celltype_col}' not found")
        print(f"Available columns: {list(adata.obs.columns)}")
        return
    
    print(f"\nCell type distribution:")
    print(adata.obs[args.celltype_col].value_counts().head(20))
    
    # Neighborhood analysis
    if not args.skip_neighborhood:
        adata = cell_neighborhood_analysis(
            adata,
            mode=args.neighborhood_mode,
            sample_col=args.sample_col,
            cluster_col=args.celltype_col,
            radius=args.radius,
            k_neighbor=args.k_neighbor,
            n_clusters=args.n_clusters,
            spatial_key=args.spatial_key,
            min_cells_per_type=args.min_cells_per_type,
            out_dir=args.output_dir
        )
    
    # Communication analysis
    if not args.skip_communication:
        adata = cell_communication_analysis(
            adata,
            sample_col=args.sample_col,
            celltype_col=args.celltype_col,
            resource_name=args.comm_resource,
            expr_prop=args.expr_prop,
            n_jobs=args.n_jobs,
            min_cells_per_type=args.min_cells_per_type,
            out_dir=args.output_dir
        )
    
    # Visualization
    if not args.skip_visualization:
        plot_spatial_distribution(
            adata,
            celltype_col=args.celltype_col,
            spatial_key=args.spatial_key,
            out_dir=args.output_dir
        )
        
        if 'cpdb_res' in adata.uns:
            visualize_communication(
                adata,
                source_labels=args.source_celltypes,
                target_labels=args.target_celltypes,
                ligands=args.ligands,
                receptors=args.receptors,
                pval_threshold=args.pval_threshold,
                top_interactions=args.top_interactions,
                figure_size=(args.fig_width, args.fig_height),
                dpi=args.dpi,
                out_dir=args.output_dir
            )
    
    # Save final object
    adata.write_h5ad(output_path / 'xenium_analyzed.h5ad')
    print(f"\nFinal AnnData saved to {output_path / 'xenium_analyzed.h5ad'}")
    
    total_time = time.time() - start_time
    print("\n" + "="*80)
    print(f"PIPELINE COMPLETED IN {total_time/60:.2f} MINUTES")
    print("="*80)


if __name__ == '__main__':
    main()
