# set the work dir
import os
os.chdir('/home/data/fhz/project/reference_map_PH')
print(os.getcwd())

import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="anndata")

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import omicverse as ov
ov.plot_set()

# set the display options
pd.set_option('display.max_columns', None)


adata = sc.read('data/00_merged_all.h5ad')

# add decontX contamination score to adata.obs
decontX = pd.read_csv('data/02_human_decontX_contamination.csv', index_col=None)
decontX.index = adata.obs_names.copy()
adata.obs = adata.obs.join(decontX)


ov.utils.store_layers(adata,layers='counts')
adata.layers['counts'] = adata.X

adata.var_names_make_unique()
adata.obs_names_make_unique()

# calculate some qc indicators
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]")) 

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)

# Original QC plot
n0 = adata.shape[0]
print(f'Original cell number: {n0}')

import numpy as np
tresh={
    'mito_perc': 10, 'nUMIs': 300, 'detected_genes_low': 200, 'detected_genes_high': 8000, 'decontX': 0.2
}

adata.obs['passing_mt'] = adata.obs['pct_counts_mt'] < tresh['mito_perc']
adata.obs['passing_nUMIs'] = adata.obs['total_counts'] > tresh['nUMIs']
adata.obs['passing_ngenes_high'] = adata.obs['n_genes_by_counts'] < tresh['detected_genes_high']
adata.obs['passing_ngenes_low'] = adata.obs['n_genes_by_counts'] > tresh['detected_genes_low']
adata.obs['passing_decontX'] = adata.obs['decontX_contamination'] < tresh['decontX']

print(f'Lower treshold, nUMIs: {tresh["nUMIs"]}; filtered-out-cells: {n0-np.sum(adata.obs["passing_nUMIs"])}')
print(f'Lower treshold, n genes high: {tresh["detected_genes_high"]}; filtered-out-cells: {n0-np.sum(adata.obs["passing_ngenes_high"])}')
print(f'Lower treshold, n genes low: {tresh["detected_genes_low"]}; filtered-out-cells: {n0-np.sum(adata.obs["passing_ngenes_low"])}')
print(f'Lower treshold, mito %: {tresh["mito_perc"]}; filtered-out-cells: {n0-np.sum(adata.obs["passing_mt"])}')
print(f'Lower treshold, decontX: {tresh["decontX"]}; filtered-out-cells: {n0-np.sum(adata.obs["passing_decontX"])}')

QC_test = (adata.obs['passing_mt']) & (adata.obs['passing_nUMIs']) & (adata.obs['passing_ngenes_high']) & (adata.obs['passing_ngenes_low'])
removed = QC_test.loc[lambda x : x == False]
print(f'Total cell filtered out with this last QC (and its chosen options): {n0-np.sum(QC_test)}')
adata_filtered = adata[QC_test, :].copy()
n2 = adata_filtered.shape[0]

# Store cleaned adata
print(f'Cells retained after scrublet and filtering: {n2}, {n0-n2} removed.')

# Last gene and cell filter
sc.pp.filter_cells(adata_filtered, min_genes=200)
sc.pp.filter_genes(adata_filtered, min_cells=5)

adata_filtered.write_h5ad('data/01_after_qc.h5ad', compression='gzip')
