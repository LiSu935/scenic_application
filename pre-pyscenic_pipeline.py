# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='pre-pyscenic pipeline, generating loom file for scenic step1')
parser.add_argument('--prefix', type=str, default='SULI', help='prefix: prefix for the output files')
parser.add_argument('--path_ad_file_from_seurat', type=str, default="/storage/htc/joshilab/Su_Li/Spencerlab/data/c1_slim.h5ad", help='SCT mtx from seurat containing all features')
parser.add_argument('--wdir', type=str, default="/storage/htc/joshilab/Su_Li/Spencerlab/scenic_application/", help='working dir of the scenic, major dir.')

args = parser.parse_args()
print(args)

# set variables for file paths to read from and write to:

# set a working directory
wdir = args.wdir
os.chdir( wdir )

# here to set a prefix for ease
prefix = args.prefix

output_dir = "results/"+prefix+"/"

if not os.path.exists('output_dir'):
   os.makedirs('output_dir')


# =====No need if already QC by Seurat================================================= #
# ===================================================================================== #
# path to unfiltered loom file (this will be created in the optional steps below)
#f_loom_path_unfilt = "pbmc10k_unfiltered.loom" # test dataset, n=500 cells
# ===================================================================================== #

# path to anndata file that has been QC and clustering by Seurat:
f_anndata_path_input = args.path_ad_file_from_seurat

# path to unfiltered loom file (this will be created in the optional steps below)
f_loom_path_unfilt = output_dir+prefix+"_unfiltered.loom"

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = output_dir+prefix+"_integrated.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = output_dir+prefix+"_anndata.h5ad"

# path to pyscenic output
f_pyscenic_output = output_dir+prefix+"_pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = output_dir+prefix+'_scenic_integrated-output.loom'

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 20


# Expression data import
adata = sc.read_h5ad( f_anndata_path_input )

row_attrs = { 
    "Gene": np.array(adata.var.index) ,
}
col_attrs = { 
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

lp.create( f_loom_path_unfilt, adata.X.transpose(), row_attrs, col_attrs )

# read unfiltered data from a loom file
adata = sc.read_loom( f_loom_path_unfilt )

nCountsPerGene = np.sum(adata.X, axis=0)
nCellsPerGene = np.sum(adata.X>0, axis=0)

# Show info
print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())

nCells=adata.X.shape[0]

# pySCENIC thresholds
minCountsPerGene=3*.01*nCells # 3 counts in 1% of cells
print("minCountsPerGene: ", minCountsPerGene)

minSamples=.01*nCells # 1% of cells
print("minSamples: ", minSamples)

# simply compute the number of genes per cell (computers 'n_genes' column)
#sc.pp.filter_cells(adata, min_genes=0)
# mito and genes/counts cuts
mito_genes = adata.var_names.str.startswith('Mt-')
# for each cell compute fraction of counts in mito genes vs. all genes
adata.obs['percent_mito'] = np.sum(
    np.matrix(adata[:, mito_genes].X), axis=1).A1 / np.sum(np.matrix(adata.X), axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = np.matrix(adata.X).sum(axis=1).A1

adata.obs["n_genes"]= np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten()

# update the anndata file:
adata.write( f_anndata_path )

# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
# here the loom file for scenic, it requires all the features not just variable features.
# Additionally, here I used the datasets from SCT assay, before integration steps in Seurat. 
lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)


# ===================================================================================== # 
# then go to SCENIC steps
# ===================================================================================== # 





