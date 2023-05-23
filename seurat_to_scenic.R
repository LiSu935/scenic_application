library(Seurat)
library(SeuratDisk)

c1 = readRDS("c1_integrated.RDS")
DefaultAssay(object = c1) <- "SCT"
c1_1 = c1
#c1_1 = DietSeurat(c1_1, assays = "SCT", dimreducs = c("umap", "pca"), graphs = TRUE)
#c1_1[['integrated']] = NULL
c1_1 = DietSeurat(c1_1, assays = "SCT",  scale.data = FALSE, dimreducs = c("umap", "pca"), graphs = TRUE)
SaveH5Seurat(c1_1, filename = "c1_slim.h5Seurat",overwrite = TRUE)
Convert("c1_slim.h5Seurat", dest = "h5ad", overwrite = TRUE)


# set variables for file paths to read from and write to:

# set a working directory
wdir = "/storage/htc/joshilab/Su_Li/Spencerlab/scenic_application/"
os.chdir( wdir )

# here to set a prefix for ease
prefix = "c1"

# =====No need if already QC by Seurat================================================= #
# ===================================================================================== #
# path to unfiltered loom file (this will be created in the optional steps below)
#f_loom_path_unfilt = "pbmc10k_unfiltered.loom" # test dataset, n=500 cells
# ===================================================================================== #

# path to loop file that has been QC and clustering by Seurat:
f_anndata_path_input = "/storage/htc/joshilab/Su_Li/Spencerlab/data/c1_slim.h5ad"

# path to unfiltered loom file (this will be created in the optional steps below)
f_loom_path_unfilt = prefix+"_unfiltered.loom"

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic = "integrated.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path = prefix+"_anndata.h5ad"

# path to pyscenic output
f_pyscenic_output = prefix+"_pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom = prefix+'_scenic_integrated-output.loom'

'/storage/htc/joshilab/Su_Li/Spencerlab/scenic_application/jobs'
