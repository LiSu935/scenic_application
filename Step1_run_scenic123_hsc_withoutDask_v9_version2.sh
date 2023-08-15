#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --job-name=scenic_withoutDask-%j.out
#SBATCH --output=scenic_hsc-%j.out  # %j is the unique jobID
echo "### Starting at: $(date) ###"

source activate pyscenic37

prefix='c_pnd_combined'
working_dir='/scratch/lisu/Spencerlab/scenic_application/'
output_dir=${working_dir}'results/'${prefix}'/'
f_anndata_path_input="/scratch/lisu/Spencerlab/scenic_application/seurat_scenicInput/c_pnd_combined_slim.h5ad"

cd ${working_dir}
mkdir ${working_dir}'results/'
mkdir ${output_dir}

# path to unfiltered loom file (this will be created in the optional steps below)
f_loom_path_unfilt=${output_dir}${prefix}"_unfiltered.loom"

# # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
f_loom_path_scenic=${output_dir}${prefix}"_integrated.loom"

# path to anndata object, which will be updated to store Scanpy results as they are generated below
f_anndata_path=${output_dir}${prefix}"_anndata.h5ad"

# path to pyscenic output
f_pyscenic_output=${output_dir}${prefix}"_pyscenic_output.loom"

# loom output, generated from a combination of Scanpy and pySCENIC results:
f_final_loom=${output_dir}${prefix}'_scenic_integrated-output.loom'



#STEP 0:
python /scratch/lisu/lisu_git/scenic_application/pre-pyscenic_pipeline.py --prefix ${prefix} --path_ad_file_from_seurat ${f_anndata_path_input} --wdir ${working_dir}

echo "step 0 finished"


#STEP 1:

# human tfs
#f_tfs='/scratch/lisu/tools_related/scenic_AuxiliaryDatasets/allTFs_hg38.txt'

# mouse tfs
f_tfs='/scratch/lisu/tools_related/scenic_AuxiliaryDatasets/allTFs_mm.txt'

#pyscenic grn ${f_loom_path_scenic} ${f_tfs} -o ${output_dir}${prefix}adj.csv --num_workers 2

arboreto_with_multiprocessing.py ${f_loom_path_scenic} ${f_tfs} -o ${output_dir}${prefix}adj.csv --num_workers 24 --method grnboost2 --seed 777


echo "step 1 finished"

# STEP 2-3:

# ranking databases
# human
#f_db_names='/scratch/lisu/tools_related/scenic_AuxiliaryDatasets/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
#hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'

# mouse
f_db_names='/scratch/lisu/tools_related/scenic_AuxiliaryDatasets/mm9-tss-centered-5kb-10species.mc9nr.genes_vs_motifs.rankings.feather'

# motif databases
# human
#f_motif_path="/scratch/lisu/tools_related/scenic_AuxiliaryDatasets/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
#motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
#motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

# mouse
f_motif_path="/scratch/lisu/tools_related/scenic_AuxiliaryDatasets/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

pyscenic ctx ${output_dir}${prefix}adj.csv \
    ${f_db_names} \
    --annotations_fname ${f_motif_path} \
    --expression_mtx_fname ${f_loom_path_scenic} \
    --output ${output_dir}${prefix}reg.csv \
    --mask_dropouts \
    --num_workers 20 \
    --mode "custom_multiprocessing"

echo "step 2-3 finished"



echo "### Ending at: $(date) ###"
