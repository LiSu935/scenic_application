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
path_ad_file_from_seurat_integrated="/scratch/lisu/Spencerlab/scenic_application/seurat_scenicInput/c_pnd_combined_integrated.h5ad"


#STEP 0:
python /scratch/lisu/lisu_git/scenic_application/pyscenic_pipeline2.py --prefix ${prefix} --wdir ${working_dir} --path_ad_file_from_seurat_integrated ${path_ad_file_from_seurat_integrated} 

echo "step 0 finished"
