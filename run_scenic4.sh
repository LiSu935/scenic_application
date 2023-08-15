#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
echo "### Starting at: $(date) ###"

source activate pyscenic37
#STEP 4:



prefix='c_pnd_combined'
working_dir='/scratch/lisu/Spencerlab/scenic_application/'
output_dir=${working_dir}'results/'${prefix}'/'
cd ${output_dir}

reg_path=${prefix}'reg.csv'
f_loom_path_scenic=${prefix}"_integrated.loom"
f_pyscenic_output=${prefix}"_pyscenic_output.loom"

pyscenic aucell \
    ${f_loom_path_scenic} \
    ${reg_path} \
    --output ${f_pyscenic_output} \
    --num_workers 20

echo "step 4 finished"



echo "### Ending at: $(date) ###"
