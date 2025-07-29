#!/bin/bash

#SBATCH --time 4-0
#SBATCH --nodes=1
#SBATCH --mem=122G
#SBATCH --account=lunadiee-epi-virtualmchine
#SBATCH --output=/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/Downstream/Methylation_around_features/Slurm_outputs/New_methylation_around_features_sbatch_%j.out
#SBATCH --error=/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/Downstream/Methylation_around_features/Slurm_outputs/New_methylation_around_features_sbatch_%j.err

# The aim of this script is to generate deeptool plots that shows methylation patterns across gene features in oak
# For publication standard plots (hopefully)

################################################################
### SBATCH header
################################################################
# Remove modules in compute node

set -e
module purge
module load bear-apps/2022a
module load deepTools/3.5.2-foss-2022a


################################################################
### Create bigwig file (signal file)
################################################################

# Filepath to bigwig files and individuals to compute
Bigwig_filepath="/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/bedgraph_files"
Bigwig_individuals=(R10 R5 R9_A R4 R3 R6_B Aa6387J Aa8673J Aa6387S Aa8673S Aa8749S P63872 P63875 P86734)

Bigwig_files=()
for individual in "${Bigwig_individuals[@]}"; do
    pattern="${Bigwig_filepath}/${individual}_${context}_min_cov_1.bw"
    # Check if the file exists and add to the list
   if [[ -f "$pattern" ]]; then
        Bigwig_files+=("$pattern")
        echo "$pattern"
    fi
done

echo ${Work_dir} Work_dir
echo ${bed_file} bed_file
echo ${context} context
echo "${Bigwig_files[@]}" Bigwig_files

# Create list of individuals to include in analysis

echo "Computing matrix for"
Out_mat_filepath=${Work_dir}/Deeptool_matrices/$(basename "${bed_file}")_${context}_matrix.tab.gz
echo $Out_mat_filepath

Option_compute_mat="true"
if [ "$Option_compute_mat" = "true" ]; then
	computeMatrix scale-regions \
		-S "${Bigwig_files[@]}" \
		-R ${bed_file} \
		--binSize 200 \
		-m 10000 \
		-p 10 \
		-b 5000 \
		-a 5000 \
		-out ${Out_mat_filepath}
fi



