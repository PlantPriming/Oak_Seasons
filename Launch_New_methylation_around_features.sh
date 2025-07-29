#!/bin/bash

# The aim of this script is to generate deeptool plots that shows methylation patterns across gene features in oak
# For publication standard plots (hopefully)


################################################################
### SBATCH header
################################################################
# Remove modules in compute node

module purge
module load bluebear


Work_dir="/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/Downstream/Methylation_around_features/New_methyaltion_around_features"
Bed_file_dir=/rds/projects/l/lunadiee-rawdata-storage/MEMBRA_Datasets/Genomics/Oak/DMR_bed_files/DMR_Overlaps_CG_Seasons
Contexts=(CG CHG CHH)

# Filepath to bigwig files and individuals to compute
Bigwig_filepath="/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Data/OakAll_eCO2/analysis/objects/bedgraph_files"
Bigwig_individuals=(R10 R5 R9_A R4 R3 R6_B Aa6387J Aa8673J Aa6387S Aa8673S Aa8749S P63872 P63875 P86734)

mkdir -p ${Work_dir}
mkdir -p ${Work_dir}/Deeptool_matrices
cd ${Work_dir}

################################################################
### Launch Deeptools for each TE/Genetic element
################################################################

# Define the list of files to include for deeptools matrix creation

for bed_file in ${Bed_file_dir}/*.bed; do
    for context in "${Contexts[@]}"; do
        echo "Launching Deeptools for ${bed_file} with context ${context}"
        sbatch --export=Work_dir=${Work_dir},bed_file=${bed_file},context=${context} \
        ${Work_dir}/New_methylation_around_features_sbatch.sh
    done
done
