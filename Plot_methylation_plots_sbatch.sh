#!/usr/bin/env bash

#SBATCH --time 1-0
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --account=lunadiee-epi-virtualmchine
#SBATCH --output=Plot_methylation_boxplots_sbatch_%j.out
#SBATCH --error=Plot_methylation_boxplots_sbatch_%j.err

################################################################
### SBATCH header
################################################################
# Remove modules in compute node
set -e
module purge
module load bluebear # this line is required
module load bear-apps/2023a
module load R-bundle-Bioconductor/3.19-foss-2023a-R-4.4.1

Rscript --vanilla ${Scripts_dir}/Methylation_calling/Plot_methylation_plots_sbatch.R ${Yaml_filepath}

