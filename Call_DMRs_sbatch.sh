#!/usr/bin/env bash

#SBATCH --time 3-0
#SBATCH --ntasks=1
#SBATCH --mem=366G
#SBATCH --constraint="sapphire|emerald"
#SBATCH --account=lunadiee-epi-virtualmchine
#SBATCH --output=Call_DMRs_sbatch_%j.out
#SBATCH --error=Call_DMRs_sbatch_%j.err

################################################################
### SBATCH header
################################################################
# Remove modules in compute node
set -e
module purge
module load bluebear # this line is required
module load bear-apps/2023a
module load R-bundle-Bioconductor/3.19-foss-2023a-R-4.4.1

echo $Control_group
echo $Test_group
echo $Yaml_filepath
echo ${Scripts_dir}

Rscript --vanilla ${Scripts_dir}/DMR_calling/Call_DMRs_sbatch.R ${Yaml_filepath} ${Control_group} ${Test_group}
