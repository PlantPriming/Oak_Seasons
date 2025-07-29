#!/bin/bash

#!/usr/bin/env bash

#SBATCH --time 3-0
#SBATCH --nodes=1
#SBATCH --mem=366G
#SBATCH --account=lunadiee-epi-virtualmchine
#SBATCH --output=Launch_DMRCaller_%j.out
#SBATCH --error=Launch_DMRCaller_%j.err

################################################################
### SBATCH header
################################################################
# Remove modules in compute node

set -e
module purge
module load bear-apps/2023a
module load R-bundle-Bioconductor/3.19-foss-2023a-R-4.4.1

export R_LIBS_USER=${HOME}/R/library/${EBVERSIONR}/${BB_APPS_BASE}
mkdir -p "${R_LIBS_USER}"

Rscript --vanilla /rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/DMR_calling/Do_DMRCaller_sbatch.R ${Yaml_filepath}
