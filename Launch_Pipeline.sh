#!/usr/bin/env bash

################################################################
### Workflow
################################################################

# 1. Define variables and functions
# 2. Launch Bismark pipeline
# 3. Launch DMR Calling pipeline

################################################################
### Define variables and functions
################################################################

# Variables
Yaml_reader_filepath="/rds/homes/h/hejq/bin/parse_yaml/src/parse_yaml.sh"
Yaml_filepath="/rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/Control/Test_Bruno_oak/Bruno_oak_test.yaml"

source ${Yaml_reader_filepath}
eval $(parse_yaml ${Yaml_filepath})

# Define variables
Specondition=${Species}${Condition}

Genome_filename=$(basename ${Genome_filepath})
if [[ $Genome_filename == *.gz ]]; then
    Genome_filename="${Genome_filename%.gz}"
fi

# Functions
is_job_in_queue() {
    squeue --noheader -o %j | grep -q "$1"
}

################################################################
### Create record of constants and scripts used in analysis
################################################################

echo ${Yaml_filepath}
echo ${Working_dir}

Record_dir=${Working_dir}/${Specondition}/Analysis_records

mkdir -p ${Record_dir}/Slurm_outputs
cp -rf ${Yaml_filepath} ${Record_dir}
# For completeness of records - but this is not portable for others to use the pipeine
cp -rf ${Scripts_dir}/Launch_Pipeline.sh ${Record_dir}

#cp -rf ${Scripts_dir}/Methylation_calling ${Working_dir}/${Specondition}/Analysis_records
#cp -rf ${Scripts_dir}/DMR_calling ${Working_dir}/${Specondition}/Analysis_records
#Scripts_dir=${Working_dir}/${Specondition}/Analysis_records

Yaml_filepath=${Record_dir}/$(basename ${Yaml_filepath})

################################################################
### Launch Bismark pipeline
################################################################

if [[ $Do_Bismark == TRUE ]]; then
    cd ${Record_dir}/Slurm_outputs
    mkdir -p ${Record_dir}/Methylation_calling
    cp -rf ${Scripts_dir}/Methylation_calling/Launch_Bismark_pipeline.sh ${Record_dir}/Methylation_calling
    source ${Record_dir}/Methylation_calling/Launch_Bismark_pipeline.sh
fi

sleep 10

################################################################
### Launch DMR Calling pipeline
################################################################

while is_job_in_queue "Coverage_sbatch.sh"; do
    echo "Coverage_sbatch.sh is still running - wait 1 minute"
    sleep 60
done

## this block needs fixing ...
#if find . -maxdepth 1 -type f \( -name "*.err" -o -name "*.out" -o -name "*.stats" \) | grep -q .; then
#    echo moving log files to ${Working_dir}/${Specondition}/Analysis_records/Slurm_outputs
#    for file in Coverage_sbatch_*.{err,out,stats}; do
#       mv "$file" ${Working_dir}/${Specondition}/Analysis_records/Slurm_outputs
#    done
#fi

sleep 10

if [[ $Do_DMRCalling == TRUE ]]; then
    cp -r /rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/DMR_calling/${Scripts_dir}
    source ${Scripts_dir}/DMR_calling/Launch_DMRCaller_pipeline.sh
fi

sleep 60

################################################################
### Launch Downstream Analyses
################################################################

while is_job_in_queue "Call_DMRs_sbatch.sh"; do
    echo "Call_DMRs_sbatch.sh is still running - wait 10 minute"
    sleep 600
done

sleep 60

while is_job_in_queue "Plot_methylation_boxplots_sbatch.sh"; do
    echo "Plot_methylation_boxplots_sbatch.sh is still running - wait 10 minute"
    sleep 600
done


if [[ $Do_Downstream == TRUE ]]; then
    cp -r /rds/projects/l/lunadiee-epi-virtualmchine/Slurm_pipeline/Scripts/Downstream ${Scripts_dir}
    source ${Scripts_dir}/Downstream/Launch_Downstream.sh
fi
