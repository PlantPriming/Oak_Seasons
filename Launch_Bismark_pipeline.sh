#!/usr/bin/env bash

################################################################
### Workflow
################################################################

# Launch Bismark section of the pipeline

# 1. Copy_untar_raw_index_genome_sbatch.sh
# 2. Concat_FastQC_trim_sbatch.sh
# 3. Bismark_alignment_sbatch.sh
# 4. Coverage_sbatch.sh
# 5. Clean_intermediate_sbatch.sh (WIP)


################################################################
### Copy raw data to working dir, concatenate multi-reads
################################################################

echo "Starting Copy_untar_raw_index_genome_sbatch.sh at:"
Start_time=$(date)
echo $Start_time

mkdir -p ${Working_dir}/${Specondition}/raw
mkdir -p ${Working_dir}/${Specondition}/annotation/bismark_index

cp -rf ${Scripts_dir}/Methylation_calling/Copy_untar_raw_index_genome_sbatch.sh ${Record_dir}/Methylation_calling
sbatch --export=Yaml_filepath=${Yaml_filepath},Yaml_reader_filepath=${Yaml_reader_filepath},Specondition=${Specondition},Record_dir=${Record_dir} \
${Record_dir}/Methylation_calling/Copy_untar_raw_index_genome_sbatch.sh

sleep 10

################################################################
### Do gzip, FastQC and trim on the reads
################################################################

while is_job_in_queue "Copy_untar_raw_index_genome_sbatch.sh"; do
    echo "Copy_untar_raw_index_genome_sbatch.sh (previous step) running, wait 1 minute"
    echo "Waiting since:"
    echo $Start_time
    sleep 60
done

echo "Starting Concat_FastQC_trim_sbatch.sh at:"
Start_time=$(date)
echo $Start_time

cp -rf ${Scripts_dir}/Methylation_calling/Concat_FastQC_trim_sbatch.sh ${Record_dir}/Methylation_calling
mkdir -p ${Working_dir}/${Specondition}/trimmed

for individual in ${All_individuals_[@]}; do
    ind_name=$(eval echo \$${individual})
    sbatch --export=Yaml_filepath=${Yaml_filepath},Individual=${ind_name},Yaml_reader_filepath=${Yaml_reader_filepath} \
    ${Record_dir}/Methylation_calling/Concat_FastQC_trim_sbatch.sh
done

sleep 10

################################################################
### Bismark alignment
################################################################

echo "Starting Bismark_alignment_sbatch.sh at:"
Start_time=$(date)
echo $Start_time

while is_job_in_queue "Concat_FastQC_trim_sbatch.sh"; do
    echo "Concat_FastQC_trim_sbatch.sh (previous step) running, wait 10 minutes"
    echo "Waiting since:"
    echo $Start_time
    sleep 600
done

cp -rf ${Scripts_dir}/Methylation_calling/Bismark_alignment_sbatch.sh ${Record_dir}/Methylation_calling
mkdir -p ${Working_dir}/${Specondition}/bismark_alignment
mkdir -p ${Working_dir}/${Specondition}/meth_call/CX_reports

# SBATCH for all individuals
for individual in ${All_individuals_[@]}; do
    ind_name=$(eval echo \$${individual})
    sbatch --export=Yaml_filepath=${Yaml_filepath},Individual=${ind_name},Yaml_reader_filepath=${Yaml_reader_filepath} \
    ${Record_dir}/Methylation_calling/Bismark_alignment_sbatch.sh
done

# Remove unecessary data files to save space
# rm -rf ${Working_dir}/${Specondition}/raw
# rm -rf ${Working_dir}/${Specondition}/trimmed

################################################################
### Estimate Coverage parameters
################################################################

echo "Starting Coverage_sbatch.sh at:"
Start_time=$(date)
echo $Start_time

cp ${Scripts_dir}/Methylation_calling/Bismark_alignment_sbatch.sh ${Record_dir}/Methylation_calling/
mkdir -p ${Working_dir}/${Specondition}/coverage

while is_job_in_queue "Bismark_alignment_sbatch.sh"; do
    echo "Bismark_alignment_sbatch.sh (previous step) running, wait 10 minutes"
    echo "Waiting since:"
    echo $Start_time
    sleep 600
done

# If the genome coverage file exists, then remove it; in case of rerunning this pipeline
if [ -f "${Working_dir}/${Specondition}/genome_coverage.txt" ]; then
    rm "${Working_dir}/${Specondition}/genome_coverage.txt"
fi

Genome_size=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "${Working_dir}/${Specondition}/annotation/bismark_index/${Genome_filename}")
# SBATCH for all individuals
for individual in ${All_individuals_[@]}; do
    ind_name=$(eval echo \$${individual})
    sbatch --export=Yaml_filepath=${Yaml_filepath},Individual=${ind_name},Genome_size=${Genome_size},Genome_filename=${Genome_filename},Yaml_reader_filepath=${Yaml_reader_filepath} \
    ${Record_dir}/Methylation_calling/Coverage_sbatch.sh
done

# This can (probably) run at the same time as DMRCaller step
echo "Starting Plot_methylation_boxplots_sbatch.sh at:"
Start_time=$(date)
echo $Start_time

while is_job_in_queue "Coverage_sbatch.sh"; do
    echo "Coverage_sbatch.sh (previous step) running, wait 1 minute"
    echo "Waiting since:"
    echo $Start_time
    sleep 60
done

cp ${Scripts_dir}/Methylation_calling/Plot_methylation_plots_sbatch.sh ${Record_dir}/Methylation_calling/Plot_methylation_plots_sbatch.sh
sbatch --export=Yaml_filepath=${Yaml_filepath},Scripts_dir=${Scripts_dir} ${Record_dir}/Methylation_calling/Plot_methylation_plots_sbatch.sh

while is_job_in_queue "Plot_methylation_boxplots_sbatch.sh"; do
    echo "Plot_methylation_boxplots_sbatch.sh is still running - wait 1 minute"
    sleep 60
done

for file in Plot_methylation_boxplots_sbatch_*.{err,out,stats}; do
    mv "$file" ${Working_dir}/${Specondition}/Analysis_records/Slurm_outputs
done

# Remove unecessary intermediate data files to save space
echo removing unneeded files in meth_call
rm -rf ${Working_dir}/${Specondition}/meth_call/*.deduplicated.txt

################################################################
### Connect to DMR calling section of pipeline
################################################################

echo "Bismark pipeline scripts submitted:"
date
