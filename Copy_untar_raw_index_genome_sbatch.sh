#!/usr/bin/env bash

#SBATCH --time 1-0
#SBATCH --ntasks=10
#SBATCH --constraint="sapphire|icelake|emerald"
#SBATCH --account=lunadiee-epi-virtualmchine
#SBATCH --output=Copy_untar_raw_index_genome_sbatch_%j.out
#SBATCH --error=Copy_untar_raw_index_genome_sbatch_%j.err

################################################################
### SBATCH header
################################################################
# Remove modules in compute node
set -e
module purge
module load bluebear # this line is required
module load bear-apps/2022b
module load SAMtools/1.17-GCC-12.2.0


# Get script directory and start time
echo $0
start=`date +%s`
echo "Start time:"
date

source ${Yaml_reader_filepath}
eval $(parse_yaml ${Yaml_filepath})

Genome_filename=$(basename ${Genome_filepath})
if [[ $Genome_filename == *.gz ]]; then
    Genome_filename="${Genome_filename%.gz}"
fi


################################################################
### Main script
################################################################

echo "Copying raw data to working directory"
for individual in ${All_individuals_[@]}; do
    eval cp -rf ${Raw_data_dir}/\$${individual} ${Working_dir}/${Specondition}/raw
    eval echo \$${individual} copied to ${Working_dir}/${Specondition}/raw
done

echo "Getting the genome and gunzipping if needed"
cp -f ${Genome_filepath} ${Working_dir}/${Specondition}/annotation/bismark_index/
if [[ ${Genome_filepath} == *.gz ]]; then
    gunzip -f ${Working_dir}/${Specondition}/annotation/bismark_index/$(basename ${Genome_filepath})
fi

echo "Indexing genome and making txt file of genome size"
mkdir -p ${Working_dir}/${Specondition}/annotation/bismark_index/
samtools faidx ${Working_dir}/${Specondition}/annotation/bismark_index/${Genome_filename}
cut -f1,2 ${Working_dir}/${Specondition}/annotation/bismark_index/${Genome_filename}.fai > ${Working_dir}/${Specondition}/annotation/bismark_index/${Genome_filename}.txt

module purge
module load bear-apps/2021b
module load Bismark/0.24.2-foss-2021b

echo "Starting Bismark genome preparation"
bismark_genome_preparation --verbose --bowtie2 --parallel 4 ${Working_dir}/${Specondition}/annotation/bismark_index


################################################################
### SBATCH footer
################################################################

# Get metrics on time taken
end=`date +%s`
echo "End time:"
date
echo "Time taken (seconds):"
expr $end - $start