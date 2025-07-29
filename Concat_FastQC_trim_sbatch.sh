#!/usr/bin/env bash

#SBATCH --time 1-0
#SBATCH --ntasks=10
#SBATCH --constraint="sapphire|icelake|emerald"
#SBATCH --account=lunadiee-epi-virtualmchine
#SBATCH --output=Concat_FastQC_trim_sbatch_%j.out
#SBATCH --error=Concat_FastQC_trim_sbatch_%j.err

################################################################
### SBATCH header
################################################################
# Remove modules in compute node
set -e
module purge
module load bluebear # this line is required
module load FastQC/0.11.9-Java-11
module load bear-apps/2021b
module load fastp/0.23.2-GCC-11.2.0 # trim illumina adapter sequences


# Get script directory and start time
echo $0
start=`date +%s`
echo "Start time:"
date

source ${Yaml_reader_filepath}
eval $(parse_yaml ${Yaml_filepath})

# Define variables
Specondition=${Species}${Condition}
Scripts_dir=${Working_dir}/${Specondition}/Analysis_records
Yaml_filepath=${Working_dir}/${Specondition}/Analysis_records/$(basename ${Yaml_filepath})


################################################################
### Main script
################################################################

echo "Working on the following indivdual:"
echo ${Individual}

echo "Concatenating multiple fastq files of the same individual and doing fastqc"
cd $(find ${Working_dir}/${Specondition}/raw -type d -name ${Individual})
gunzip *.fq.gz
for read_direction in $(seq 1 2); do
    # Concatenate reads and gzip to save space
    files=$(ls *_${read_direction}.fq | sort) # important to sort to make sure concatenated files are in matching order
    cat ${files} > ${Working_dir}/${Specondition}/raw/${Individual}_${read_direction}.fq
    gzip -c ${Working_dir}/${Specondition}/raw/${Individual}_${read_direction}.fq > ${Working_dir}/${Specondition}/raw/${Individual}_${read_direction}.fq.gz
    rm ${Working_dir}/${Specondition}/raw/${Individual}_${read_direction}.fq

    # do fastQC
    fastqc ${Working_dir}/${Specondition}/raw/${Individual}_${read_direction}.fq.gz
done

echo "Trimming and doing fastQC after"
mkdir -p ${Working_dir}/${Specondition}/trimmed/
fastp -w 8 -i ${Working_dir}/${Specondition}/raw/${Individual}_1.fq.gz -I ${Working_dir}/${Specondition}/raw/${Individual}_2.fq.gz \
-o ${Working_dir}/${Specondition}/trimmed/${Individual}_trimmed_1.fq.gz -O ${Working_dir}/${Specondition}/trimmed/${Individual}_trimmed_2.fq.gz

for read_direction in $(seq 1 2); do
    fastqc ${Working_dir}/${Specondition}/trimmed/${Individual}_trimmed_${read_direction}.fq.gz
done


################################################################
### SBATCH footer
################################################################

# Get metrics on time taken
end=`date +%s`
echo "End time:"
date
echo "Time taken (seconds):"
expr $end - $start