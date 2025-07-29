#!/usr/bin/env bash

#SBATCH --time 2-0
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --constraint="sapphire|icelake|emerald"
#SBATCH --account=lunadiee-epi-virtualmchine
#SBATCH --output=Bismark_alignment_sbatch_%j.out
#SBATCH --error=Bismark_alignment_sbatch_%j.err

################################################################
### SBATCH header
################################################################
# Remove modules in compute node
set -e
module purge
module load bluebear # this line is required
module load Python/2.7.16-GCCcore-8.3.0
module load Bismark/0.22.3-foss-2019b

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

# may need to change ulimit if there are lots of contigs in the chromosome file

# This step writes substantial temporary files - one solution to use scratch? (https://docs.bear.bham.ac.uk/bluebear/advanced_jobs/#broadwell)
# Alternatively, make sure that there is free space in the working directory, around 1TB maybe usually enough?
echo "Doing Bismark alignment with Bowtie2"
bismark --bowtie2 \
-N 1 -L 20 -p 4 -X 1000 -score_min L,0,-0.6 --parallel 4 \
-o ${Working_dir}/${Specondition}/bismark_alignment \
${Working_dir}/${Specondition}/annotation/bismark_index \
-1 ${Working_dir}/${Specondition}/trimmed/${Individual}_trimmed_1.fq.gz \
-2 ${Working_dir}/${Specondition}/trimmed/${Individual}_trimmed_2.fq.gz

sleep 10

echo "Doing Bismark deduplication"
deduplicate_bismark \
-p --output_dir ${Working_dir}/${Specondition}/bismark_alignment \
--bam ${Working_dir}/${Specondition}/bismark_alignment/${Individual}_trimmed_1_bismark_bt2_pe.bam

sleep 10

echo "Extracting cytosine methylation contexts"
## Here, consider removing the --comprehensive flag. Its in the original pipeline, but generates huge files not in the output of Marco's example oak outputs
bismark_methylation_extractor -p --parallel 6 --cytosine_report --CX --comprehensive \
--bedGraph --buffer_size 90% \
--genome_folder ${Working_dir}/${Specondition}/annotation/bismark_index \
${Working_dir}/${Specondition}/bismark_alignment/${Individual}_trimmed_1_bismark_bt2_pe.deduplicated.bam \
--output ${Working_dir}/${Specondition}/meth_call/
date

sleep 10

mv ${Working_dir}/${Specondition}/meth_call/${Individual}_trimmed_1_bismark_bt2_pe.deduplicated.CX_report.txt \
${Working_dir}/${Specondition}/meth_call/CX_reports


################################################################
### SBATCH footer
################################################################

# Get metrics on time taken
end=`date +%s`
echo "End time:"
date
echo "Time taken (seconds):"
expr $end - $start
