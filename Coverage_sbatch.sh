#!/usr/bin/env bash

#SBATCH --time 1-0
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --constraint="sapphire|icelake|emerald"
#SBATCH --account=lunadiee-epi-virtualmchine
#SBATCH --output=Coverage_sbatch_%j.out
#SBATCH --error=Coverage_sbatch_%j.err

################################################################
### SBATCH header
################################################################
# Remove modules in compute node
set -e
module purge
module load bluebear # this line is required
module load bear-apps/2022b
module load SAMtools/1.17-GCC-12.2.0
module load bear-apps/2021b
module load picard/2.25.1-Java-11

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

# Sort (~ 5 minutes)
echo "Sorting .bam file for ${Individual}"
samtools sort --threads 16 ${Working_dir}/${Specondition}/bismark_alignment/${Individual}_trimmed_1_bismark_bt2_pe.bam \
-o ${Working_dir}/${Specondition}/bismark_alignment/${Individual}_sorted.bam

# indexing (~ 3 minutes)
echo "Indexing ${Individual}"
samtools index ${Working_dir}/${Specondition}/bismark_alignment/${Individual}_sorted.bam

### Estimate average reads across the genome (~ 3 minutes)
echo "Estimate average reads across genome for ${Individual}"
t=${Working_dir}/${Specondition}/bismark_alignment/${Individual}_sorted.bam
c="$(samtools depth $t | awk -v var=${Genome_size} '{sum+=$3;cnt++}END{print sum/var}')" # this is an average read per bp
echo -e "${Individual}\t$c" >> ${Working_dir}/${Specondition}/genome_coverage.txt

module purge
module load bluebear
module load bear-apps/2022b
module load BEDTools/2.30.0-GCC-12.2.0

### Create a Bedgraph summarising mapped reads at each bp (there is not longer the need of the chromosome file; ~ 3 minutes)
echo "Create the Bedgraph file for ${Individual}"
bedtools genomecov -split -bg -ibam ${Working_dir}/${Specondition}/bismark_alignment/${Individual}_sorted.bam > ${Working_dir}/${Specondition}/coverage/${Individual}.bedgraph

### Use a simple Perl command to normalise read counts by genome-wide average counts (simple coverage; ~ 3 minutes)
echo "Normalise read counts by genome-wide average counts for ${Individual}"
perl -ne 'chomp($_); @a=split(/\t/,$_);print $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]/'$c'."\t"."\n";' "${Working_dir}/${Specondition}/coverage/${Individual}"".bedgraph" > "${Working_dir}/${Specondition}/coverage/${Individual}""_norm.bedgraph"

### transform in Bigwig (~ 2 minutes)
echo "Sorting bedGraph by chromosome"
bedtools sort -i ${Working_dir}/${Specondition}/coverage/${Individual}_norm.bedgraph > ${Working_dir}/${Specondition}/coverage/${Individual}_sorted.bedGraph

### Generate Bigwig files (~ 2 minutes)
echo "Generating BigWig files"
bedGraphToBigWig ${Working_dir}/${Specondition}/coverage/${Individual}_sorted.bedGraph \
${Working_dir}/${Specondition}/annotation/bismark_index/${Genome_filename}.txt \
${Working_dir}/${Specondition}/coverage/${Individual}.bw

module purge
module load bluebear
module load bear-apps/2022b
module load R/4.3.1-foss-2022b
module load picard/2.27.5-Java-11

### calculate average insertion size (~ 3 minutes)
echo "Calculating average insertion sizes with Picard CollectInsertSizeMetrics"
java -Xmx60g -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
I=${Working_dir}/${Specondition}/bismark_alignment/${Individual}_sorted.bam \
O=${Working_dir}/${Specondition}/bismark_alignment/${Individual}_insert_size_metrics.txt \
H=${Working_dir}/${Specondition}/bismark_alignment/${Individual}_insert_size_hist.pdf


################################################################
### SBATCH footer
################################################################

# Get metrics on time taken
end=`date +%s`
echo "End time:"
date
echo "Time taken (seconds):"
expr $end - $start

