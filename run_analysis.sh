#!/bin/bash 

#SBATCH --job-name=mafft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=devel
#SBATCH --time=0-00:10:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/GeneScanner/mafft_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/GeneScanner/mafft_output.txt

#module load
module purge
module load Anaconda3/2024.02-1
source activate $DATA/python3_12_2
module load snp-sites/2.5.1-GCCcore-11.2.0
module load MAFFT/7.490-GCC-11.2.0-with-extensions

alignment_file="icaA_extracted_sequences_aligned_reference.fasta"

awk 'BEGIN{srand()} /^>/{ sub(/^>/,""); print $0 "," int(rand()*2) }' "$alignment_file" > group.csv

# python mutation_analysis.py \
#   -a "$alignment_file" -t both --vcf \
#   -o results/ \
#   --temp \
#   --tmp_dir ./temp_results/ \
#   --min_coverage 85 --min_identity 90 \
#   --job_id runA \

python mutation_analysis.py \
 -a "$alignment_file" -t both --vcf \
 -o grouping_results/ \
 --group group.csv \
 --min_coverage 85 --min_identity 90 \
 --temp \
 --tmp_dir ./temp_results/ \
 --job_id groupA \

#python mutation_analysis.py --help

#python mutation_analysis.py \
#  -a "$alignment_file" -t both --vcf \
#  -o results_quiet/ \
#  --quiet \
#  --min_coverage 85 --min_identity 90 \
#  --job_id runA


