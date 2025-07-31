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

# awk 'BEGIN{srand()} /^>/{ sub(/^>/,""); print $0 "," int(rand()*2) }' "$alignment_file" > group.csv

# alignment_file="/data/biol-micro-genomics/kell7366/GeneScanner/simulation_data/syn_nonsyn_stop_isolates.fasta"

# python mutation_analysis.py \
#   -a "$alignment_file" -t both   \
#   -o results/ \
#   --temp \
#   --tmp_dir ./temp_results/ \
#   --min_coverage 85 --min_identity 90 \
#   --job_id simulation_syn_nonsyn \

# alignment_file="/data/biol-micro-genomics/kell7366/GeneScanner/simulation_data/indel_isolates.fasta"

#   python mutation_analysis.py \
#   -a "$alignment_file" -t both   \
#   -o results/ \
#   --temp \
#   --tmp_dir ./temp_results/ \
#   --min_coverage 85 --min_identity 90 \
#   --job_id simulation_indel \

# All mutations example

alignment_file="/data/biol-micro-genomics/kell7366/GeneScanner/simulation_data/all_mut_isolates.fasta"

  python mutation_analysis.py \
  -a "$alignment_file" -t both   \
  -o results/ \
  --temp \
  --tmp_dir ./temp_results/ \
  --min_coverage 85 --min_identity 90 \
  --job_id simulation_all \

python genescanner_viewer.py -i ./results/simulation_all_both_mutation_analysis.xlsx -o ./results/all_mut_figure.png -t 1

# groupming exmaple

# alignment_file="/data/biol-micro-genomics/kell7366/GeneScanner/simulation_data/all_isolates.mut.group.fa"

# awk '
#   BEGIN { count = 0 }
#   /^>/ {
#     sub(/^>/, "")
#     count++
#     grp = (count <= 100 ? 1 : 2)
#     print $0 "," grp
#   }
# ' "$alignment_file" > group.csv

# python mutation_analysis.py \
#  -a "$alignment_file" -t both --vcf \
#  -o grouping_results/ \
#  --group group.csv \
#  --min_coverage 0 --min_identity 0 \
#  --temp \
#  --tmp_dir ./temp_results/ \
#  --job_id group \

#python mutation_analysis.py --help

#python mutation_analysis.py \
#  -a "$alignment_file" -t both --vcf \
#  -o results_quiet/ \
#  --quiet \
#  --min_coverage 85 --min_identity 90 \
#  --job_id runA


