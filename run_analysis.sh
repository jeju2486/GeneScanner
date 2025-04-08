#!/bin/bash 

#SBATCH --job-name=mafft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=devel
#SBATCH --time=0-00:10:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/caro_project/mafft_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/caro_project/mafft_output.txt

#module load
module purge
module load Anaconda3/2024.02-1
source activate $DATA/python3_12_2
module load snp-sites/2.5.1-GCCcore-11.2.0
module load MAFFT/7.490-GCC-11.2.0-with-extensions

#alignment_file="/data/biol-micro-genomics/kell7366/caro_project/data/BIGSdb_1462203_9129229247_35550_aligned.fas"
alignment_file="/data/biol-micro-genomics/kell7366/caro_project/aligned_lessgap.fasta"
#alignment_file="/data/biol-micro-genomics/kell7366/caro_project/caro_test_proteins_aligned.fas"

#awk '/^>/{if(seq && seq != "-") {print header; print seq}; header=$0; seq=""} !/^>/{seq=seq $0} END {if(seq && seq != "-") {print header; print seq}}' BIGSdb_1062038_2648453179_69396.fas > filtered.fasta

#awk '/^>/{header=$0; getline seq; if(substr(seq,1,3)=="ATG") {print header; print seq}}' filtered.fasta > filtered_trimmed.fasta

#filter_blank: 7238 filter_incompete: 7221 original: 7293


#python mutation_analysis.py --help

#mafft --auto --leavegappyregion BIGSdb_1062038_2648453179_69396.fas > aligned.fasta

#mafft --thread 8 filtered.fasta > aligned.fasta

#mafft --auto --leavegappyregion --thread 8 filtered_trimmed.fasta > aligned_lessgap.fasta

python mutation_analysis.py -a "$alignment_file" -t "both" --vcf --temp --job_id "caro_test"

#module load 
module load HMMER/3.3.2-gompic-2020b