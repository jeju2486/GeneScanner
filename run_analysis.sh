#/bin/bash 

#module load
module purge
module load Anaconda3/2024.02-1
source activate $DATA/python3_12_2

alignment_file="/data/biol-micro-genomics/kell7366/caro_project/BACT000002.aln"
vcf_file="/data/biol-micro-genomics/kell7366/caro_project/BACT000002.vcf"

python mutation_analysis.py -a "$alignment_file" -t "n" -d