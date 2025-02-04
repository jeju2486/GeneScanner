#/bin/bash 

#module load
module purge
module load Anaconda3/2024.02-1
source activate $DATA/python3_12_2
module load snp-sites/2.5.1-GCCcore-11.2.0
module load MAFFT/7.490-GCC-11.2.0-with-extensions

alignment_file="/data/biol-micro-genomics/kell7366/caro_project/data/BIGSdb_1462935_6708830946_63728_aligned.fas"
vcf_file="/data/biol-micro-genomics/kell7366/caro_project/data/BACT000002.vcv"


#python mutation_analysis.py --help
python mutation_analysis.py -a "$alignment_file" -t "both" -d --vcf