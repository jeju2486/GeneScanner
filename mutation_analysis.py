#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description='Mutation Analysis Script')
    parser.add_argument('-a', '--alignment', required=True, help='Path to the alignment file in FASTA format')
    parser.add_argument('-v', '--vcf', required=True, help='Path to the VCF file')
    parser.add_argument('-o', '--output', default='mutation_matrix.tsv', help='Output file name (default: mutation_matrix.tsv)')
    args = parser.parse_args()
    return args

def parse_alignment(alignment_file):
    sequences = {}
    for record in SeqIO.parse(alignment_file, 'fasta'):
        sequences[record.id] = str(record.seq).upper().replace('-', 'N')

    # Verify that all sequences are of the same length
    seq_lengths = set(len(seq) for seq in sequences.values())
    if len(seq_lengths) > 1:
        sys.exit("Error: Sequences are not aligned. All sequences must be of the same length.")
    
    return sequences

def translate_sequences(sequences):
    amino_acid_sequences = {}
    for sample_id, nucleotide_seq in sequences.items():
        # Ensure the sequence length is a multiple of 3 for translation
        seq_length = len(nucleotide_seq)
        trimmed_length = seq_length - (seq_length % 3)
        nucleotide_seq = nucleotide_seq[:trimmed_length]
        amino_acid_seq = str(Seq(nucleotide_seq).translate())
        amino_acid_sequences[sample_id] = amino_acid_seq
    return amino_acid_sequences

def parse_vcf(vcf_file, sample_ids):
    variant_positions = defaultdict(dict)  # {sample_id: {position: alt_base}}
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                    vcf_sample_ids = headers[9:]
                    # Map VCF sample IDs to alignment sample IDs
                    # Assuming they are the same; modify if different
                    sample_id_map = {vcf_id: aln_id for vcf_id, aln_id in zip(vcf_sample_ids, sample_ids)}
                continue
            parts = line.strip().split('\t')
            chrom, pos, _, ref, alt, _, _, _, _, *genotypes = parts
            pos = int(pos) - 1  # VCF positions are 1-based; convert to 0-based

            for vcf_id, genotype in zip(vcf_sample_ids, genotypes):
                if genotype.strip() == '1' or genotype.strip() == '1/1' or genotype.strip() == '0/1':
                    sample_id = sample_id_map.get(vcf_id)
                    if sample_id:
                        variant_positions[sample_id][pos] = alt.upper()
    return variant_positions

def categorize_samples(sample_ids):
    # Define your own logic here based on your data
    # For example, you can categorize based on sample ID patterns
    # Below is a placeholder categorization
    # Modify this function to suit your classification criteria
    resistant_samples = set()
    susceptible_samples = set()

    # Example logic:
    for sid in sample_ids:
        if 'Saitama11' in sid or 'ERR' in sid:
            resistant_samples.add(sid)
        else:
            susceptible_samples.add(sid)
    
    return resistant_samples, susceptible_samples

def analyze_mutations(sequences, amino_acid_sequences, variant_positions, resistant_samples, susceptible_samples):
    position_data = []

    sample_ids = list(sequences.keys())
    total_samples = len(sample_ids)
    amino_acid_length = len(next(iter(amino_acid_sequences.values())))
    
    # Use the first sequence as the reference
    reference_sample = next(iter(sequences))
    reference_aa_seq = amino_acid_sequences[reference_sample]
    
    for aa_pos in range(amino_acid_length):
        # Alignment position is aa_pos +1
        alignment_position = aa_pos + 1
        wildtype_aa = reference_aa_seq[aa_pos]
        
        # Initialize mutation categories
        resistant_mutations = set()
        susceptible_mutations = set()
        resistant_stop = 0
        susceptible_stop = 0
        
        for sample_id in sample_ids:
            sample_aa = amino_acid_sequences[sample_id][aa_pos]
            if sample_aa != wildtype_aa:
                if sample_id in resistant_samples:
                    resistant_mutations.add(sample_aa)
                    if sample_aa == '*':
                        resistant_stop += 1
                elif sample_id in susceptible_samples:
                    susceptible_mutations.add(sample_aa)
                    if sample_aa == '*':
                        susceptible_stop += 1
        
        # Calculate totals
        total_mutations = len(resistant_mutations) + len(susceptible_mutations)
        total_stop = resistant_stop + susceptible_stop
        
        # Calculate frequencies
        res_mut_freq = len(resistant_mutations) / len(resistant_samples) if resistant_samples else 0
        sus_mut_freq = len(susceptible_mutations) / len(susceptible_samples) if susceptible_samples else 0
        res_stop_freq = resistant_stop / len(resistant_samples) if resistant_samples else 0
        sus_stop_freq = susceptible_stop / len(susceptible_samples) if susceptible_samples else 0
        total_freq = (len(resistant_mutations) + len(susceptible_mutations)) / total_samples if total_samples else 0
        
        # Calculate percentages
        res_mut_percent = (len(resistant_mutations) / len(resistant_samples) * 100) if resistant_samples else 0
        sus_mut_percent = (len(susceptible_mutations) / len(susceptible_samples) * 100) if susceptible_samples else 0
        total_mut_percent = (total_freq * 100) if total_freq else 0
        
        # Append data
        position_data.append({
            'alignment position': alignment_position,
            'wildtype amino acid': wildtype_aa,
            'resistant/intermediate mutations': ','.join(sorted(resistant_mutations)),
            'susceptible mutations': ','.join(sorted(susceptible_mutations)),
            'resistant/intermediate stop frequency': f"{res_stop_freq:.6f}",
            'susceptible stop frequency': f"{susus_stop_freq:.6f}" if 'susus_stop_freq' in locals() else f"{sus_stop_freq:.6f}",
            'total': total_mutations,
            'resistant/intermediate': len(resistant_mutations),
            'susceptible': len(susceptible_mutations),
            'total_mutations': len(resistant_mutations) + len(susceptible_mutations),
            'resistant/intermediate_freq': f"{res_mut_freq:.6f}",
            'susceptible_freq': f"{sus_mut_freq:.6f}",
            'total_freq': f"{total_freq:.6f}",
            'resistant/intermediate_percent': f"{res_mut_percent:.2f}",
            'susceptible_percent': f"{sus_mut_percent:.2f}",
            'total_percent': f"{total_mut_percent:.2f}",
        })

    return position_data

def write_output(position_data, output_file):
    fieldnames = [
        'alignment position',
        'wildtype amino acid',
        'resistant/intermediate mutations',
        'susceptible mutations',
        'resistant/intermediate stop frequency',
        'susceptible stop frequency',
        'total',
        'resistant/intermediate',
        'susceptible',
        'total_mutations',
        'resistant/intermediate_freq',
        'susceptible_freq',
        'total_freq',
        'resistant/intermediate_percent',
        'susceptible_percent',
        'total_percent'
    ]
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for data in position_data:
            writer.writerow(data)

def main():
    args = parse_arguments()

    # Step 1: Parse the Alignment File
    sequences = parse_alignment(args.alignment)
    sample_ids = list(sequences.keys())

    # Step 2: Translate Nucleotide Sequences to Amino Acid Sequences
    amino_acid_sequences = translate_sequences(sequences)

    # Step 3: Parse the VCF File
    variant_positions = parse_vcf(args.vcf, sample_ids)

    # Step 4: Categorize Samples
    resistant_samples, susceptible_samples = categorize_samples(sample_ids)

    # Step 5: Analyze Mutations
    position_data = analyze_mutations(
        sequences, amino_acid_sequences, variant_positions,
        resistant_samples, susceptible_samples
    )

    # Step 6: Write Output
    write_output(position_data, args.output)
    print(f"Mutation matrix has been written to {args.output}")

if __name__ == '__main__':
    main()

