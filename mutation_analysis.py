#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
import sys
import logging

def parse_arguments():
    parser = argparse.ArgumentParser(description='Mutation Analysis Script')
    parser.add_argument('-a', '--alignment', required=True, help='Path to the alignment file in FASTA format')
    parser.add_argument('-v', '--vcf', required=True, help='Path to the VCF file')
    parser.add_argument('-o', '--output', default='mutation_matrix.tsv', help='Output file name (default: mutation_matrix.tsv)')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug mode for detailed logging')
    args = parser.parse_args()
    return args

def setup_logging(debug=False):
    """Set up the logging configuration."""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

def parse_alignment(alignment_file):
    logging.info(f"Parsing alignment file: {alignment_file}")
    sequences = {}
    try:
        for record in SeqIO.parse(alignment_file, 'fasta'):
            sequences[record.id] = str(record.seq).upper().replace('-', 'N')
            logging.debug(f"Loaded sequence for sample '{record.id}' with length {len(record.seq)}")
    except Exception as e:
        logging.error(f"Error parsing alignment file: {e}")
        sys.exit(1)

    # Verify that all sequences are of the same length
    seq_lengths = set(len(seq) for seq in sequences.values())
    if len(seq_lengths) > 1:
        logging.error("Sequences are not aligned. All sequences must be of the same length.")
        sys.exit("Error: Sequences are not aligned. All sequences must be of the same length.")
    logging.info(f"All sequences are aligned with length {seq_lengths.pop()} nucleotides.")
    return sequences

def translate_sequences(sequences):
    logging.info("Translating nucleotide sequences to amino acid sequences.")
    amino_acid_sequences = {}
    for sample_id, nucleotide_seq in sequences.items():
        # Ensure the sequence length is a multiple of 3 for translation
        seq_length = len(nucleotide_seq)
        trimmed_length = seq_length - (seq_length % 3)
        if trimmed_length != seq_length:
            logging.warning(f"Sequence '{sample_id}' length {seq_length} is not a multiple of 3. Trimming to {trimmed_length}.")
            nucleotide_seq = nucleotide_seq[:trimmed_length]
        try:
            amino_acid_seq = str(Seq(nucleotide_seq).translate())
            amino_acid_sequences[sample_id] = amino_acid_seq
            logging.debug(f"Translated sequence for sample '{sample_id}' to amino acids with length {len(amino_acid_seq)}")
        except Exception as e:
            logging.error(f"Error translating sequence for sample '{sample_id}': {e}")
            sys.exit(1)
    logging.info("Translation of all sequences completed.")
    return amino_acid_sequences

def parse_vcf(vcf_file, sample_ids):
    logging.info(f"Parsing VCF file: {vcf_file}")
    variant_positions = defaultdict(dict)  # {sample_id: {position: alt_base}}
    try:
        with open(vcf_file, 'r') as vcf:
            for line in vcf:
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        headers = line.strip().split('\t')
                        vcf_sample_ids = headers[9:]
                        # Map VCF sample IDs to alignment sample IDs
                        sample_id_map = {vcf_id: aln_id for vcf_id, aln_id in zip(vcf_sample_ids, sample_ids)}
                        logging.debug("Sample ID mapping between VCF and alignment established.")
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    logging.warning(f"Skipping malformed VCF line: {line.strip()}")
                    continue
                chrom, pos, _, ref, alt, _, _, _, _, *genotypes = parts
                pos = int(pos) - 1  # VCF positions are 1-based; convert to 0-based

                for vcf_id, genotype in zip(vcf_sample_ids, genotypes):
                    genotype = genotype.strip()
                    if genotype in {'1', '1/1', '0/1', '1|0', '0|1', '1|1'}:
                        sample_id = sample_id_map.get(vcf_id)
                        if sample_id:
                            variant_positions[sample_id][pos] = alt.upper()
                            logging.debug(f"Recorded variant for sample '{sample_id}' at position {pos + 1}: {ref} -> {alt}")
    except FileNotFoundError:
        logging.error(f"VCF file not found: {vcf_file}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error parsing VCF file: {e}")
        sys.exit(1)

    logging.info(f"Completed parsing VCF file. Total samples with variants: {len(variant_positions)}")
    return variant_positions

def categorize_samples(sample_ids):
    logging.info("Categorizing samples into resistant/intermediate and susceptible groups.")
    # Define your own logic here based on your data
    # For this example, we'll assume sample IDs containing 'Saitama11' are resistant
    resistant_samples = set()
    susceptible_samples = set()

    for sid in sample_ids:
        if 'Saitama11' in sid or 'ERR' in sid:
            resistant_samples.add(sid)
            logging.debug(f"Sample '{sid}' categorized as Resistant/Intermediate.")
        else:
            susceptible_samples.add(sid)
            logging.debug(f"Sample '{sid}' categorized as Susceptible.")

    logging.info(f"Total Resistant/Intermediate samples: {len(resistant_samples)}")
    logging.info(f"Total Susceptible samples: {len(susceptible_samples)}")
    return resistant_samples, susceptible_samples

def analyze_mutations(sequences, amino_acid_sequences, variant_positions, resistant_samples, susceptible_samples):
    logging.info("Starting mutation analysis.")
    position_data = []

    sample_ids = list(sequences.keys())
    total_samples = len(sample_ids)
    amino_acid_length = len(next(iter(amino_acid_sequences.values())))
    logging.debug(f"Amino acid sequence length: {amino_acid_length}")

    # Use the first sequence as the reference
    reference_sample = sample_ids[0]
    reference_aa_seq = amino_acid_sequences[reference_sample]
    logging.debug(f"Reference sample: '{reference_sample}'")

    for aa_pos in range(amino_acid_length):
        # Alignment position is aa_pos +1
        alignment_position = aa_pos + 1
        wildtype_aa = reference_aa_seq[aa_pos]
        logging.debug(f"Analyzing position {alignment_position}: Wildtype AA = '{wildtype_aa}'")

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
                        logging.debug(f"Stop codon detected in Resistant sample '{sample_id}' at position {alignment_position}")
                elif sample_id in susceptible_samples:
                    susceptible_mutations.add(sample_aa)
                    if sample_aa == '*':
                        susceptible_stop += 1
                        logging.debug(f"Stop codon detected in Susceptible sample '{sample_id}' at position {alignment_position}")

        # Calculate totals
        total_mutations = len(resistant_mutations) + len(susceptible_mutations)
        logging.debug(f"Position {alignment_position}: Total mutations = {total_mutations}")

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
            'susceptible stop frequency': f"{sus_stop_freq:.6f}",
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

        # Optionally, log progress every 100 positions
        if (aa_pos + 1) % 100 == 0:
            logging.info(f"Processed {aa_pos + 1}/{amino_acid_length} amino acid positions.")

    logging.info("Mutation analysis completed.")
    return position_data

def write_output(position_data, output_file):
    logging.info(f"Writing output to file: {output_file}")
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
    try:
        with open(output_file, 'w', newline='') as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for data in position_data:
                writer.writerow(data)
        logging.info("Output file written successfully.")
    except Exception as e:
        logging.error(f"Error writing output file: {e}")
        sys.exit(1)

def main():
    args = parse_arguments()
    setup_logging(debug=args.debug)

    logging.info("=== Mutation Analysis Script Started ===")

    # Step 1: Parse the Alignment File
    sequences = parse_alignment(args.alignment)
    sample_ids = list(sequences.keys())
    logging.debug(f"Total samples parsed from alignment: {len(sample_ids)}")

    # Step 2: Translate Nucleotide Sequences to Amino Acid Sequences
    amino_acid_sequences = translate_sequences(sequences)
    logging.debug(f"Translated amino acid sequences for {len(amino_acid_sequences)} samples.")

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

    logging.info("=== Mutation Analysis Script Completed Successfully ===")

if __name__ == '__main__':
    main()
