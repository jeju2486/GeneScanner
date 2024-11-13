#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import sys
import logging

def parse_arguments():
    parser = argparse.ArgumentParser(description='Mutation Analysis Script')
    parser.add_argument('-a', '--alignment', required=True, help='Path to the alignment file in FASTA format')
    parser.add_argument('-t', '--type', required=True, choices=['p', 'n'], help='Alignment type: "p" for protein or "n" for nucleotide')
    parser.add_argument('-o', '--output', help='Output file name (default based on input file and type)')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug mode for detailed logging')
    args = parser.parse_args()
    return args

def setup_logging(debug=False):
    """Set up the logging configuration."""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

def parse_alignment(alignment_file, aln_type):
    logging.info(f"Parsing {aln_type} alignment file: {alignment_file}")
    sequences = {}
    try:
        for record in SeqIO.parse(alignment_file, 'fasta'):
            seq = str(record.seq).upper()
            sequences[record.id] = seq
            logging.debug(f"Loaded sequence for sample '{record.id}' with length {len(seq)}")
    except Exception as e:
        logging.error(f"Error parsing alignment file: {e}")
        sys.exit(1)
    
    # Verify that all sequences are of the same length
    seq_lengths = set(len(seq) for seq in sequences.values())
    if len(seq_lengths) > 1:
        logging.error("Sequences are not aligned. All sequences must be of the same length.")
        sys.exit("Error: Sequences are not aligned. All sequences must be of the same length.")
    alignment_length = seq_lengths.pop()
    logging.info(f"All sequences are aligned with length {alignment_length} positions.")
    return sequences, alignment_length

def analyze_nucleotide_alignment(sequences, alignment_length):
    logging.info("Starting nucleotide mutation analysis.")
    position_data = []

    sample_ids = list(sequences.keys())
    total_samples = len(sample_ids)

    # Use the first sequence as the reference
    reference_sample = sample_ids[0]
    reference_seq = sequences[reference_sample]
    logging.debug(f"Reference sample: '{reference_sample}'")

    for pos in range(alignment_length):
        alignment_position = pos + 1
        wildtype_nt = reference_seq[pos]

        # Skip positions with gaps in the reference
        if wildtype_nt == '-':
            logging.debug(f"Gap in reference at position {alignment_position}, skipping.")
            continue

        mutations = set()
        for sample_id in sample_ids:
            sample_nt = sequences[sample_id][pos]
            if sample_nt != wildtype_nt and sample_nt != '-':
                mutations.add(sample_nt)

        mutation_count = len(mutations)
        mutation_rate = mutation_count / total_samples if total_samples else 0

        position_data.append({
            'alignment_position': alignment_position,
            'wildtype': wildtype_nt,
            'mutations': ','.join(sorted(mutations)),
            'mutation_count': mutation_count,
            'mutation_rate': f"{mutation_rate:.6f}"
        })

    logging.info("Nucleotide mutation analysis completed.")
    return position_data

def analyze_protein_alignment(sequences, alignment_length):
    logging.info("Starting protein mutation analysis.")
    position_data = []

    sample_ids = list(sequences.keys())
    total_samples = len(sample_ids)

    # Use the first sequence as the reference
    reference_sample = sample_ids[0]
    reference_seq = sequences[reference_sample]
    logging.debug(f"Reference sample: '{reference_sample}'")

    for pos in range(alignment_length):
        alignment_position = pos + 1
        wildtype_aa = reference_seq[pos]

        # Skip positions with gaps in the reference
        if wildtype_aa == '-':
            logging.debug(f"Gap in reference at position {alignment_position}, skipping.")
            continue

        mutations = set()
        nonsense_count = 0
        for sample_id in sample_ids:
            sample_aa = sequences[sample_id][pos]
            if sample_aa != wildtype_aa and sample_aa != '-':
                mutations.add(sample_aa)
                if sample_aa == '*':
                    nonsense_count += 1

        mutation_count = len(mutations)
        mutation_rate = mutation_count / total_samples if total_samples else 0
        nonsense_rate = nonsense_count / total_samples if total_samples else 0

        position_data.append({
            'alignment_position': alignment_position,
            'wildtype': wildtype_aa,
            'mutations': ','.join(sorted(mutations)),
            'mutation_count': mutation_count,
            'mutation_rate': f"{mutation_rate:.6f}",
            'nonsense_count': nonsense_count,
            'nonsense_rate': f"{nonsense_rate:.6f}"
        })

    logging.info("Protein mutation analysis completed.")
    return position_data

def write_output(position_data, output_file, aln_type):
    logging.info(f"Writing output to file: {output_file}")

    if aln_type == 'n':
        fieldnames = ['alignment_position', 'wildtype', 'mutations', 'mutation_count', 'mutation_rate']
    elif aln_type == 'p':
        fieldnames = ['alignment_position', 'wildtype', 'mutations', 'mutation_count', 'mutation_rate', 'nonsense_count', 'nonsense_rate']

    try:
        with open(output_file, 'w', newline='') as tsvfile:
            header = '\t'.join(fieldnames)
            tsvfile.write(header + '\n')
            for data in position_data:
                row = '\t'.join(str(data[field]) for field in fieldnames)
                tsvfile.write(row + '\n')
        logging.info("Output file written successfully.")
    except Exception as e:
        logging.error(f"Error writing output file: {e}")
        sys.exit(1)

def main():
    args = parse_arguments()
    setup_logging(debug=args.debug)

    logging.info("=== Mutation Analysis Script Started ===")

    aln_type = args.type.lower()
    if aln_type == 'n':
        aln_type_full = 'nucleotide'
    elif aln_type == 'p':
        aln_type_full = 'protein'
    else:
        logging.error("Invalid alignment type. Use 'p' for protein or 'n' for nucleotide.")
        sys.exit(1)

    # Parse the alignment file
    sequences, alignment_length = parse_alignment(args.alignment, aln_type_full)

    # Analyze mutations based on alignment type
    if aln_type == 'n':
        position_data = analyze_nucleotide_alignment(sequences, alignment_length)
    elif aln_type == 'p':
        position_data = analyze_protein_alignment(sequences, alignment_length)

    # Set default output file name if not provided
    if args.output:
        output_file = args.output
    else:
        output_file = f"{aln_type_full}_mutation_analysis.tsv"

    # Write output
    write_output(position_data, output_file, aln_type)

    logging.info("=== Mutation Analysis Script Completed Successfully ===")

if __name__ == '__main__':
    main()
