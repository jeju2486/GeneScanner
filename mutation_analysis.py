import argparse
import logging
import sys
import subprocess
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
import pandas as pd
import xlsxwriter

def parse_arguments():
    parser = argparse.ArgumentParser(description='Mutation Analysis Script')
    parser.add_argument(
        '-a', '--alignment', required=True,
        help='FASTA alignment file'
    )
    parser.add_argument(
        '-t', '--type', required=True, choices=['p', 'n', 'both'],
        help='Alignment type: "p"=protein, "n"=nucleotide, "both"=both'
    )
    parser.add_argument(
        '-f', '--frame', default=1, type=int,
        help='Reading frame for nucleotide alignment'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file name'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Enable debug logging (also writes logs to output.log)'
    )
    # Optional arguments for SNP-sites/MAFFT
    parser.add_argument(
        '--vcf', action='store_true',
        help='Generate VCF file from the input nucleotide alignment using snp-sites'
    )
    parser.add_argument(
        '--snp_sites_path',
        default='snp-sites',
        help='Path or command for snp-sites (default: snp-sites in PATH)'
    )
    parser.add_argument(
        '--mafft_path',
        default='mafft',
        help='Path or command for mafft (default: mafft in PATH)'
    )
    parser.add_argument(
        '--temp', 
        action='store_true',
        help='Keep temporary files generated during the analysis'
    )
    parser.add_argument(
        '--job_id', default='output',
        help='Prefix for output files to avoid clashes in multiple runs'
    )
    return parser.parse_args()


def setup_logging(debug=False):
    """
    Sets up logging to both stdout and to 'output.log' if debug is True.
    """
    level = logging.DEBUG if debug else logging.INFO

    # Always log to console
    handlers = [logging.StreamHandler(sys.stdout)]

    # If debug is enabled, also log to a file named 'output.log'
    if debug:
        file_handler = logging.FileHandler('output.log')
        file_handler.setLevel(level)
        handlers.append(file_handler)

    logging.basicConfig(
        level=level,
        format='%(message)s',
        handlers=handlers
    )


def run_snp_sites(alignment_file, snp_sites_path='snp-sites', output_prefix='out'):
    """
    Runs snp-sites on the provided alignment file to generate a VCF file.
    """
    vcf_file = f'{output_prefix}.vcf'
    cmd = [snp_sites_path, '-v', alignment_file, '-o', vcf_file]
    logging.info(f'Running SNP-sites command: {" ".join(cmd)}')
    try:
        subprocess.run(cmd, check=True)
        logging.info(f'VCF file generated: {vcf_file}')
    except subprocess.CalledProcessError as e:
        logging.error(f'Error running snp-sites: {e}')
        sys.exit(1)
    return vcf_file


def run_mafft(input_fasta, mafft_path='mafft', output_fasta='aligned_proteins.fasta'):
    """
    Runs MAFFT on the provided amino acid FASTA to get a realigned protein alignment.
    MAFFT's log output (stderr) is silenced.
    """
    cmd = [mafft_path, '--auto', '--anysymbol', '--leavegappyregion', input_fasta]
    logging.info(f'Running MAFFT command: {" ".join(cmd)}')
    with open(output_fasta, 'w') as out_f:
        try:
            subprocess.run(cmd, stdout=out_f, stderr=subprocess.DEVNULL, check=True)
            logging.info(f'MAFFT alignment finished: {output_fasta}')
        except subprocess.CalledProcessError as e:
            logging.error(f'Error running MAFFT: {e}')
            sys.exit(1)
    return output_fasta

def parse_alignment(alignment_file, aln_type_str):
    logging.info(f'Parsing {aln_type_str} alignment from: {alignment_file}')
    sequences = {}
    for record in SeqIO.parse(alignment_file, 'fasta'):
        seq = str(record.seq).upper()
        sequences[record.id] = seq
        logging.debug(f'Loaded {record.id} (length={len(seq)})')

    # Check alignment length consistency
    lengths = {len(s) for s in sequences.values()}
    if len(lengths) > 1:
        sys.exit('Error: Sequences are not aligned (unequal lengths).')

    aln_length = lengths.pop()
    logging.info(f'Alignment length: {aln_length}')
    return sequences, aln_length


def translate_codon(codon):
    try:
        return str(Seq(codon).translate())
    except TranslationError:
        return ''


def analyze_nucleotide_alignment(sequences, aln_length, frame=1):
    logging.info('Analyzing nucleotide alignment (syn/non-syn)...')
    data = []
    samples = list(sequences.keys())
    ref_seq = sequences[samples[0]]
    total_samples = len(samples)

    # Adjust to reading frame
    start_pos = frame - 1
    usable_length = aln_length - start_pos
    codon_count = usable_length // 3

    # Prepare mutation matrix
    mut_matrix = pd.DataFrame('', index=samples, columns=range(1, aln_length + 1))

    for pos in range(start_pos, aln_length):
        aln_pos = pos + 1
        wt_nt = ref_seq[pos]
        codon_index = (pos - start_pos) // 3
        codon_pos_in_codon = (pos - start_pos) % 3
        codon_number = codon_index + 1

        ref_codon_raw = ref_seq[start_pos + codon_index * 3 : start_pos + codon_index * 3 + 3]
        ref_codon = ref_codon_raw.replace('-', '')
        ref_aa = translate_codon(ref_codon) if len(ref_codon) == 3 else ''

        syn_count = 0
        nonsyn_count = 0
        ins_count = 0
        del_count = 0
        stop_count = 0
        mutations = {}
        # Track mutated AAs in a set
        mutated_aa_set = set()

        for sid in samples:
            sample_nt = sequences[sid][pos]
            if sample_nt == wt_nt:
                continue

            if wt_nt == '-' and sample_nt != '-':
                mtype = 'Insertion'
                ins_count += 1
            elif wt_nt != '-' and sample_nt == '-':
                mtype = 'Deletion'
                del_count += 1
            else:
                sample_codon_raw = sequences[sid][
                    start_pos + codon_index * 3 : start_pos + codon_index * 3 + 3
                ]
                sample_codon = sample_codon_raw.replace('-', '')
                sample_aa = translate_codon(sample_codon) if len(sample_codon) == 3 else ''

                if sample_aa == '*':
                    mtype = 'Stop Codon'
                    stop_count += 1
                    mutated_aa_set.add('*')
                elif sample_aa == ref_aa and sample_aa != '':
                    mtype = 'Synonymous'
                    syn_count += 1
                elif sample_aa != ref_aa and sample_aa != '' and ref_aa != '':
                    mtype = 'Non-synonymous'
                    nonsyn_count += 1
                    mutated_aa_set.add(sample_aa)
                else:
                    # e.g., partial codon or unknown scenario
                    mtype = 'Unknown'
                    if sample_aa and sample_aa != ref_aa:
                        mutated_aa_set.add(sample_aa)

            mutations[sample_nt] = mtype
            mut_matrix.at[sid, aln_pos] = sample_nt

        indel_stop = ins_count + del_count + stop_count
        total_mut = syn_count + nonsyn_count + indel_stop
        mut_rate = total_mut / total_samples if total_samples else 0

        # Reorder keys so 'Mutated AA(s)' is right after 'Mutations'
        data.append({
            'Alignment Position': aln_pos,
            'Codon Number': codon_number,
            'Codon Position in Codon': codon_pos_in_codon + 1,
            'Wildtype NT': wt_nt,
            'Ref Codon': ref_codon_raw if len(ref_codon_raw) == 3 else '',
            'Ref AA': ref_aa,
            'Mutations': ','.join(sorted(mutations.keys())),
            'Mutated AA(s)': ','.join(sorted(mutated_aa_set)), 
            'Mutation Types': ','.join([mutations[nt] for nt in sorted(mutations.keys())]),
            'Synonymous Mutations': syn_count,
            'Non-synonymous Mutations': nonsyn_count,
            'Insertions': ins_count,
            'Deletions': del_count,
            'Stop Codons': stop_count,
            'Indels/Stop Codon Mutations': indel_stop,
            'Total Mutations': total_mut,
            'Mutation Rate': mut_rate,
        })

    return data, mut_matrix

def remove_reference_fill(protein_seqs):
    new_seqs = {}
    for sid, prot in protein_seqs.items():
        if '*' in prot:
            stop_index = prot.index('*')
            new_prot = prot[:stop_index+1] + '-' * (len(prot) - (stop_index+1))
        else:
            new_prot = prot
        new_seqs[sid] = new_prot
    return new_seqs

def analyze_protein_alignment(sequences, aln_length, effective_aln_lengths=None):
    logging.info('Analyzing protein alignment with real amino acid counts (stop codons end survival)...')

    # Precompute earliest stop position for each sequence (0-indexed).
    # If no stop codon, store None.
    stop_positions = {}
    for sid, seq in sequences.items():
        idx = seq.find('*')
        stop_positions[sid] = idx if idx != -1 else None
    
    data = []
    samples = list(sequences.keys())

    # Mutation matrix
    mut_matrix = pd.DataFrame('', index=samples, columns=range(1, aln_length + 1))
    
    for pos in range(aln_length):
        aln_pos = pos + 1
        # Wildtype amino acid = reference's amino acid at this position
        # (Assume the first sample in the list is reference.)
        ref_seq_id = samples[0]
        wt_aa = sequences[ref_seq_id][pos]

        ins_count = 0
        del_count = 0
        stop_count = 0
        sub_count = 0
        mutations = {}
        mutated_aa_set = set()

        # Count how many sequences are still "alive" at this position.
        survived_count = 0

        for sid in samples:
            # If effective_aln_lengths is provided, skip sequences that
            # are shorter than this position.
            if effective_aln_lengths is not None:
                if pos >= effective_aln_lengths[sid]:
                    continue

            # If this sequence encountered a stop codon previously,
            # do not count it for survival at or after that position.
            if stop_positions[sid] is not None and pos >= stop_positions[sid]:
                continue

            survived_count += 1
            sample_aa = sequences[sid][pos]

            # No mutation if it's identical to the wildtype
            if sample_aa == wt_aa:
                continue
            
            # Classify the mutation type
            if wt_aa == '-' and sample_aa != '-':
                mtype = 'Insertion'
                ins_count += 1
                if sample_aa not in ['-', '*']:
                    mutated_aa_set.add(sample_aa)
            elif wt_aa != '-' and sample_aa == '-':
                mtype = 'Deletion'
                del_count += 1
            else:
                if sample_aa == '*':
                    mtype = 'Stop Codon'
                    stop_count += 1
                    mutated_aa_set.add('*')
                else:
                    mtype = 'Substitution'
                    sub_count += 1
                    if sample_aa != '-':
                        mutated_aa_set.add(sample_aa)

            mutations[sample_aa] = mtype
            mut_matrix.at[sid, aln_pos] = sample_aa

        total_mut = ins_count + del_count + stop_count + sub_count
        mut_rate = (total_mut / survived_count) if survived_count > 0 else 0
        
        data.append({
            'Alignment Position': aln_pos,
            'Wildtype AA': wt_aa,
            'Mutations': ','.join(sorted(mutations.keys())),
            'Mutated AA(s)': ','.join(sorted(mutated_aa_set)),
            'Mutation Types': ','.join([mutations[aa] for aa in sorted(mutations.keys())]),
            'Insertions': ins_count,
            'Deletions': del_count,
            'Stop Codons': stop_count,
            'Substitutions': sub_count,
            'Total Mutations': total_mut,
            'Mutation Rate': mut_rate,
            'Survived Count': survived_count
        })
    
    return data, mut_matrix


def translate_nucleotide_to_protein(sequences, frame=1):
    logging.info('Translating nucleotide sequences (ungapped) to protein and replacing truncated area with reference...')
    protein_seqs = {}
    effective_lengths = {}

    # Use the first sequence as the reference.
    ref_id = next(iter(sequences))
    ref_nt_seq = sequences[ref_id].replace('-', '')
    ref_adjusted_seq = ref_nt_seq[frame - 1:]
    remainder = len(ref_adjusted_seq) % 3
    if remainder != 0:
        ref_adjusted_seq += 'N' * (3 - remainder)
    ref_full_prot = str(Seq(ref_adjusted_seq).translate())

    # Process each sequence.
    for sid, nt_seq in sequences.items():
        nt_nogap = nt_seq.replace('-', '')
        adjusted_seq = nt_nogap[frame - 1:]
        remainder = len(adjusted_seq) % 3
        if remainder != 0:
            adjusted_seq += 'N' * (3 - remainder)
        full_prot = str(Seq(adjusted_seq).translate())

        # If a stop codon is found, keep everything up to and including the first stop
        # and then fill in the remainder from the reference protein.
        if '*' in full_prot:
            # Get the translation up to (and including) the first stop.
            prefix = full_prot.split('*')[0] + '*'
            # Fill with the corresponding region from the reference, if available.
            fill = ref_full_prot[len(prefix):] if len(ref_full_prot) > len(prefix) else ''
            prot = prefix + fill
        else:
            prot = full_prot

        protein_seqs[sid] = prot
        effective_lengths[sid] = len(prot)

    return protein_seqs, effective_lengths


def format_excel_columns(df, worksheet):
    for i, col in enumerate(df.columns):
        max_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
        worksheet.set_column(i, i, max_len)


def write_matrix_sheet(writer, df, sheet_name):
    df.to_excel(writer, sheet_name=sheet_name)
    ws = writer.sheets[sheet_name]
    # Format headers
    for col_num, value in enumerate(df.columns.values):
        ws.write(0, col_num + 1, value)  # +1 offset for index column
    ws.set_column(0, 0, max(df.index.astype(str).map(len).max(), 12))
    for i, col in enumerate(df.columns):
        max_len = max(df[col].astype(str).map(len).max(), len(str(col))) + 2
        ws.set_column(i + 1, i + 1, max_len)

def compute_effective_alignment_lengths(aligned_protein_seqs, original_effective_lengths):
    effective_aln_lengths = {}
    for sid, seq in aligned_protein_seqs.items():
        count = 0
        orig_len = original_effective_lengths[sid]
        for i, char in enumerate(seq):
            if char != '-':
                count += 1
            if count == orig_len:
                effective_aln_lengths[sid] = i + 1  # using 1-indexing for alignment positions
                break
        if sid not in effective_aln_lengths:
            effective_aln_lengths[sid] = len(seq)
    return effective_aln_lengths


def write_analysis_sheet(writer, data, sheet_name, ref_seq):
    df = pd.DataFrame(data)
    if df.empty:
        return None
    df.to_excel(writer, sheet_name=sheet_name, index=False, startrow=2)
    ws = writer.sheets[sheet_name]
    ws.write(0, 0, f'Reference sequence: {ref_seq}')
    for col_num, value in enumerate(df.columns.values):
        ws.write(2, col_num, value)
    format_excel_columns(df, ws)
    return df


def summarize_data(df, analysis_type):
    total_pos = len(df)
    mutated = df[df['Total Mutations'] > 0]
    num_mut = len(mutated)
    perc_mut = (num_mut / total_pos * 100) if total_pos else 0

    high_mut = mutated[mutated['Mutation Rate'] > 0.2]
    num_high = len(high_mut)
    perc_high = (num_high / total_pos * 100) if total_pos else 0

    summary = {
        'Analysis Type': analysis_type.capitalize(),
        'Total Positions': total_pos,
        'Positions with Mutations': num_mut,
        'Percent Positions with Mutations': f'{perc_mut:.2f}',
        'Positions with >20% Mutations': num_high,
        'Percent Positions with >20% Mutations': f'{perc_high:.2f}'
    }

    if analysis_type == 'n':
        for c in [
            'Synonymous Mutations', 'Non-synonymous Mutations',
            'Insertions', 'Deletions', 'Stop Codons', 'Total Mutations'
        ]:
            df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)
        summary.update({
            'Total Synonymous Mutations': df['Synonymous Mutations'].sum(),
            'Total Non-synonymous Mutations': df['Non-synonymous Mutations'].sum(),
            'Total Indels': df[['Insertions', 'Deletions']].sum().sum(),
            'Total Stop Codons': df['Stop Codons'].sum(),
            'Total Mutations': df['Total Mutations'].sum()
        })
    else:
        for c in [
            'Insertions', 'Deletions', 'Stop Codons', 'Substitutions', 'Total Mutations'
        ]:
            df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)
        summary.update({
            'Total Insertions': df['Insertions'].sum(),
            'Total Deletions': df['Deletions'].sum(),
            'Total Stop Codons': df['Stop Codons'].sum(),
            'Total Substitutions': df['Substitutions'].sum(),
            'Total Mutations': df['Total Mutations'].sum()
        })

    return summary

def write_fasta(sequences, output_file):
    with open(output_file, 'w') as f:
        for sid, seq in sequences.items():
            f.write(f'>{sid}\n{seq}\n')

def write_summary_sheet(writer, summaries):
    if not summaries:
        return
    df = pd.DataFrame(summaries)
    df.to_excel(writer, sheet_name='Summary Statistics', index=False)
    ws = writer.sheets['Summary Statistics']
    for col_num, value in enumerate(df.columns.values):
        ws.write(0, col_num, value)
    for i, col in enumerate(df.columns):
        max_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
        ws.set_column(i, i, max_len)


def write_output(output_file,
                 nuc_data, prot_data,
                 sequences, prot_sequences,
                 nuc_matrix, prot_matrix):
    logging.info(f'Writing output to {output_file}')
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        summaries = []
        # Nucleotide sheets
        if nuc_data:
            df_nuc = write_analysis_sheet(
                writer, nuc_data, 'Nucleotide Analysis',
                sequences[list(sequences.keys())[0]]
            )
            if df_nuc is not None and nuc_matrix is not None:
                write_matrix_sheet(writer, nuc_matrix, 'Nucleotide Mutation Matrix')
                summaries.append(summarize_data(df_nuc, 'n'))

        # Protein sheets
        if prot_data:
            df_prot = write_analysis_sheet(
                writer, prot_data, 'Protein Analysis',
                prot_sequences[list(prot_sequences.keys())[0]] if prot_sequences else ''
            )
            if df_prot is not None and prot_matrix is not None:
                write_matrix_sheet(writer, prot_matrix, 'Protein Mutation Matrix')
                summaries.append(summarize_data(df_prot, 'p'))

        # Summary
        write_summary_sheet(writer, summaries)

    logging.info('Output file written successfully.')


def main():
    args = parse_arguments()
    setup_logging(debug=args.debug)
    logging.info('Starting Mutation Analysis')

    aln_type = args.type.lower()
    if aln_type == 'n':
        aln_type_full = 'nucleotide'
    elif aln_type == 'p':
        aln_type_full = 'protein'
    elif aln_type == 'both':
        aln_type_full = 'both'
    else:
        sys.exit('Invalid alignment type.')

    sequences, aln_length = parse_alignment(args.alignment, aln_type_full)

    # If requested, generate VCF for NT alignments
    if args.vcf and (aln_type in ['n', 'both']):
        prefix = os.path.splitext(os.path.basename(args.alignment))[0]
        run_snp_sites(args.alignment, snp_sites_path=args.snp_sites_path, output_prefix=prefix)

    nuc_data = None
    prot_data = None
    nuc_matrix = None
    prot_matrix = None
    prot_sequences = None

    if aln_type == 'n':
        # NT analysis only
        nuc_data, nuc_matrix = analyze_nucleotide_alignment(
            sequences, aln_length, frame=args.frame
        )
    elif aln_type == 'p':
        # Protein analysis only (assumes the input file is already a protein alignment)
        prot_data, prot_matrix = analyze_protein_alignment(sequences, aln_length)
    elif aln_type == 'both':
        # 1) Nucleotide analysis on the current alignment
        nuc_data, nuc_matrix = analyze_nucleotide_alignment(
            sequences, aln_length, frame=args.frame
        )

        # 2) Translate original NT -> realign in protein space -> analyze protein
        logging.info("Resetting to original ungapped NT -> translating -> re-aligning with MAFFT...")

        # The new translate function returns both the protein sequences (truncated at the first stop codon)
        # and a dictionary of the original effective lengths (the count of amino acids before the stop).
        prot_sequences, original_effective_lengths = translate_nucleotide_to_protein(sequences, frame=args.frame)
        
        # Write the translated protein sequences to a temporary FASTA file
        tmp_prot_in = f'{args.job_id}_proteins_in.tmp'
        with open(tmp_prot_in, 'w') as f:
            for sid, prot_seq in prot_sequences.items():
                f.write(f'>{sid}\n{prot_seq}\n')

        # Run MAFFT to realign the protein sequences
        # Run MAFFT to realign the translated protein sequences.
        aligned_prot_fasta = f'{args.job_id}_proteins_aligned.tmp'
        run_mafft(tmp_prot_in, mafft_path=args.mafft_path, output_fasta=aligned_prot_fasta)

        # Parse the newly aligned protein sequences.
        prot_seqs_aligned, prot_aln_length = parse_alignment(aligned_prot_fasta, 'protein')

        # Remove the reference fill after the stop codon (replace characters after '*' with dashes).
        prot_seqs_aligned = remove_reference_fill(prot_seqs_aligned)

        # Overwrite the temporary FASTA file with the modified sequences.
        write_fasta(prot_seqs_aligned, aligned_prot_fasta)

        # Continue with your downstream processing...
        effective_aln_lengths = compute_effective_alignment_lengths(prot_seqs_aligned, original_effective_lengths)
        prot_data, prot_matrix = analyze_protein_alignment(prot_seqs_aligned, prot_aln_length,
                                                            effective_aln_lengths=effective_aln_lengths)
        prot_sequences = prot_seqs_aligned

        if not args.temp:
            os.remove(tmp_prot_in)
            os.remove(aligned_prot_fasta)
            
    # Write final output (Excel workbook)
    output_file = args.output if args.output else f'{args.job_id}_{aln_type_full}_mutation_analysis.xlsx'
    write_output(
        output_file,
        nuc_data, prot_data,
        sequences, prot_sequences or {},
        nuc_matrix, prot_matrix
    )

    logging.info('Mutation Analysis Completed Successfully')

if __name__ == '__main__':
    main()
