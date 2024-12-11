import argparse
import logging
import sys

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
        help='Enable debug logging'
    )
    return parser.parse_args()


def setup_logging(debug=False):
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )


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

        ref_codon_raw = ref_seq[
            start_pos + codon_index * 3:start_pos + codon_index * 3 + 3
        ]
        ref_codon = ref_codon_raw.replace('-', '')
        ref_aa = translate_codon(ref_codon) if len(ref_codon) == 3 else ''

        syn_count = 0
        nonsyn_count = 0
        ins_count = 0
        del_count = 0
        stop_count = 0
        mutations = {}

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
                    start_pos + codon_index * 3:start_pos + codon_index * 3 + 3
                ]
                sample_codon = sample_codon_raw.replace('-', '')
                sample_aa = translate_codon(sample_codon) if len(sample_codon) == 3 else ''

                if sample_aa == '*':
                    mtype = 'Stop Codon'
                    stop_count += 1
                elif sample_aa == ref_aa and sample_aa != '':
                    mtype = 'Synonymous'
                    syn_count += 1
                elif (sample_aa != ref_aa and sample_aa != '' and
                        ref_aa != ''):
                    mtype = 'Non-synonymous'
                    nonsyn_count += 1
                else:
                    mtype = 'Unknown'

            mutations[sample_nt] = mtype
            mut_matrix.at[sid, aln_pos] = sample_nt

        indel_stop = ins_count + del_count + stop_count
        total_mut = syn_count + nonsyn_count + indel_stop
        mut_rate = total_mut / total_samples if total_samples else 0

        data.append({
            'Alignment Position': aln_pos,
            'Codon Number': codon_number,
            'Codon Position in Codon': codon_pos_in_codon + 1,
            'Wildtype NT': wt_nt,
            'Ref Codon': ref_codon_raw if len(ref_codon_raw) == 3 else '',
            'Ref AA': ref_aa,
            'Mutations': ','.join(sorted(mutations.keys())),
            'Mutation Types': ','.join(
                [mutations[nt] for nt in sorted(mutations.keys())]
            ),
            'Synonymous Mutations': syn_count,
            'Non-synonymous Mutations': nonsyn_count,
            'Insertions': ins_count,
            'Deletions': del_count,
            'Stop Codons': stop_count,
            'Indels/Stop Codon Mutations': indel_stop,
            'Total Mutations': total_mut,
            'Mutation Rate': mut_rate,
            'Mutation Percentage': mut_rate * 100
        })

    return data, mut_matrix


def analyze_protein_alignment(sequences, aln_length):
    logging.info('Analyzing protein alignment...')
    data = []
    samples = list(sequences.keys())
    ref_seq = sequences[samples[0]]
    total_samples = len(samples)

    mut_matrix = pd.DataFrame('', index=samples, columns=range(1, aln_length + 1))

    for pos in range(aln_length):
        aln_pos = pos + 1
        wt_aa = ref_seq[pos]

        ins_count = 0
        del_count = 0
        stop_count = 0
        sub_count = 0
        mutations = {}

        for sid in samples:
            sample_aa = sequences[sid][pos] if pos < len(sequences[sid]) else '-'
            if sample_aa == wt_aa:
                continue

            if wt_aa == '-' and sample_aa != '-':
                mtype = 'Insertion'
                ins_count += 1
            elif wt_aa != '-' and sample_aa == '-':
                mtype = 'Deletion'
                del_count += 1
            else:
                if sample_aa == '*':
                    mtype = 'Stop Codon'
                    stop_count += 1
                else:
                    mtype = 'Substitution'
                    sub_count += 1

            mutations[sample_aa] = mtype
            mut_matrix.at[sid, aln_pos] = sample_aa

        total_mut = ins_count + del_count + stop_count + sub_count
        mut_rate = total_mut / total_samples if total_samples else 0

        data.append({
            'Alignment Position': aln_pos,
            'Wildtype AA': wt_aa,
            'Mutations': ','.join(sorted(mutations.keys())),
            'Mutation Types': ','.join(
                [mutations[aa] for aa in sorted(mutations.keys())]
            ),
            'Insertions': ins_count,
            'Deletions': del_count,
            'Stop Codons': stop_count,
            'Substitutions': sub_count,
            'Total Mutations': total_mut,
            'Mutation Rate': mut_rate,
            'Mutation Percentage': mut_rate * 100
        })

    return data, mut_matrix


def translate_nucleotide_to_protein(sequences, frame=1):
    logging.info('Translating nucleotide sequences to protein...')
    protein_seqs = {}
    for sid, nt_seq in sequences.items():
        nt_nogap = nt_seq.replace('-', '')
        protein = str(Seq(nt_nogap[frame - 1:]).translate())
        protein_seqs[sid] = protein
    return protein_seqs


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
    # Generic summary for both nucleotide and protein
    total_pos = len(df)
    mutated = df[df['Total Mutations'] > 0]
    num_mut = len(mutated)
    perc_mut = (num_mut / total_pos * 100) if total_pos else 0

    high_mut = mutated[mutated['Mutation Percentage'] > 20]
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

    # Additional analysis-type-specific metrics
    if analysis_type == 'n':
        # Convert numeric columns for easy sum
        for c in [
            'Synonymous Mutations', 'Non-synonymous Mutations',
            'Insertions', 'Deletions', 'Stop Codons', 'Total Mutations'
        ]:
            df[c] = pd.to_numeric(df[c])
        summary.update({
            'Total Synonymous Mutations': df['Synonymous Mutations'].sum(),
            'Total Non-synonymous Mutations': df['Non-synonymous Mutations'].sum(),
            'Total Indels': df[['Insertions', 'Deletions']].sum().sum(),
            'Total Stop Codons': df['Stop Codons'].sum(),
            'Total Mutations': df['Total Mutations'].sum()
        })
    else:  # protein
        for c in [
            'Insertions', 'Deletions', 'Stop Codons', 'Substitutions', 'Total Mutations'
        ]:
            df[c] = pd.to_numeric(df[c])
        summary.update({
            'Total Insertions': df['Insertions'].sum(),
            'Total Deletions': df['Deletions'].sum(),
            'Total Stop Codons': df['Stop Codons'].sum(),
            'Total Substitutions': df['Substitutions'].sum(),
            'Total Mutations': df['Total Mutations'].sum()
        })

    return summary


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


def write_output(output_file, nuc_data, prot_data, sequences, prot_sequences, nuc_matrix, prot_matrix):
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
                prot_sequences[list(prot_sequences.keys())[0]]
                if prot_sequences else ''
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

    nuc_data = None
    prot_data = None
    nuc_matrix = None
    prot_matrix = None
    prot_sequences = None

    # Perform analyses
    if aln_type == 'n':
        nuc_data, nuc_matrix = analyze_nucleotide_alignment(
            sequences, aln_length, frame=args.frame
        )
    elif aln_type == 'p':
        prot_data, prot_matrix = analyze_protein_alignment(sequences, aln_length)
    else:  # both
        # Nucleotide analysis
        nuc_data, nuc_matrix = analyze_nucleotide_alignment(
            sequences, aln_length, frame=args.frame
        )
        # Translate and do protein analysis
        prot_sequences = translate_nucleotide_to_protein(sequences, frame=args.frame)
        prot_data, prot_matrix = analyze_protein_alignment(
            prot_sequences,
            len(prot_sequences[list(prot_sequences.keys())[0]])
        )

    # Default output if none provided
    output_file = args.output if args.output else f'{aln_type_full}_mutation_analysis.xlsx'
    write_output(output_file, nuc_data, prot_data, sequences, prot_sequences or {},
                 nuc_matrix, prot_matrix)

    logging.info('Mutation Analysis Completed Successfully')


if __name__ == '__main__':
    main()
