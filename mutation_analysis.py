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

###############################################################################
#                          UTILITY & PARSING FUNCTIONS
###############################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description='Mutation Analysis Script')
    parser.add_argument('-a', '--alignment', required=True,
                        help='FASTA alignment file')
    parser.add_argument('-t', '--type', required=True, choices=['p', 'n', 'both'],
                        help='Alignment type: "p"=protein, "n"=nucleotide, "both"=both')
    parser.add_argument('-f', '--frame', default=1, type=int,
                        help='Reading frame for nucleotide alignment. Default=1')
    parser.add_argument('-o', '--output',
                        help='Output Excel file name')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Enable debug logging (writes logs to output.log)')

    # Optional arguments for SNP-sites/MAFFT
    parser.add_argument('--vcf', action='store_true',
                        help='Generate VCF file from the input nucleotide alignment using snp-sites')
    parser.add_argument('--snp_sites_path', default='snp-sites',
                        help='Path or command for snp-sites (default: snp-sites in PATH)')
    parser.add_argument('--mafft_path', default='mafft',
                        help='Path or command for mafft (default: mafft in PATH)')

    # Housekeeping
    parser.add_argument('--temp', action='store_true',
                        help='Keep temporary files generated during the analysis')
    parser.add_argument('--tmp_dir', default='./',
                        type=lambda p: p if p.endswith('/') else p + '/',
                        help='Directory for temporary files (default=./)')
    parser.add_argument('--job_id', default='output',
                        help='Prefix for output files (default=output)')

    # Grouping & References
    parser.add_argument('--groups',
                        help='Optional CSV file with columns [Isolate,Category] for grouping')
    parser.add_argument('--reference',
                        help='Optional reference isolate ID from the FASTA file')

    # NEW: Strict type validation
    parser.add_argument('--strict_validation', action='store_true',
                        help='If set, error out if the input alignment does not match the chosen -t (n/p/both)')

    return parser.parse_args()


def setup_logging(debug=False):
    level = logging.DEBUG if debug else logging.INFO
    handlers = [logging.StreamHandler(sys.stdout)]
    if debug:
        file_handler = logging.FileHandler('output.log')
        file_handler.setLevel(level)
        handlers.append(file_handler)
    logging.basicConfig(level=level, format='%(message)s', handlers=handlers)


def run_snp_sites(alignment_file, snp_sites_path='snp-sites', output_prefix='out'):
    """Runs snp-sites on the provided alignment file to generate a VCF file."""
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

###############################################################################
#                          TYPE VALIDATION
###############################################################################

def guess_sequence_type(sequences):
    """
    A simple heuristic to guess if the alignment is mostly nucleotide or protein.
    If 85%+ of the non-gap chars are in A,C,G,T,N, we call it "nucleotide".
    Otherwise, we call it "protein".
    """
    if not sequences:
        return "unknown"

    valid_nuc = set("ACGTN-")
    total_chars = 0
    nuc_chars = 0

    for sid, seq in sequences.items():
        for c in seq:
            if c == '-':
                continue
            total_chars += 1
            if c in valid_nuc:
                nuc_chars += 1

    if total_chars == 0:
        return "unknown"
    frac_nuc = nuc_chars / total_chars
    if frac_nuc >= 0.85:
        return "nucleotide"
    else:
        return "protein"

###############################################################################
#                         MAIN ANALYSIS FUNCTIONS
###############################################################################

def translate_codon(codon):
    try:
        return str(Seq(codon).translate())
    except TranslationError:
        return ''


def analyze_nucleotide_alignment(sequences, aln_length, frame=1, reference_id=None):
    """
    Performs a codon-aware nucleotide analysis, restoring:
      - "Codon Number"
      - "Mutations"
      - "Mutation Types"
      - "Reference NT"
      - "Syn/Non-syn" classification using the entire codon
      - "Mutation Frequency" (instead of "Mutation Rate")
    """

    import pandas as pd
    from Bio.Seq import Seq
    from Bio.Data.CodonTable import TranslationError

    logging.info('Analyzing nucleotide alignment (syn/non-syn) with codon logic...')

    # Reorder sequences so the requested reference is first
    all_ids = list(sequences.keys())
    if reference_id and reference_id in all_ids:
        all_ids.remove(reference_id)
        samples = [reference_id] + all_ids
        actual_ref_id = reference_id
    else:
        samples = all_ids
        actual_ref_id = samples[0] if samples else "NoIsolates"

    if not samples:
        return [], pd.DataFrame(), actual_ref_id

    ref_seq = sequences[samples[0]]
    total_samples = len(samples)

    # set up
    start_pos = frame - 1
    mut_matrix = pd.DataFrame('', index=samples, columns=range(1, aln_length + 1))

    data_rows = []

    for pos in range(start_pos, aln_length):
        aln_pos = pos + 1
        ref_nt = ref_seq[pos]  # "Reference NT"
        
        # Identify which codon we are in, using the reading frame
        codon_index = (pos - start_pos) // 3  # 0-based
        codon_number = codon_index + 1        # 1-based for user

        # We'll extract the full reference codon (3 nts) for translation
        ref_codon_start = start_pos + codon_index * 3
        ref_codon = ref_seq[ref_codon_start : ref_codon_start + 3].replace('-', '')
        try:
            ref_aa = str(Seq(ref_codon).translate()) if len(ref_codon) == 3 else ''
        except TranslationError:
            ref_aa = ''

        syn_count = 0
        nonsyn_count = 0
        ins_count = 0
        del_count = 0
        stop_count = 0

        # We'll store each mutated base with its type (Syn, Non-syn, etc.)
        base2mutation = {}

        for sid in samples:
            sample_nt = sequences[sid][pos]
            if sample_nt == ref_nt:
                continue  # no difference at this position

            # Distinguish insertion, deletion, or partial codon difference
            if ref_nt == '-' and sample_nt != '-':
                mtype = 'Insertion'
                ins_count += 1
            elif ref_nt != '-' and sample_nt == '-':
                mtype = 'Deletion'
                del_count += 1
            else:
                # We do a full codon translation check
                sample_codon = sequences[sid][ref_codon_start : ref_codon_start + 3].replace('-', '')
                try:
                    sample_aa = str(Seq(sample_codon).translate()) if len(sample_codon) == 3 else ''
                except TranslationError:
                    sample_aa = ''

                if sample_aa == '*':
                    mtype = 'Stop Codon'
                    stop_count += 1
                elif sample_aa == ref_aa and sample_aa != '':
                    mtype = 'Synonymous'
                    syn_count += 1
                elif sample_aa != ref_aa and sample_aa != '' and ref_aa != '':
                    mtype = 'Non-synonymous'
                    nonsyn_count += 1
                else:
                    mtype = 'Unknown'

            base2mutation[sample_nt] = mtype
            mut_matrix.at[sid, aln_pos] = sample_nt

        # Summarize for this position
        indel_stop = ins_count + del_count + stop_count
        total_mut = syn_count + nonsyn_count + indel_stop
        mut_freq = (total_mut / total_samples) if total_samples else 0

        row_data = {
            'Alignment Position': aln_pos,
            'Codon Number': codon_number,  # restored
            'Reference NT': ref_nt,
            'Mutations': ','.join(sorted(base2mutation.keys())),
            'Mutation Types': ','.join(base2mutation[b] for b in sorted(base2mutation.keys())),
            'Synonymous Mutations': syn_count,
            'Non-synonymous Mutations': nonsyn_count,
            'Insertions': ins_count,
            'Deletions': del_count,
            'Stop Codons': stop_count,
            'Total Mutations': total_mut,
            'Mutation Frequency': mut_freq
        }
        data_rows.append(row_data)

    # Convert to DataFrame
    data_df = pd.DataFrame(data_rows)
    return data_df.to_dict('records'), mut_matrix, actual_ref_id


def remove_reference_fill(protein_seqs):
    """
    If there's a stop codon '*', replace everything after the first '*' with '-'.
    """
    new_seqs = {}
    for sid, prot in protein_seqs.items():
        if '*' in prot:
            stop_index = prot.index('*')
            new_prot = prot[:stop_index+1] + '-' * (len(prot) - (stop_index+1))
        else:
            new_prot = prot
        new_seqs[sid] = new_prot
    return new_seqs


def analyze_protein_alignment(sequences, aln_length, effective_aln_lengths=None, reference_id=None):
    """
    Returns (list_of_analysis_dicts, mutation_matrix_dataframe, actual_reference_id).
      - rename "Wildtype AA" -> "Reference AA"
      - rename "Mutations" -> "Substitutions"
      - rename "Mutation Types" -> "Substitution Types"
      - rename "Mutation Rate" -> "Mutation Frequency"
      - rename "Survived Count" -> "Remaining sequences before encountering stop codon"
      - reorder "Substitutions" column before "Insertions"
    """
    logging.info('Analyzing protein alignment...')

    all_ids = list(sequences.keys())
    if reference_id and reference_id in all_ids:
        all_ids.remove(reference_id)
        samples = [reference_id] + all_ids
        actual_ref_id = reference_id
    else:
        samples = all_ids
        actual_ref_id = samples[0] if samples else "NoIsolates"

    if not samples:
        return [], pd.DataFrame(), actual_ref_id

    ref_seq_id = samples[0]
    data = []

    stop_positions = {}
    for sid, seq in sequences.items():
        idx = seq.find('*')
        stop_positions[sid] = idx if idx != -1 else None

    mut_matrix = pd.DataFrame('', index=samples, columns=range(1, aln_length + 1))
    
    for pos in range(aln_length):
        aln_pos = pos + 1
        wt_aa = sequences[ref_seq_id][pos]

        ins_count = 0
        del_count = 0
        stop_count = 0
        sub_count = 0
        mutations = {}
        survived_count = 0

        for sid in samples:
            if effective_aln_lengths is not None:
                if pos >= effective_aln_lengths[sid]:
                    continue

            if stop_positions[sid] is not None and pos > stop_positions[sid]:
                # Past the stop codon
                continue

            survived_count += 1
            sample_aa = sequences[sid][pos]

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
        mut_freq = (total_mut / survived_count) if survived_count > 0 else 0

        row_dict = {
            'Alignment Position': aln_pos,
            'Reference AA': wt_aa,
            'Substitutions': ','.join(
                sorted([aa for aa in mutations if mutations[aa] == 'Substitution'])
            ),
            'Substitution Types': ','.join([mutations[aa] for aa in sorted(mutations.keys())]),
            'Substitution Count': sub_count,
            'Insertions': ins_count,
            'Deletions': del_count,
            'Stop Codons': stop_count,
            'Total Mutations': total_mut,
            'Mutation Frequency': mut_freq,
            'Remaining sequences before encountering stop codon': survived_count
        }
        data.append(row_dict)

    return data, mut_matrix, actual_ref_id


def translate_nucleotide_to_protein(sequences, frame=1, reference_id=None):
    """
    Translates each ungapped nucleotide sequence into protein.
    Returns (dict_of_prot_sequences, dict_of_effective_lengths, actual_ref_id).
    """
    logging.info('Translating nucleotide sequences (ungapped) to proteins...')

    all_ids = list(sequences.keys())
    if reference_id and reference_id in all_ids:
        ref_id = reference_id
        actual_ref_id = reference_id
    else:
        ref_id = all_ids[0] if all_ids else None
        actual_ref_id = ref_id if ref_id else "NoIsolates"

    protein_seqs = {}
    effective_lengths = {}
    if ref_id is None:
        return protein_seqs, effective_lengths, actual_ref_id

    ref_nt_seq = sequences[ref_id].replace('-', '')
    ref_adjusted_seq = ref_nt_seq[frame - 1:]
    remainder = len(ref_adjusted_seq) % 3
    if remainder != 0:
        ref_adjusted_seq += 'N' * (3 - remainder)
    ref_full_prot = str(Seq(ref_adjusted_seq).translate())

    for sid, nt_seq in sequences.items():
        nt_nogap = nt_seq.replace('-', '')
        adjusted_seq = nt_nogap[frame - 1:]
        remainder = len(adjusted_seq) % 3
        if remainder != 0:
            adjusted_seq += 'N' * (3 - remainder)
        full_prot = str(Seq(adjusted_seq).translate())

        if '*' in full_prot:
            prefix = full_prot.split('*')[0] + '*'
            fill = ref_full_prot[len(prefix):] if len(ref_full_prot) > len(prefix) else ''
            prot = prefix + fill
        else:
            prot = full_prot

        protein_seqs[sid] = prot
        effective_lengths[sid] = len(prot)

    return protein_seqs, effective_lengths, actual_ref_id


def compute_effective_alignment_lengths(aligned_protein_seqs, original_effective_lengths):
    """
    Once protein sequences have been realigned, figure out how many aligned columns
    belong to each sequence before it reaches its 'true' length.
    """
    effective_aln_lengths = {}
    for sid, seq in aligned_protein_seqs.items():
        count = 0
        orig_len = original_effective_lengths[sid]
        for i, char in enumerate(seq):
            if char != '-':
                count += 1
            if count == orig_len:
                effective_aln_lengths[sid] = i + 1  # 1-based
                break
        if sid not in effective_aln_lengths:
            effective_aln_lengths[sid] = len(seq)
    return effective_aln_lengths


def write_fasta(sequences, output_file):
    with open(output_file, 'w') as f:
        for sid, seq in sequences.items():
            f.write(f'>{sid}\n{seq}\n')

###############################################################################
#                              EXCEL FORMATTING
###############################################################################

def safe_sheet_name(name, max_len=31):
    """Truncate the sheet name to <= 31 characters, required by Excel."""
    if len(name) <= max_len:
        return name
    else:
        return name[:max_len]


def format_excel_columns(writer, df, worksheet, start_header_row=3):
    """
    Auto-fit column widths based on content and bold the column headers.
    'writer' is the ExcelWriter, so we can do 'writer.book.add_format'.
    'start_header_row' is which row the DF column headers occupy.
    """
    workbook = writer.book
    header_format = workbook.add_format({'bold': True})

    # Bold each header cell
    for col_num in range(len(df.columns)):
        worksheet.write(start_header_row, col_num, df.columns[col_num], header_format)

    # Auto-fit columns
    for i, col in enumerate(df.columns):
        if df.empty:
            max_len = len(col) + 2
        else:
            max_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
        worksheet.set_column(i, i, max_len)


def write_matrix_sheet(writer, df, sheet_name):
    """
    Writes the mutation matrix to its own sheet.
    Bold the first row as headers, freeze top row + first column, etc.
    """
    safe_name = safe_sheet_name(sheet_name)
    df.to_excel(writer, sheet_name=safe_name)
    ws = writer.sheets[safe_name]

    # If not empty, bold the header row
    if not df.empty:
        workbook = writer.book
        header_format = workbook.add_format({'bold': True})
        for col_num, value in enumerate(df.columns.values):
            ws.write(0, col_num + 1, value, header_format)

        idx_width = max(df.index.astype(str).map(len).max(), 12)
    else:
        idx_width = 12

    ws.set_column(0, 0, idx_width)
    ws.freeze_panes(1, 1)


def write_analysis_sheet(writer, data, sheet_name, ref_id_used, ref_seq_string, total_seqs):
    """
    Writes the main analysis data to a sheet.
    Also adds lines: "Reference sequence: <ref_id>" and "Total sequences: X"
    and "Reference sequence bases/AA: <ACTUAL_SEQUENCE>" below that.
    """
    safe_name = safe_sheet_name(sheet_name)
    df = pd.DataFrame(data)

    # If it's a protein DataFrame with "Substitution Count", reorder columns
    if 'Substitution Count' in df.columns:
        col_list = list(df.columns)
        if 'Substitution Count' in col_list and 'Insertions' in col_list:
            col_list.remove('Substitution Count')
            insert_index = col_list.index('Insertions')
            col_list.insert(insert_index, 'Substitution Count')
            df = df[col_list]

    df.to_excel(writer, sheet_name=safe_name, index=False, startrow=4)

    ws = writer.sheets[safe_name]
    ws.write(0, 0, f"Reference sequence: {ref_id_used or 'N/A'}")
    ws.write(1, 0, f"Total number of sequences: {total_seqs}")
    ws.write(2, 0, f"Reference sequence bases/AA: {ref_seq_string or ''}")

    if not df.empty:
        format_excel_columns(writer, df, ws, start_header_row=4)
    else:
        ws.write(4, 0, "No mutation data found (empty).")

    ws.freeze_panes(5, 0)

###############################################################################
#                      SUMMARY: NUC AND PROT STATS
###############################################################################

def summarize_for_nuc(df, group_name):
    """Make a summary row for a NUC analysis DataFrame (Nucleotide)."""
    if df.empty:
        return None
    label = "Nucleotide"
    if group_name and group_name != "All":
        label += f" - {group_name}"

    total_pos = len(df)
    mutated = df[df['Total Mutations'] > 0]
    num_mut = len(mutated)
    perc_mut = (num_mut / total_pos * 100) if total_pos else 0
    high_mut = mutated[mutated['Mutation Frequency'] > 0.2]
    num_high = len(high_mut)
    perc_high = (num_high / total_pos * 100) if total_pos else 0

    syn_total = df['Synonymous Mutations'].sum()
    nonsyn_total = df['Non-synonymous Mutations'].sum()
    ins_total = df['Insertions'].sum()
    del_total = df['Deletions'].sum()
    stop_total = df['Stop Codons'].sum()
    mut_total = df['Total Mutations'].sum()

    return {
        'Analysis Type': label,
        'Total Positions': total_pos,
        'Positions with Mutations': num_mut,
        'Percent Positions with Mutations': f"{perc_mut:.2f}",
        'Positions with >20% Mutations': num_high,
        'Percent >20%': f"{perc_high:.2f}",
        'Total Synonymous': syn_total,
        'Total Non-synonymous': nonsyn_total,
        'Total Insertions': ins_total,
        'Total Deletions': del_total,
        'Total Stop Codons': stop_total,
        'Total Mutations': mut_total
    }


def summarize_for_prot(df, group_name):
    """Make a summary row for PROT analysis DataFrame (Protein)."""
    if df.empty:
        return None
    label = "Protein"
    if group_name and group_name != "All":
        label += f" - {group_name}"

    total_pos = len(df)
    mutated = df[df['Total Mutations'] > 0]
    num_mut = len(mutated)
    perc_mut = (num_mut / total_pos * 100) if total_pos else 0
    high_mut = mutated[mutated['Mutation Frequency'] > 0.2]
    num_high = len(high_mut)
    perc_high = (num_high / total_pos * 100) if total_pos else 0

    ins_total = df['Insertions'].sum()
    del_total = df['Deletions'].sum()
    stop_total = df['Stop Codons'].sum()
    sub_count = df['Substitution Count'].sum() if 'Substitution Count' in df.columns else 0
    total_mut = df['Total Mutations'].sum()

    return {
        'Analysis Type': label,
        'Total Positions': total_pos,
        'Positions with Mutations': num_mut,
        'Percent Positions with Mutations': f"{perc_mut:.2f}",
        'Positions with >20% Mutations': num_high,
        'Percent >20%': f"{perc_high:.2f}",
        'Total Substitutions': sub_count,
        'Total Insertions': ins_total,
        'Total Deletions': del_total,
        'Total Stop Codons': stop_total,
        'Total Mutations': total_mut
    }


def write_summary_sheet(writer, summary_nuc, summary_prot):
    """
    We create two separate tables for Nucleotide vs Protein summary,
    each with its own columns and header formatting.
    Fix the KeyError by storing the newly created Worksheet in writer.sheets.
    """
    sheet_name = safe_sheet_name("Summary Statistics")
    ws = writer.book.add_worksheet(sheet_name)
    # Important: store reference in writer.sheets so we can use it
    writer.sheets[sheet_name] = ws

    row_cursor = 0
    bold_format = writer.book.add_format({'bold': True})

    # Nucleotide Summary
    if not summary_nuc.empty:
        ws.write(row_cursor, 0, "Nucleotide Summary Statistics", bold_format)
        row_cursor += 1
        summary_nuc.to_excel(writer, sheet_name=sheet_name, index=False, startrow=row_cursor)
        # Bold the header
        for col_num, value in enumerate(summary_nuc.columns.values):
            ws.write(row_cursor, col_num, value, bold_format)
        row_cursor += len(summary_nuc) + 2

    # Protein Summary
    if not summary_prot.empty:
        ws.write(row_cursor, 0, "Protein Summary Statistics", bold_format)
        row_cursor += 1
        summary_prot.to_excel(writer, sheet_name=sheet_name, index=False, startrow=row_cursor)
        for col_num, value in enumerate(summary_prot.columns.values):
            ws.write(row_cursor, col_num, value, bold_format)
        row_cursor += len(summary_prot) + 2

    # Adjust column widths based on the bigger of the two frames
    if not summary_nuc.empty or not summary_prot.empty:
        # pick whichever has more columns
        if len(summary_nuc.columns) >= len(summary_prot.columns):
            big_df = summary_nuc
        else:
            big_df = summary_prot

        for i, col in enumerate(big_df.columns):
            max_len = max(big_df[col].astype(str).map(len).max(), len(col)) + 2
            ws.set_column(i, i, max_len)


###############################################################################
#                               GROUPING
###############################################################################

def load_group_assignments(csv_file):
    """
    Reads a CSV file with header [Isolate,Category].
    Returns a dict: { 'isolate_id': 'category_value', ... }
    """
    df = pd.read_csv(csv_file)
    group_dict = {}
    for _, row in df.iterrows():
        isolate = str(row[0]).strip()
        category = str(row[1]).strip()
        group_dict[isolate] = category
    return group_dict


###############################################################################
#                                 MAIN
###############################################################################

def main():
    args = parse_arguments()
    setup_logging(debug=args.debug)
    logging.info('Starting Mutation Analysis')

    aln_type = args.type.lower()
    if aln_type not in ['n', 'p', 'both']:
        sys.exit('Invalid alignment type.')

    # Parse alignment
    sequences, aln_length = parse_alignment(args.alignment, aln_type)

    # Possibly detect input type if user wants strict validation
    input_type_guess = guess_sequence_type(sequences)
    logging.info(f"Guessed input type is: {input_type_guess}")

    if args.strict_validation:
        if aln_type in ['n', 'both'] and input_type_guess == 'protein':
            sys.exit("Error: You chose -t n or both, but input looks like protein!")
        if aln_type == 'p' and input_type_guess == 'nucleotide':
            logging.info("Note: Input looks like nucleotides, but -t p chosen. Will auto-translate -> protein...")

    # If user requested a VCF from a presumably nucleotide alignment
    if args.vcf and aln_type in ['n', 'both'] and input_type_guess != 'protein':
        prefix = os.path.splitext(os.path.basename(args.alignment))[0]
        run_snp_sites(args.alignment, snp_sites_path=args.snp_sites_path, output_prefix=prefix)

    # Possibly load groups
    group_dict = {}
    if args.groups:
        logging.info(f"Loading group assignments from {args.groups}")
        group_dict = load_group_assignments(args.groups)

    def run_analysis_for_subset(subset_name, subset_seqs):
        """
        Returns a list of summary dictionaries, plus a list of
        (sheet_name, (df, analysis_type, matrix_df, ref_id, ref_seq, totalSeqs)).
        """
        local_summaries = []
        local_dataframes = []
        totalSeqs = len(subset_seqs)

        # CASE 1: N or BOTH => run Nuc
        if aln_type in ['n', 'both']:
            nuc_data, nuc_matrix, nuc_ref_id = analyze_nucleotide_alignment(
                subset_seqs, aln_length, frame=args.frame,
                reference_id=args.reference
            )
            df_nuc = pd.DataFrame(nuc_data)
            suffix = f"({subset_name})" if subset_name != "All" else ""
            sheet_name = f"NucAnalysis{suffix}"
            matrix_name = f"NucMatrix{suffix}"
            ref_seq_str = subset_seqs.get(nuc_ref_id, "")

            local_dataframes.append((
                sheet_name,
                (df_nuc, 'n', nuc_matrix, nuc_ref_id, ref_seq_str, totalSeqs)
            ))
            summary = summarize_for_nuc(df_nuc, group_name=subset_name)
            if summary:
                local_summaries.append(summary)

        # Decide if we do protein analysis
        do_prot = (aln_type == 'p') or (aln_type == 'both')
        if do_prot:
            # If user chose p but input is nucleotides => we do "translate, realign, analyze"
            if (aln_type == 'p') and (input_type_guess == 'nucleotide'):
                prot_sequences, original_lengths, prot_ref_id = translate_nucleotide_to_protein(
                    subset_seqs, frame=args.frame, reference_id=args.reference
                )
                tmp_prot_in = f"{args.tmp_dir}{args.job_id}_{subset_name}_proteins_in.tmp"
                with open(tmp_prot_in, 'w') as f:
                    for sid, prot_seq in prot_sequences.items():
                        f.write(f'>{sid}\n{prot_seq}\n')

                aligned_prot_fasta = f"{args.tmp_dir}{args.job_id}_{subset_name}_proteins_aligned.tmp"
                run_mafft(tmp_prot_in, mafft_path=args.mafft_path, output_fasta=aligned_prot_fasta)

                prot_seqs_aligned, prot_aln_length = parse_alignment(aligned_prot_fasta, 'protein')
                prot_seqs_aligned = remove_reference_fill(prot_seqs_aligned)
                write_fasta(prot_seqs_aligned, aligned_prot_fasta)

                eff_lengths = compute_effective_alignment_lengths(prot_seqs_aligned, original_lengths)
                prot_data, prot_matrix, prot_actual_ref = analyze_protein_alignment(
                    prot_seqs_aligned, prot_aln_length,
                    effective_aln_lengths=eff_lengths,
                    reference_id=args.reference
                )
                df_prot = pd.DataFrame(prot_data)
                ref_seq_str = prot_seqs_aligned.get(prot_actual_ref, "")
                suffix = f"({subset_name})" if subset_name != "All" else ""
                sheet_name = f"ProtAnalysis{suffix}"
                matrix_name = f"ProtMatrix{suffix}"

                local_dataframes.append((
                    sheet_name,
                    (df_prot, 'p', prot_matrix, prot_actual_ref, ref_seq_str, totalSeqs)
                ))
                summary = summarize_for_prot(df_prot, group_name=subset_name)
                if summary:
                    local_summaries.append(summary)

                if not args.temp:
                    os.remove(tmp_prot_in)
                    os.remove(aligned_prot_fasta)

            elif (aln_type == 'both'):
                # "both": we've already done the N step. Now do protein
                prot_sequences, original_lengths, prot_ref_id = translate_nucleotide_to_protein(
                    subset_seqs, frame=args.frame, reference_id=args.reference
                )
                tmp_prot_in = f"{args.tmp_dir}{args.job_id}_{subset_name}_proteins_in.tmp"
                with open(tmp_prot_in, 'w') as f:
                    for sid, prot_seq in prot_sequences.items():
                        f.write(f'>{sid}\n{prot_seq}\n')

                aligned_prot_fasta = f"{args.tmp_dir}{args.job_id}_{subset_name}_proteins_aligned.tmp"
                run_mafft(tmp_prot_in, mafft_path=args.mafft_path, output_fasta=aligned_prot_fasta)

                prot_seqs_aligned, prot_aln_length = parse_alignment(aligned_prot_fasta, 'protein')
                prot_seqs_aligned = remove_reference_fill(prot_seqs_aligned)
                write_fasta(prot_seqs_aligned, aligned_prot_fasta)

                eff_lengths = compute_effective_alignment_lengths(prot_seqs_aligned, original_lengths)
                prot_data, prot_matrix, prot_actual_ref = analyze_protein_alignment(
                    prot_seqs_aligned, prot_aln_length,
                    effective_aln_lengths=eff_lengths,
                    reference_id=args.reference
                )
                df_prot = pd.DataFrame(prot_data)
                ref_seq_str = prot_seqs_aligned.get(prot_actual_ref, "")
                suffix = f"({subset_name})" if subset_name != "All" else ""
                sheet_name = f"ProtAnalysis{suffix}"
                matrix_name = f"ProtMatrix{suffix}"
                local_dataframes.append((
                    sheet_name,
                    (df_prot, 'p', prot_matrix, prot_actual_ref, ref_seq_str, totalSeqs)
                ))
                summary = summarize_for_prot(df_prot, group_name=subset_name)
                if summary:
                    local_summaries.append(summary)

                if not args.temp:
                    os.remove(tmp_prot_in)
                    os.remove(aligned_prot_fasta)

            else:
                # pure protein analysis if input is actually protein
                prot_data, prot_matrix, prot_ref_id = analyze_protein_alignment(
                    subset_seqs, aln_length,
                    effective_aln_lengths=None,
                    reference_id=args.reference
                )
                df_prot = pd.DataFrame(prot_data)
                ref_seq_str = subset_seqs.get(prot_ref_id, "")
                suffix = f"({subset_name})" if subset_name != "All" else ""
                sheet_name = f"ProtAnalysis{suffix}"
                matrix_name = f"ProtMatrix{suffix}"

                local_dataframes.append((
                    sheet_name,
                    (df_prot, 'p', prot_matrix, prot_ref_id, ref_seq_str, totalSeqs)
                ))
                summary = summarize_for_prot(df_prot, group_name=subset_name)
                if summary:
                    local_summaries.append(summary)

        return local_summaries, local_dataframes

    # Summaries storage
    all_nuc_summaries = []
    all_prot_summaries = []

    # Group or single?
    if group_dict:
        categories = sorted(set(group_dict.values()))
        dataframes = []
        for cat in categories:
            cat_seqs = {k: v for k, v in sequences.items() if group_dict[k] == cat}
            if not cat_seqs:
                logging.warning(f"No sequences found for category {cat}, skipping.")
                continue
            sums, frames = run_analysis_for_subset(cat, cat_seqs)
            for s in sums:
                if s.get('Analysis Type','').startswith("Nucleotide"):
                    all_nuc_summaries.append(s)
                else:
                    all_prot_summaries.append(s)
            dataframes.extend(frames)
    else:
        # Single subset => "All"
        sums, frames = run_analysis_for_subset('All', sequences)
        for s in sums:
            if s.get('Analysis Type','').startswith("Nucleotide"):
                all_nuc_summaries.append(s)
            else:
                all_prot_summaries.append(s)
        dataframes = frames

    # Write final Excel
    output_file = args.output if args.output else f"{args.job_id}_{aln_type}_mutation_analysis.xlsx"
    logging.info(f"Writing output to {output_file}")
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        # Write each analysis + matrix
        for sheet_name, (df, analysis_type, matrix_df, ref_id, ref_seq, totalSeqs) in dataframes:
            # Write analysis
            write_analysis_sheet(
                writer,
                df.to_dict('records'),
                sheet_name,
                ref_id_used=ref_id,
                ref_seq_string=ref_seq,
                total_seqs=totalSeqs
            )
            # Write matrix
            if sheet_name.startswith("NucAnalysis"):
                matrix_sheet = sheet_name.replace("NucAnalysis", "NucMatrix")
            elif sheet_name.startswith("ProtAnalysis"):
                matrix_sheet = sheet_name.replace("ProtAnalysis", "ProtMatrix")
            else:
                matrix_sheet = sheet_name + " Matrix"
            write_matrix_sheet(writer, matrix_df, matrix_sheet)

        # Summaries
        df_nuc = pd.DataFrame(all_nuc_summaries)
        df_prot = pd.DataFrame(all_prot_summaries)
        write_summary_sheet(writer, df_nuc, df_prot)

    logging.info("Mutation Analysis Completed Successfully")


if __name__ == '__main__':
    main()
