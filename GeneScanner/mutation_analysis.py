import argparse
import logging
import sys
import subprocess
import os

from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from typing import List, Tuple  
from Bio.Data.CodonTable import TranslationError
import numpy as np
import pandas as pd
import re 
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
        '-o', '--output', default='./', type=ensure_trailing_slash,
        help=('Output directory (default: current directory). '
              'All result files will be created inside this directory.')
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
        '--min_coverage', default=80.0, type=float,
        help='Minimum per-sequence alignment coverage percentage to keep '
             '(default: 80)'
    )
    parser.add_argument(
        '--min_identity', default=80.0, type=float,
        help='Minimum identity percentage to the reference sequence to keep '
             '(default: 80)'
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
        '--tmp_dir',
        default=None,
        type=ensure_trailing_slash,
        help='Directory used for temporary files (defaults to same as --output)'
    )
    parser.add_argument(
        '--job_id', default='output',
        help='Prefix for output files to avoid clashes in multiple runs'
    )
    # New arguments for grouping and reference selection
    parser.add_argument(
        '--groups',
        help='Optional CSV file: two columns [Isolate,Category] to split analysis by category'
    )
    parser.add_argument(
        '--reference',
        help='Optional reference sequence ID from the FASTA file'
    )
    parser.add_argument(
        '--quiet', action='store_true',
        help='Suppress console output (only errors shown); full log is still '
             'written to output.log'
    )
    
    parser.add_argument(
        '--strict_validation', action='store_true',
        help='If set, error out if the input alignment does not match the chosen -t (n/p/both)'
    )

    parser.add_argument(
        '--reference_csv',
        help=('CSV file with two columns [Isolate,Category]. '
              'Each row chooses the reference isolate for that category.')
    )

    
    return parser.parse_args()

# Compile the forbidden‐character pattern once, at module load:
invalid_xlsx_chars = re.compile(r'[:\\/?*\[\]]')

def ensure_trailing_slash(path):
    if not path.endswith('/'):
        path += '/'
    return path

def safe_sheet_name(name, max_len=31):
    """Truncate the sheet name to <= 31 characters, required by Excel."""
    if len(name) <= max_len:
        return name
    else:
        return name[:max_len]

def setup_logging(debug=False, output_dir='./', quiet=False, job_id='output'):
    """
    Sets up logging to both stdout and to 'output.log' if debug is True.
    """
    level = logging.DEBUG if debug else logging.INFO
        
    handlers = []
    if not quiet:
        handlers.append(logging.StreamHandler(sys.stdout))
    # Always write <job_id>_output.log so removed-sequence names are recorded
    log_filename = f"{job_id}_output.log"
    file_handler = logging.FileHandler(os.path.join(output_dir, log_filename))
    file_handler.setLevel(level)          # INFO or DEBUG, depending on flag
    handlers.append(file_handler)
        
    logging.basicConfig(
        level=level,
        format='%(message)s',
        handlers=handlers
    )

def run_snp_sites(alignment_file, output_dir='', snp_sites_path='snp-sites',output_prefix='out'):
    """
    Runs snp-sites on the provided alignment file to generate a VCF file.
    """
    vcf_file = os.path.join(output_dir, f'{output_prefix}.vcf')
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
    cmd = [mafft_path, '--auto', '--anysymbol', input_fasta]
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

def filter_sequences(
        sequences: dict, aln_length: int,
        *, min_cov: float = 80.0, min_ident: float = 80.0,
        reference_id: str | None = None
) -> Tuple[dict, List[Tuple[str, float, float]]]:

    if not sequences:
        return sequences, []

    ids = list(sequences)
    N, L = len(ids), aln_length
    # Build contiguous bytes then reshape → (N, L) view of dtype='S1'
    arr = np.frombuffer(''.join(sequences[i] for i in ids)
                        .encode('ascii'), dtype='S1').reshape(N, L)

    # Coverage per isolate (% non-gap)
    non_gap = arr != b'-'
    coverage = non_gap.sum(axis=1) / L * 100.0

    # Choose reference row
    ref_idx = ids.index(reference_id) if reference_id in sequences else 0
    ref_row = arr[ref_idx]

    # Identity per isolate (ignore gap positions)
    cmp_mask = non_gap & (ref_row != b'-')
    matches  = (arr == ref_row) & cmp_mask
    identity = matches.sum(axis=1) / cmp_mask.sum(axis=1).clip(min=1) * 100.0

    keep = (coverage >= min_cov) & (identity >= min_ident)
    kept = {ids[i]: sequences[ids[i]] for i in np.where(keep)[0]}
    removed = [(ids[i], float(coverage[i]), float(identity[i]))
               for i in np.where(~keep)[0]]

    for sid, cov, ident in removed:
        logging.info(f'Removed {sid}: coverage={cov:.2f}% identity={ident:.2f}%')

    return kept, removed

def analyze_nucleotide_alignment(sequences, aln_length, frame=1, reference_id=None):
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
    
    # net INDEL balance (+1 for insertion, -1 for deletion) per isolate
    offset       = {sid: 0    for sid in samples}
    window_start = {sid: None for sid in samples}
    fs_windows   = {sid: []   for sid in samples}
    
    ref_seq = sequences[samples[0]]
    total_samples = len(samples)

    # Prepare mutation matrix (rows = sample IDs, columns = 1..aln_length)
    mut_matrix = pd.DataFrame('', index=samples, columns=range(1, aln_length + 1))

    data_rows = []
    # per-isolate mutation counters
    iso_stats = {sid: dict(syn=0, nonsyn=0, ins=0, del_=0, stop=0)
                 for sid in samples}

    # we track how many non-gap reference nt's we've seen so far (nt_count).
    nt_count = 0

    for pos in range(frame - 1, aln_length):
        aln_pos = pos + 1
        ref_nt = ref_seq[pos]

        for sid in samples:
            sample_nt = sequences[sid][pos]
            if ref_nt == '-' and sample_nt != '-':
                offset[sid] = (offset[sid] + 1) % 3          # insertion
            elif ref_nt != '-' and sample_nt == '-':
                offset[sid] = (offset[sid] - 1) % 3          # deletion

        # Determine codon number for this position:
        if ref_nt == '-':
            # If reference has a gap here, codon number is blank.
            current_codon_number = ''
        else:
            # Increment non-gap counter, then compute codon index as 1-based.
            nt_count += 1
            current_codon_number = ((nt_count - 1) // 3) + 1

        if current_codon_number != '':
            for sid in samples:
                # open: we are off-frame and no window yet
                if offset[sid] != 0 and window_start[sid] is None:
                    window_start[sid] = current_codon_number
                # close: frame restored
                elif offset[sid] == 0 and window_start[sid] is not None:
                    fs_windows[sid].append(
                        (window_start[sid], current_codon_number))
                    window_start[sid] = None
                    
            # Calculate the alignment index of the first base of this codon.
            codon_nts = []
            back_pos = pos
            needed = ((nt_count - 1) % 3) + 1  # 1,2 or 3 within the codon
            while back_pos >= 0 and len(codon_nts) < needed:
                if ref_seq[back_pos] != '-':
                    codon_nts.append((back_pos, ref_seq[back_pos]))
                back_pos -= 1

            # Now collect exactly 3 non-gap nt's (if available) to form the codon.
            codon_positions = []
            codon_bases = []
            codon_nts = list(reversed(codon_nts))
            for (idx, base) in codon_nts:
                codon_positions.append(idx)
                codon_bases.append(base)

            # If we only have < 3 bases so far, keep scanning forward from pos+1.
            forward_pos = pos + 1
            while len(codon_bases) < 3 and forward_pos < aln_length:
                if ref_seq[forward_pos] != '-':
                    codon_positions.append(forward_pos)
                    codon_bases.append(ref_seq[forward_pos])
                forward_pos += 1

            # If we still don't have 3, we won't translate (partial codon).
            if len(codon_bases) == 3:
                ref_codon = ''.join(codon_bases)
                try:
                    ref_aa = str(Seq(ref_codon).translate())
                except TranslationError:
                    ref_aa = ''
            else:
                ref_aa = ''
        else:
            # For gaps, leave codon data empty:
            ref_codon = ''
            ref_aa = ''

        syn_count = 0
        nonsyn_count = 0
        ins_count = 0
        del_count = 0
        stop_count = 0

        base2mutation = {}

        for sid in samples:
            sample_nt = sequences[sid][pos]
            if sample_nt == ref_nt:
                continue  # no difference at this position

            # Determine mutation type
            if ref_nt == '-' and sample_nt != '-':
                mtype = 'Insertion'
                ins_count += 1
                display_nt = sample_nt
                iso_stats[sid]['ins'] += 1
                
            elif ref_nt != '-' and sample_nt == '-':
                mtype = 'Deletion'
                del_count += 1
                display_nt = '-'
                iso_stats[sid]['del_'] += 1
                
            else:
                # Both are non-gaps: check codon-level effect if we have ref_codon
                if ref_codon and len(ref_codon) == 3:
                    # Build the sample codon by fetching the same three positions
                    sample_codon_bases = []
                    for cpos in codon_positions:
                        sample_base = sequences[sid][cpos]
                        if sample_base != '-':
                            sample_codon_bases.append(sample_base)
                        else:
                            # If the sample has a gap inside the codon, skip translation
                            sample_codon_bases = []
                            break

                    if len(sample_codon_bases) == 3:
                        sample_codon = ''.join(sample_codon_bases)
                        try:
                            sample_aa = str(Seq(sample_codon).translate())
                        except TranslationError:
                            sample_aa = ''
                    else:
                        sample_aa = ''
                else:
                    sample_aa = ''

                if sample_aa == '*':
                    mtype = 'Non-synonymous Stop Codon'
                    stop_count += 1
                    display_nt = sample_nt
                    iso_stats[sid]['stop'] += 1
                    # premature stop closes any open window
                    if window_start[sid] is not None:
                        fs_windows[sid].append(
                            (window_start[sid], current_codon_number))
                        window_start[sid] = None
                    
                elif sample_aa and ref_aa and sample_aa == ref_aa:
                    mtype = 'Synonymous'
                    syn_count += 1
                    display_nt = sample_nt
                    iso_stats[sid]['syn'] += 1
                    
                elif sample_aa and ref_aa and sample_aa != ref_aa:
                    mtype = 'Non-synonymous'
                    nonsyn_count += 1
                    display_nt = sample_nt
                    iso_stats[sid]['nonsyn'] += 1
                    
                else:
                    mtype = 'Unknown'
                    display_nt = sample_nt

            base2mutation[display_nt] = mtype
            mut_matrix.at[sid, aln_pos] = display_nt

        # Summarize for this position
        indel_stop = ins_count + del_count + stop_count
        total_mut = syn_count + nonsyn_count + indel_stop
        mut_freq = (total_mut / total_samples) if total_samples else 0

        row_data = {
            'Alignment Position': aln_pos,
            'Codon Number': current_codon_number,
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

    # ── close any FS window still open at alignment end ─────────────
    for sid in samples:
        if window_start[sid] is not None:
            fs_windows[sid].append((window_start[sid], nt_count // 3))

    data_df = pd.DataFrame(data_rows)
    return (data_df.to_dict('records'),
            mut_matrix,
            actual_ref_id,
            fs_windows,
            iso_stats)

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


def analyze_protein_alignment(sequences, aln_length, *, effective_aln_lengths=None, reference_id=None, frameshift_windows=None):
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
    column2codon = []
    codon_ctr = 0
    for pos in range(aln_length):
        ref_aa_here = sequences[ref_seq_id][pos]
        if ref_aa_here != '-':
            codon_ctr += 1
            column2codon.append(codon_ctr)
        else:
            column2codon.append('')
            
    iso_stats = {sid: dict(sub=0, sub_fs=0, ins=0, del_=0,
                           stop=0, stop_fs=0)
                 for sid in samples}
    
    for sid, seq in sequences.items():
        idx = seq.find('*')
        stop_positions[sid] = idx if idx != -1 else None

    mut_matrix = pd.DataFrame('', index=samples, columns=range(1, aln_length + 1))
    
    for pos in range(aln_length):
        aln_pos = pos + 1
        wt_aa = sequences[ref_seq_id][pos]
        current_codon_number = column2codon[pos]   # ← needed for FS test

        ins_count = 0
        del_count = 0
        stop_count = 0
        stop_fs_count = 0
        sub_reg_count = 0
        sub_fs_count = 0
        mutations = {}
        survived_count = 0

        for sid in samples:
            if effective_aln_lengths is not None and pos >= effective_aln_lengths[sid]:
                continue
            if stop_positions[sid] is not None and pos > stop_positions[sid]:
                continue

            survived_count += 1
            sample_aa = sequences[sid][pos]

            if sample_aa == wt_aa:
                continue

            # decide mutation type and what to display
            if wt_aa == '-' and sample_aa != '-':
                mtype = 'Insertion'
                ins_count += 1
                display_aa = sample_aa
                iso_stats[sid]['ins'] += 1
                
            elif wt_aa != '-' and sample_aa == '-':
                mtype = 'Deletion'
                del_count += 1
                display_aa = '-'
                iso_stats[sid]['del_'] += 1
            
            else:
                if sample_aa == '*':
                    fs_active = (
                        frameshift_windows and
                        any(start <= current_codon_number <= end
                            for start, end in frameshift_windows.get(sid, []))
                    )
                    if fs_active:
                        mtype = 'Stop Codon(Frameshift)'
                        stop_fs_count += 1
                        iso_stats[sid]['stop_fs'] += 1
                    else:
                        mtype = 'Stop Codon'
                        stop_count += 1
                        iso_stats[sid]['stop'] += 1
                    display_aa = '*'
                else:
                    fs_active = (
                        frameshift_windows and
                        any(start <= current_codon_number <= end
                            for start, end in frameshift_windows.get(sid, []))
                    )
                    if fs_active:
                        mtype = 'Substitution(Frameshift)'
                        sub_fs_count += 1
                        iso_stats[sid]['sub_fs'] += 1
                    else:
                        mtype = 'Substitution'
                        sub_reg_count += 1
                        iso_stats[sid]['sub'] += 1
                    display_aa = sample_aa
                
            mutations[display_aa] = mtype
            mut_matrix.at[sid, aln_pos] = display_aa

        total_mut = ins_count + del_count + stop_count + sub_reg_count + sub_fs_count
        mut_freq = (total_mut / survived_count) if survived_count > 0 else 0

        row_dict = {
            'Alignment Position': aln_pos,
            'Reference AA':        wt_aa,
            'Mutations':           ','.join(sorted(mutations.keys())),
            'Mutation Types':      ','.join(mutations[aa] for aa in sorted(mutations.keys())),
            'Substitution Count':  sub_reg_count,
            'Substitution(Frameshift) Count': sub_fs_count,
            'Insertions':          ins_count,
            'Deletions':           del_count,
            'Stop Codons':         stop_count,
            'Stop Codons(Frameshift)': stop_fs_count,
            'Total Mutations':     total_mut,
            'Mutation Frequency':  mut_freq,
            'Remaining sequences before encountering stop codon': survived_count
        }
        data.append(row_dict)

    return data, mut_matrix, actual_ref_id, iso_stats
    

def translate_nucleotide_to_protein(sequences, frame=1, reference_id=None):
    """
    Translates each ungapped nucleotide sequence into protein.
    Returns a triple:
      - dict_of_prot_sequences,
      - dict_of_effective_lengths,
      - actual_ref_id_used
    """
    logging.info('Translating nucleotide sequences (ungapped) to proteins...')
    protein_seqs = {}
    effective_lengths = {}

    # Decide which sequence ID to treat as reference
    all_ids = list(sequences.keys())
    if reference_id and reference_id in all_ids:
        actual_ref_id = reference_id
    else:
        actual_ref_id = all_ids[0] if all_ids else None

    # Build the reference translation for filling
    ref_nt = sequences.get(actual_ref_id, '').replace('-', '')
    ref_adj = ref_nt[frame-1:]
    rem = len(ref_adj) % 3
    if rem:
        ref_adj += 'N' * (3 - rem)
    ref_full_prot = str(Seq(ref_adj).translate())

    for sid, nt_seq in sequences.items():
        raw = nt_seq.replace('-', '')
        seq_adj = raw[frame-1:]
        rem = len(seq_adj) % 3
        if rem:
            seq_adj += 'N' * (3 - rem)
        prot_full = str(Seq(seq_adj).translate())

        if '*' in prot_full:
            prefix = prot_full.split('*')[0] + '*'
            fill = ref_full_prot[len(prefix):] if len(ref_full_prot) > len(prefix) else ''
            prot = prefix + fill
        else:
            prot = prot_full

        protein_seqs[sid] = prot
        effective_lengths[sid] = len(prot)

    return protein_seqs, effective_lengths, actual_ref_id


def compute_effective_alignment_lengths(aligned_protein_seqs, original_effective_lengths):
    effective_aln_lengths = {}
    for sid, seq in aligned_protein_seqs.items():
        count = 0
        orig_len = original_effective_lengths[sid]
        for i, char in enumerate(seq):
            if char != '-':
                count += 1
            if count == orig_len:
                effective_aln_lengths[sid] = i + 1  # using 1-based indexing in alignment
                break
        if sid not in effective_aln_lengths:
            effective_aln_lengths[sid] = len(seq)
    return effective_aln_lengths


def write_fasta(sequences, output_file):
    with open(output_file, 'w') as f:
        for sid, seq in sequences.items():
            f.write(f'>{sid}\n{seq}\n')


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
    df.to_excel(writer, sheet_name=sheet_name)
    ws = writer.sheets[sheet_name]
    # Format headers
    for col_num, value in enumerate(df.columns.values):
        ws.write(0, col_num + 1, value)  # +1 offset for index column
    ws.set_column(0, 0, max(df.index.astype(str).map(len).max(), 12))
    for i, col in enumerate(df.columns):
        max_len = max(df[col].astype(str).map(len).max(), len(str(col))) + 2
        ws.set_column(i + 1, i + 1, max_len)


def write_analysis_sheet(writer, data, sheet_name, ref_id_used, ref_seq_string, total_seqs):
    """
    Writes the main analysis data to a sheet.
    Also adds lines: "Reference sequence: <ref_id>" and "Total sequences: X"
    and "Reference sequence bases/AA: <ACTUAL_SEQUENCE>" below that.
    """
    df = pd.DataFrame(data)

    # If it's a protein DataFrame with "Substitution Count", reorder columns
    if 'Substitution Count' in df.columns:
        col_list = list(df.columns)
        if 'Substitution Count' in col_list and 'Insertions' in col_list:
            col_list.remove('Substitution Count')
            insert_index = col_list.index('Insertions')
            col_list.insert(insert_index, 'Substitution Count')
            df = df[col_list]

    df.to_excel(writer, sheet_name=sheet_name, index=False, startrow=4)

    ws = writer.sheets[sheet_name]
    ws.write(0, 0, f"Reference sequence: {ref_id_used or 'N/A'}")
    ws.write(1, 0, f"Total number of sequences: {total_seqs}")
    ws.write(2, 0, f"Reference sequence bases/AA: {ref_seq_string or ''}")

    if not df.empty:
        format_excel_columns(writer, df, ws, start_header_row=4)
    else:
        ws.write(4, 0, "No mutation data found (empty).")

    ws.freeze_panes(5, 0)


def summarize_for_nuc(df, group_name, iso_stats):
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
    stop_total     = df['Stop Codons'].sum()
    stop_fs_total  = df['Stop Codon(Frameshift)'].sum() \
                     if 'Stop Codon(Frameshift)' in df.columns else 0
    mut_total = df['Total Mutations'].sum()
    
    # isolate-level additions
    n_iso = len(iso_stats)
    no_syn   = sum(1 for st in iso_stats.values() if st['syn'   ] == 0)
    no_nsyn  = sum(1 for st in iso_stats.values() if st['nonsyn'] == 0)
    no_ins   = sum(1 for st in iso_stats.values() if st['ins'   ] == 0)
    no_del   = sum(1 for st in iso_stats.values() if st['del_'  ] == 0)
    no_stop  = sum(1 for st in iso_stats.values() if st['stop'  ] == 0)


    return {
        'Analysis Type': label,
        'Total Positions': total_pos,
        'Total Isolates': n_iso,
        'Positions with Mutations': num_mut,
        'Percent Positions with Mutations': f"{perc_mut:.2f}",
        'Positions with >20% Mutations': num_high,
        'Percent >20%': f"{perc_high:.2f}",
        'Total Synonymous': syn_total,
        'Total Non-synonymous': nonsyn_total,
        'Total Insertions': ins_total,
        'Total Deletions': del_total,
        'Total Stop Codons': stop_total,
        'Total Stop Codons(Frameshift)':   stop_fs_total,
        'Total Mutations': mut_total,
        'Isolate w/o synonymous mutation':     no_syn,
        'Isolate w/o non-synonymous mutation': no_nsyn,
        'Isolate w/o Insertion':               no_ins,
        'Isolate w/o Deletion':                no_del,
        'Isolate w/o Stop codons':             no_stop
    }


def summarize_for_prot(df, group_name, iso_stats):
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
    reg_sub_total    = df['Substitution Count'].sum() if 'Substitution Count' in df.columns else 0
    total_frameshift = df['Substitution(Frameshift) Count'].sum() \
                       if 'Substitution(Frameshift) Count' in df.columns else 0
    total_mut        = df['Total Mutations'].sum()

    n_iso   = len(iso_stats)
    no_sub  = sum(1 for st in iso_stats.values() if st['sub'    ] == 0)
    no_fs   = sum(1 for st in iso_stats.values() if st['sub_fs' ] == 0)
    no_ins  = sum(1 for st in iso_stats.values() if st['ins'    ] == 0)
    no_del  = sum(1 for st in iso_stats.values() if st['del_'   ] == 0)
    no_stop    = sum(1 for st in iso_stats.values() if st['stop'    ] == 0)
    no_stop_fs = sum(1 for st in iso_stats.values() if st['stop_fs'] == 0)

    return {
        'Analysis Type': label,
        'Total Isolates': n_iso,
        'Total Positions': total_pos,
        'Positions with Mutations': num_mut,
        'Percent Positions with Mutations': f"{perc_mut:.2f}",
        'Positions with >20% Mutations': num_high,
        'Percent >20%': f"{perc_high:.2f}",
        'Total Substitutions': sub_count,
        'Total Insertions': ins_total,
        'Total Deletions': del_total,
        'Total Stop Codons': stop_total,
        'Total Substitutions(Regular)': reg_sub_total,
        'Total Substitution(Frameshift)': total_frameshift,
        'Total Mutations':                total_mut,
        'Isolate w/o Substitution':      no_sub,
        'Isolate w/o Frameshift':        no_fs,
        'Isolate w/o Insertion':         no_ins,
        'Isolate w/o Deletion':          no_del,
        'Isolate w/o Stop codons':       no_stop,
        'Isolate w/o Stop codon(Frameshift)': no_stop_fs
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
            
def load_reference_assignments(csv_file):
    """
    CSV with header [Isolate, Category]; returns
        { 'category_value': 'isolate_id',  ... }
    If multiple isolates are marked for the same category, the first wins.
    """
    df = pd.read_csv(csv_file)
    ref_dict = {}
    for _, row in df.iterrows():
        isolate  = str(row.iloc[0]).strip()
        category = str(row.iloc[1]).strip()
        # keep only the first isolate per category
        ref_dict.setdefault(category, isolate)
    return ref_dict

def load_group_assignments(csv_file):
    """
    Reads a CSV file with header [Isolate, Category].
    Returns a dict: { 'isolate_id': 'category_value', ... }
    """
    df = pd.read_csv(csv_file)
    group_dict = {}
    for _, row in df.iterrows():
        isolate = str(row.iloc[0]).strip()
        category = str(row.iloc[1]).strip()
        group_dict[isolate] = category
    return group_dict

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

def main():
    args = parse_arguments()
    # ensure output directory exists and configure logging there
    os.makedirs(args.output, exist_ok=True)
    setup_logging(
        debug=args.debug,
        output_dir=args.output,
        quiet=args.quiet,
        job_id=args.job_id
    )

    # if --tmp_dir was not given, default it to the output directory
    if not args.tmp_dir:
        args.tmp_dir = args.output
    else:
        # make sure any custom path ends in a slash
        args.tmp_dir = ensure_trailing_slash(args.tmp_dir)
    # ensure temp directory exists
    os.makedirs(args.tmp_dir, exist_ok=True)
    logging.info('Starting Mutation Analysis')

    aln_type = args.type.lower()
    if aln_type not in ['n', 'p', 'both']:
        sys.exit('Invalid alignment type.')

    # Parse alignment
    all_sequences, aln_length = parse_alignment(args.alignment, aln_type)

    all_sequences, removed_seqs = filter_sequences(
        all_sequences, aln_length,
        min_cov=args.min_coverage, min_ident=args.min_identity,
        reference_id=args.reference
    )
    if not all_sequences or len(all_sequences) < 2:
        sys.exit('No sequences passed the coverage/identity filters.')

    if args.reference and args.reference in all_sequences:
        global_ref_id = args.reference
    else:
        global_ref_id = next(iter(all_sequences))   # first sequence that survived

    # Possibly detect input type if user wants strict validation
    input_type_guess = guess_sequence_type(all_sequences)
    logging.info(f"Guessed input type is: {input_type_guess}")

    ref_csv_dict = {}
    if args.reference_csv:
        logging.info(f"Loading per-category references from {args.reference_csv}")
        ref_csv_dict = load_reference_assignments(args.reference_csv)
    
    if args.strict_validation:
        if aln_type in ['n', 'both'] and input_type_guess == 'protein':
            sys.exit("Error: You chose -t n or both, but input looks like protein!")
        if aln_type == 'p' and input_type_guess == 'nucleotide':
            logging.info("Note: Input looks like nucleotides, but -t p chosen. Will auto-translate -> protein...")

    # If requested, generate VCF for NT alignments
    if args.vcf and (aln_type in ['n', 'both']):
        run_snp_sites(
            args.alignment,
            output_dir=args.output,
            snp_sites_path=args.snp_sites_path,
            output_prefix=args.job_id
        )

    # Possibly load group assignments
    group_dict = {}
    if args.groups:
        logging.info(f"Loading group assignments from {args.groups}")
        group_dict = load_group_assignments(args.groups)
        # prune assignments that correspond to filtered-out isolates
        group_dict = {sid: cat for sid, cat in group_dict.items()
                      if sid in all_sequences}
    
    if ref_csv_dict:
        if args.groups:
            cats = sorted(set(group_dict.values()), key=str)
        else:
            cats = sorted(ref_csv_dict.keys(), key=str)

        parts = [
            f"{ref_csv_dict.get(cat, global_ref_id)} for group {cat}"
            for cat in cats
        ]
        logging.info("Different references are given: " + " and ".join(parts))
    else:
        logging.info(f"Global reference isolate for all analyses: {global_ref_id}")
        
    # --------------------------------------------------------------
    # Build a single protein alignment 
    # --------------------------------------------------------------
    need_translate_global = (
        (args.type == "p" and input_type_guess == "nucleotide")
        or (args.type == "both")
    )

    if args.type in ("p", "both"):
        if need_translate_global:
            # 1. translate every nucleotide sequence to protein
            prot_raw_all, eff_len_all, _ = translate_nucleotide_to_protein(
                all_sequences, frame=args.frame, reference_id=global_ref_id
            )

            # 2. run MAFFT **once** on the full protein set
            tmp_in  = f"{args.tmp_dir}{args.job_id}_prot_in.tmp"
            out_fa  = f"{args.output}{args.job_id}_prot_aligned.fasta"
            write_fasta(prot_raw_all, tmp_in)
            run_mafft(tmp_in, mafft_path=args.mafft_path, output_fasta=out_fa)

            # 3. parse & post-process
            prot_aln_all, prot_aln_len_all = parse_alignment(out_fa, "protein")
            prot_aln_all = remove_reference_fill(prot_aln_all)
            write_fasta(prot_aln_all, out_fa)          # keep a cleaned copy

            # 4. effective alignment length per isolate
            eff_aln_len_all = compute_effective_alignment_lengths(
                prot_aln_all, eff_len_all
            )

            if not args.temp:
                os.remove(tmp_in)
        else:
            # input is already protein; no realignment needed
            prot_aln_all       = all_sequences
            prot_aln_len_all   = aln_length
            eff_aln_len_all    = {sid: aln_length for sid in all_sequences}

    # ------------------------------------------------------------------
    # Build a list of (subset_name, seq_dict, reference_id)
    # ------------------------------------------------------------------
    subset_specs = []
    if group_dict:
        for cat in sorted(set(group_dict.values())):
            cat_seqs = {k: v for k, v in all_sequences.items()
                        if group_dict.get(k) == cat}

            ref_id = ref_csv_dict.get(cat, global_ref_id)

            # make sure reference(s) are present
            if ref_id in all_sequences:
                cat_seqs.setdefault(ref_id, all_sequences[ref_id])
            cat_seqs.setdefault(global_ref_id, all_sequences[global_ref_id])

            if cat_seqs:
                subset_specs.append((cat, cat_seqs, ref_id))
    else:
        subset_specs.append(("All", all_sequences, global_ref_id))

    # ------------------------------------------------------------------
    # Run analyses for each subset
    # ------------------------------------------------------------------
    all_nuc_summaries, all_prot_summaries = [], []
    dataframes = []   # (sheet, (df, flag, matrix, ref, ref_seq, nSeq))

    for subset_name, subset_seqs, subset_ref in subset_specs:
        totalSeqs   = len(subset_seqs)
        tag         = "" if subset_name == "All" else f"({subset_name})"

        fs_windows = set()   # default

        # ── nucleotide analysis ───────────────────────────────────────
        if args.type in ("n", "both"):
            nuc_rows, nuc_matrix, nuc_ref, fs_windows, iso_stats_nuc = \
                analyze_nucleotide_alignment(
                subset_seqs, aln_length, frame=args.frame,
                reference_id=subset_ref
            )
            df_nuc = pd.DataFrame(nuc_rows)
            sheet  = f"NucAnalysis{tag}"

            dataframes.append((sheet, (df_nuc, 'n', nuc_matrix,
                                nuc_ref, subset_seqs.get(nuc_ref, ""), totalSeqs)))

            summary = summarize_for_nuc(df_nuc, subset_name, iso_stats_nuc)
            if summary:
                all_nuc_summaries.append(summary)

        # ── protein analysis ──────────────────────────────────────────        
        if args.type in ("p", "both"):
            # carve out just the isolates that belong to this subset
            prot_aln_subset = {k: prot_aln_all[k] for k in subset_seqs}
            eff_len_subset  = {k: eff_aln_len_all[k] for k in subset_seqs}

            prot_rows, prot_matrix, prot_ref, iso_stats_prot = \
                analyze_protein_alignment(
                    prot_aln_subset, prot_aln_len_all,
                    effective_aln_lengths=eff_len_subset,
                    reference_id=subset_ref,
                    frameshift_windows=fs_windows
                )

            df_prot = pd.DataFrame(prot_rows)
            sheet   = f"ProtAnalysis{tag}"

            dataframes.append((sheet, (df_prot, 'p', prot_matrix,
                                prot_ref, subset_seqs.get(prot_ref, ""), totalSeqs)))

            summary = summarize_for_prot(df_prot, subset_name, iso_stats_prot)
            if summary:
                all_prot_summaries.append(summary)

    # Write final Excel output
    output_file = os.path.join(
        args.output,
        f'{args.job_id}_{aln_type}_mutation_analysis.xlsx'
    )
    logging.info(f'Writing output to {output_file}')
    
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        used_sheet_names = set()
        # Write each analysis data frame & matrix
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
            used_sheet_names.add(sheet_name)
            
        #summaries  
        df_nuc = pd.DataFrame(all_nuc_summaries)
        df_prot = pd.DataFrame(all_prot_summaries)
        write_summary_sheet(writer, df_nuc, df_prot)

    logging.info("Mutation Analysis Completed Successfully")


if __name__ == '__main__':
    main()
