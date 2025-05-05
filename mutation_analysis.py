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
    
    return parser.parse_args()

# Compile the forbidden‐character pattern once, at module load:
invalid_xlsx_chars = re.compile(r'[:\\/?*\[\]]')

def ensure_trailing_slash(path):
    if not path.endswith('/'):
        path += '/'
    return path

def safe_sheet_name(name: str, existing: set) -> str:
    name = invalid_xlsx_chars.sub('', name).strip()
    if len(name) > 31:
        name = name[:30] + '…'          # keeps total length 31
    base = name
    i = 1
    while name in existing:
        suffix = f'_{i}'
        if len(base) + len(suffix) > 31:
            base = base[:31 - len(suffix)]
        name = base + suffix
        i += 1
    existing.add(name)
    return name

def setup_logging(debug=False, output_dir='./', quiet=False):
    """
    Sets up logging to both stdout and to 'output.log' if debug is True.
    """
    level = logging.DEBUG if debug else logging.INFO
        
    handlers = []
    if not quiet:
        handlers.append(logging.StreamHandler(sys.stdout))
    # Always write output.log so removed-sequence names are recorded
    file_handler = logging.FileHandler(os.path.join(output_dir, 'output.log'))
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
    logging.info('Analyzing nucleotide alignment (syn/non-syn)...')
    data = []

    # The sample order matters (we want reference first, if specified)
    all_ids = list(sequences.keys())
    if reference_id and reference_id in all_ids:
        all_ids.remove(reference_id)
        samples = [reference_id] + all_ids
    else:
        samples = all_ids

    ref_seq = sequences[samples[0]]
    total_samples = len(samples)

    start_pos = frame - 1
    usable_length = aln_length - start_pos
    codon_count = usable_length // 3

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
                    mtype = 'Unknown'
                    if sample_aa and sample_aa != ref_aa:
                        mutated_aa_set.add(sample_aa)

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


def analyze_protein_alignment(sequences, aln_length, effective_aln_lengths=None, reference_id=None):
    logging.info('Analyzing protein alignment...')

    # The sample order matters (we want reference first, if specified)
    all_ids = list(sequences.keys())
    if reference_id and reference_id in all_ids:
        all_ids.remove(reference_id)
        samples = [reference_id] + all_ids
    else:
        samples = all_ids

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
        mut_rate = (total_mut / survived_count) if survived_count > 0 else 0
        
        data.append({
            'Alignment Position': aln_pos,
            'Wildtype AA': wt_aa,
            'Mutations': ','.join(sorted(mutations.keys())),
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
    

def translate_nucleotide_to_protein(sequences, frame=1, reference_id=None):
    logging.info('Translating nucleotide sequences (ungapped) to proteins...')
    protein_seqs = {}
    effective_lengths = {}

    # If a reference was requested, we want to use that as the "ref" for the fill
    all_ids = list(sequences.keys())
    if reference_id and reference_id in all_ids:
        ref_id = reference_id
    else:
        ref_id = next(iter(sequences))  # fallback if the requested reference isn't found

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

    return protein_seqs, effective_lengths


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


def format_excel_columns(df, worksheet):
    # Auto-size columns based on content
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


def write_analysis_sheet(writer, data, sheet_name, ref_id):
    df = pd.DataFrame(data)
    if df.empty:
        return None
    df.to_excel(writer, sheet_name=sheet_name, index=False, startrow=2)
    ws = writer.sheets[sheet_name]
    ws.write(0, 0, f'Reference sequence: {ref_id}')
    for col_num, value in enumerate(df.columns.values):
        ws.write(2, col_num, value)
    format_excel_columns(df, ws)
    return df


def summarize_data(df, analysis_type):
    total_pos = len(df)
    mutated = df[df['Total Mutations'] > 0]
    num_mut = len(mutated)
    perc_mut = (num_mut / total_pos * 100) if total_pos else 0

    # Here, we just demonstrate one possible "high mutation" threshold
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
    else:  # p
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


def load_group_assignments(csv_file):
    """
    Reads a CSV file with header [Isolate, Category].
    Returns a dict: { 'isolate_id': 'category_value', ... }
    """
    df = pd.read_csv(csv_file)
    group_dict = {}
    for _, row in df.iterrows():
        isolate = str(row[0]).strip()
        category = str(row[1]).strip()
        group_dict[isolate] = category
    return group_dict


def main():
    args = parse_arguments()
    # ensure output directory exists and configure logging there
    os.makedirs(args.output, exist_ok=True)
    setup_logging(debug=args.debug, output_dir=args.output, quiet=args.quiet)

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
    if aln_type == 'n':
        aln_type_full = 'nucleotide'
    elif aln_type == 'p':
        aln_type_full = 'protein'
    elif aln_type == 'both':
        aln_type_full = 'both'
    else:
        sys.exit('Invalid alignment type.')

    # Parse alignment
    all_sequences, aln_length = parse_alignment(args.alignment, aln_type_full)

    all_sequences, removed_seqs = filter_sequences(
        all_sequences, aln_length,
        min_cov=args.min_coverage, min_ident=args.min_identity,
        reference_id=args.reference
    )
    if not all_sequences or len(all_sequences) < 2:
        sys.exit('No sequences passed the coverage/identity filters.')

    # If requested, generate VCF for NT alignments
    if args.vcf and (aln_type in ['n', 'both']):
        prefix = os.path.splitext(os.path.basename(args.alignment))[0]
        run_snp_sites(args.alignment, output_dir=args.output,snp_sites_path=args.snp_sites_path, output_prefix=prefix)

    # Possibly load group assignments
    group_dict = {}
    if args.groups:
        logging.info(f"Loading group assignments from {args.groups}")
        group_dict = load_group_assignments(args.groups)
        # prune assignments that correspond to filtered-out isolates
        group_dict = {sid: cat for sid, cat in group_dict.items()
                      if sid in all_sequences}

    # A helper function to gather analyses from a subset of sequences
    def run_analysis_for_subset(subset_name, subset_sequences):
        """
        Returns:
          - summaries list
          - all dataframes for final writing
        """
        local_summaries = []
        local_dataframes = {}  # For writing separate sheets, e.g. {"Nucleotide Analysis (GroupX)": df, ...}

        nuc_data, nuc_matrix = None, None
        prot_data, prot_matrix = None, None
        prot_sequences = {}

        # Determine which sequence ID to use as reference for this subset
        if args.reference and args.reference in subset_sequences:
            subset_ref_id = args.reference
        else:
            subset_ref_id = next(iter(subset_sequences))
        # Pull out the actual ungapped reference sequence string
        subset_ref_seq = subset_sequences[subset_ref_id].replace('-', '')


        if aln_type == 'n':
            nuc_data, nuc_matrix = analyze_nucleotide_alignment(
                subset_sequences, aln_length, frame=args.frame,
                reference_id=args.reference
            )

        elif aln_type == 'p':
            prot_data, prot_matrix = analyze_protein_alignment(
                subset_sequences, aln_length,
                reference_id=args.reference
            )
            prot_sequences = subset_sequences

        elif aln_type == 'both':
            # 1) Nucleotide
            nuc_data, nuc_matrix = analyze_nucleotide_alignment(
                subset_sequences, aln_length, frame=args.frame,
                reference_id=args.reference
            )
            # 2) Translate -> realign in protein space -> analyze
            prot_sequences, original_effective_lengths = translate_nucleotide_to_protein(
                subset_sequences, frame=args.frame, reference_id=args.reference
            )
            tmp_prot_in = f'{args.tmp_dir}{args.job_id}_{subset_name}_proteins_in.tmp'
            with open(tmp_prot_in, 'w') as f:
                for sid, prot_seq in prot_sequences.items():
                    f.write(f'>{sid}\n{prot_seq}\n')

            aligned_prot_fasta = f'{args.tmp_dir}{args.job_id}_{subset_name}_proteins_aligned.tmp'
            run_mafft(tmp_prot_in, mafft_path=args.mafft_path, output_fasta=aligned_prot_fasta)

            prot_seqs_aligned, prot_aln_length = parse_alignment(aligned_prot_fasta, 'protein')
            prot_seqs_aligned = remove_reference_fill(prot_seqs_aligned)

            write_fasta(prot_seqs_aligned, aligned_prot_fasta)
            effective_aln_lengths = compute_effective_alignment_lengths(prot_seqs_aligned, original_effective_lengths)

            prot_data, prot_matrix = analyze_protein_alignment(
                prot_seqs_aligned, prot_aln_length,
                effective_aln_lengths=effective_aln_lengths,
                reference_id=args.reference
            )
            prot_sequences = prot_seqs_aligned

            if not args.temp:
                os.remove(tmp_prot_in)
                os.remove(aligned_prot_fasta)

        # Convert results into DataFrames and store them
        if nuc_data:
            df_nuc = pd.DataFrame(nuc_data)
            if not df_nuc.empty:
                sheet_name = f"Nucleotide Analysis ({subset_name})"
                # store the ungapped sequence, not the ID
                local_dataframes[sheet_name] = (df_nuc, 'n', nuc_matrix, subset_ref_seq)
        if prot_data:
            df_prot = pd.DataFrame(prot_data)
            if not df_prot.empty:
                sheet_name = f"Protein Analysis ({subset_name})"
                # store the ungapped sequence, not the ID
                local_dataframes[sheet_name] = (df_prot, 'p', prot_matrix, subset_ref_seq)

        # Summaries
        #   unpack the extra ref_id but we only need df + analysis_type
        for name, (df, analysis_type, matrix, _ref_id) in local_dataframes.items():
            local_summaries.append(summarize_data(df, analysis_type))

        return local_summaries, local_dataframes

    # Decide on subsets
    if group_dict:
        # Build subsets by category
        categories = sorted(set(group_dict.values()))
        # We will produce a combined analysis across categories and add sheets for each category
        # or only by category? The code below runs for each category only. 
        # If you also want an "Overall" sheet, you can add that logic too.

        # We'll gather all results in these accumulators
        overall_summaries = []
        overall_dataframes = {}

        for cat in categories:
            # Filter sequences
            cat_sequences = {
                sid: seq for sid, seq in all_sequences.items() if group_dict.get(sid) == cat
            }
            if not cat_sequences:
                logging.warning(f"No sequences found for category {cat}, skipping.")
                continue

            cat_summaries, cat_dataframes = run_analysis_for_subset(cat, cat_sequences)
            overall_summaries.extend(cat_summaries)
            overall_dataframes.update(cat_dataframes)

    else:
        # Just one big subset = entire dataset
        overall_summaries, overall_dataframes = run_analysis_for_subset('All', all_sequences)

    # Write final Excel output
    output_file = os.path.join(
        args.output,
        f'{args.job_id}_{aln_type_full}_mutation_analysis.xlsx'
    )
    logging.info(f'Writing output to {output_file}')
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        used_sheet_names = set()
        # Write each analysis data frame & matrix
        for raw_name, (df, analysis_type, matrix_df, ref_id) in overall_dataframes.items():
            sheet_name = safe_sheet_name(raw_name, used_sheet_names)
            # write with the explicit ref_id captured earlier
            df_out = write_analysis_sheet(
                writer,
                df.to_dict('records'),
                sheet_name,
                ref_id=ref_id
            )
            # The above line sets "reference sequence" text, though not super-critical
            # If matrix_df is not None, write that on another sheet
            if matrix_df is not None:
                matrix_sheet = safe_sheet_name(sheet_name + " Matrix", used_sheet_names)
                write_matrix_sheet(writer, matrix_df, matrix_sheet)

        # Summaries
        write_summary_sheet(writer, overall_summaries)

    logging.info("Mutation Analysis Completed Successfully")


if __name__ == '__main__':
    main()
