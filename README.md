# GeneScanner v1.0

Mutation Analysis Pipeline
==========================

This pipeline analyzes either **nucleotide** or **protein** aligned
FASTA files (or both), identifies and classifies mutations, produces
summary statistics, and writes a **multi-sheet Excel** report. It
supports:

1.  **Codon-based Nucleotide Analysis** (detects synonymous vs.
    non-synonymous changes).

2.  **Protein Analysis** (detects substitutions, insertions, deletions,
    and stop codons).

3.  **Dual (“both”) Analysis** (translates nucleotides → realigns in
    protein space → runs protein analysis).

4.  **Custom Reference Selection** (via a command-line flag).

5.  **Grouping** (optionally analyze different subsets/categories of
    isolates).

6.  **Summary Statistics** (aggregated in a final multi-section sheet).

The script can also produce **mutation matrices** (which track the exact
base or residue each isolate has at each position) and supports extra
features like generating a VCF file via snp-sites.

Table of Contents
-----------------

1.  [System Requirements](#system-requirements)

2.  [Installation](#installation)

3.  [Usage](#usage)

    -   3.1 [Command-Line Flags](#command-line-flags)

    -   3.2 [Nucleotide vs Protein vs
        Both](#nucleotide-vs-protein-vs-both)

    -   3.3 [Reference Selection](#reference-selection)

    -   3.4 [Grouping](#grouping)

    -   3.5 [Strict Validation Mode](#strict-validation-mode)

4.  [Input Files](#input-files)

5.  [Output Files](#output-files)

    -   5.1 [Analysis Sheets](#analysis-sheets)

    -   5.2 [Mutation Matrix Sheets](#mutation-matrix-sheets)

    -   5.3 [Summary Statistics](#summary-statistics)

6.  [Nucleotide Analysis Columns](#nucleotide-analysis-columns)

7.  [Protein Analysis Columns](#protein-analysis-columns)

8.  [How Codon-Based Nucleotide Logic
    Works](#how-codon-based-nucleotide-logic-works)

9.  [Examples](#examples)

10. [Frequently Asked Questions (FAQ)](#frequently-asked-questions)

11. [Contact / Issues](#contact--issues)

--------
# System Requirements

-   **Python 3.6+**

-   **Biopython** (for FASTA parsing and codon translation)

-   **pandas** and **XlsxWriter** (for Excel writing and data
    manipulation)

-   **MAFFT** if you do “both” (nucleotide + protein realignment) or if
    you choose -t p but your FASTA actually contains nucleotides (the
    script will translate and then re-align them).

-   **snp-sites** (optional) if you want the script to generate
    a .vcf file from the nucleotide alignment (`--vcf` flag).

--------
# Installation

1.  **Install Python 3** (if not installed).

2.  Install Python dependencies:

`pip install biopython pandas xlsxwriter`

1.  Install **MAFFT** (e.g., apt-get install mafft on Ubuntu, or from
    the MAFFT website).

2.  Optionally install **snp-sites** if you need `--vcf`.

--------
# Usage

Simply run:

```bash
python mutation\_analysis.py \

-a <your\_alignment.fasta> \

-t <n|p|both> \

[other flags...]
```

### Command-Line Flags

<table>
<thead>
<tr class="header">
<th><strong>Flag / Option</strong></th>
<th><strong>Meaning</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>-a, --alignment</td>
<td><strong>Required.</strong> Path to the aligned FASTA file.</td>
</tr>
<tr class="even">
<td>-t, --type</td>
<td><strong>Required.</strong> One of: <br />
<strong>n</strong> = Nucleotide analysis only <br />
<strong>p</strong> = Protein analysis only <br />
<strong>both</strong> = Nucleotide + Protein analysis</td>
</tr>
<tr class="odd">
<td>-f, --frame</td>
<td>Reading frame (1, 2, or 3) for nucleotide analysis. Default=1.</td>
</tr>
<tr class="even">
<td>-o, --output</td>
<td>Directory of output file</td>
</tr>
<tr class="odd">
<td>-d, --debug</td>
<td>Enables debug-level logging (messages will also appear in output.log).</td>
</tr>
<tr class="even">
<td>--vcf</td>
<td>If analyzing nucleotides, runs snp-sites to produce a .vcf file from your alignment.</td>
</tr>
<tr class="odd">
<td>--snp_sites_path</td>
<td>Path or command for snp-sites. Default: snp-sites.</td>
</tr>
<tr class="even">
<td>--mafft_path</td>
<td>Path or command for mafft. Default: mafft.</td>
</tr>
<tr class="odd">
<td>--temp</td>
<td>Keep intermediate/temporary files (e.g., the translated protein FASTA). By default, these files are deleted upon completion.</td>
</tr>
<tr class="even">
<td>--tmp_dir</td>
<td>Directory for temporary files. Default = current folder (./).</td>
</tr>
<tr class="odd">
<td>--job_id</td>
<td>Prefix for output files (helps to avoid overwriting). Default = output.</td>
</tr>
<tr class="even">
<td>--groups</td>
<td>CSV file [Isolate,Category] for grouping.</td>
</tr>
<tr class="odd">
<td>--reference</td>
<td>A specific isolate ID in the FASTA to treat as reference (otherwise, the first isolate in the file is used by default).</td>
</tr>
<tr class="even">
<td>--strict_validation</td>
<td>If set, the script checks if your input is “nucleotide” or “protein.” If you pick -t n/both but the input is protein, it errors out, and vice versa.</td>
</tr>
<tr class="odd">
<td>--quiet</td>
<td>If set, make terminal quiet unless the error happens. Recommended for autoruns</td>
</tr>
<tr class="even">
<td>--min_identity</td>
<td>Minimum identity threshold for alignment to improve the alignment quality. Default 80%</td>
</tr>
<tr class="odd">
<td>--min_coverage</td>
<td>Minimum coverage threshold for alignment to improve the alignment quality. Default 80%</td>
</tr>    
</tbody>
</table>

### Nucleotide vs Protein vs Both

-   **-t n (Nucleotide)**: The script checks codon changes (synonymous
    vs. non-synonymous).

-   **-t p (Protein)**: The script either:

    1.  If your alignment is **already protein**, it does direct protein
        analysis.

    2.  If your file is **actually nucleotides**, it automatically
        translates and realigns with MAFFT.

-   **-t both**: The script first does **Nucleotide** analysis, then
    translates → re-aligns (MAFFT) → **Protein** analysis.

### Reference Selection

By default, the **first** isolate in your FASTA becomes the reference.
If you want a specific reference:

`--reference some_isolate_id`

-   That isolate is placed first in the analysis, and all mutation calls
    are relative to its sequence.

### Grouping

If you provide `--groups groupfile.csv` with `[Isolate,Category]`, each
category is analyzed **separately**, producing separate sheets named
like:

-   `NucAnalysis(groupA) / NucMatrix(groupA)`

-   `NucAnalysis(groupB) / NucMatrix(groupB)`

-   etc.

Each category’s results appear in the final Excel (plus a combined
summary sheet). If an isolate isn’t listed in the CSV, it won’t be
included in any group.

### Strict Validation Mode

If you set `--strict_validation`, the code tries to guess whether your
alignment is “nucleotide” or “protein” (85% or more of non-gap
characters must be A/C/G/T/N for it to be considered “nucleotide”). If
your guess doesn’t match the requested analysis, the script terminates
with an error message.

--------------
# Input Files


1.  **FASTA Alignment** (Nucleotide or Protein). All sequences must have
    the **same length** (aligned).

2.  **Groups CSV** (optional):

> Isolate,Category
>
> iso1,groupA
>
> iso2,groupA
>
> iso3,groupB
>
> ...

1.  **Temporary/Intermediate** files are automatically created for
    certain steps (especially `both`, which realigns proteins).
    If --temp is **not** set, these are deleted at the end.

---------------
# Output Files

The main output is a **multi-sheet Excel** file. By default, it’s named:

`<job_id>_<type>_mutation_analysis.xlsx`

(e.g., `output_n_mutation_analysis.xlsx` if `-t n`). You can customize the name with `-o`.

### Analysis Sheets

-   **Nucleotide Analysis** (if `-t n` or `-t both`): Lists each position,
    reference NT, codon info, how many isolates have
    insertions/deletions/stop codons, etc.

-   **Protein Analysis** (if `-t p` or `-t both`): Lists each aligned amino acid position, reference AA, how many sequences have substitutions, insertions, or stops, etc.

If grouping is used, each group has its own analysis sheet with the suffix (`groupName`).

### Mutation Matrix Sheets

Alongside each Analysis sheet is a “Matrix” sheet
(e.g., **NucMatrix**, **ProtMatrix**) that displays:

-   Rows = isolate IDs

-   Columns = alignment positions

-   Cells = the observed base (or residue) if it differs from the
    reference.

If grouping is used, you’ll see separate matrix sheets for each group.

### Summary Statistics

Finally, a **Summary Statistics** sheet includes:

-   One table for **Nucleotide** results (if any).

-   Another table for **Protein** results (if any).

Each table shows total positions, how many were mutated, how many had
&gt;20% mutation frequency, total synonymous vs. non-syn changes (for
nucleotides), total insertions/deletions, etc.

------------------------------
# Nucleotide Analysis Columns

During **NucAnalysis**, each row represents a single **alignment
position**. Key columns:

1.  **Alignment Position**: 1-based index in the alignment.

2.  **Codon Number**: Which codon (1-based) this position belongs to, considering the `-f/--frame.`

3.  **Reference NT**: The reference nucleotide (- if a gap).

4.  **Mutations**: The mutated bases found among isolates
    (comma-separated).

5.  **Mutation Types**: The classification for each mutated base
    (Synonymous, Non-synonymous, Stop Codon, Insertion, Deletion, etc.).

6.  **Synonymous Mutations**: How many total times a base was changed
    but the resulting codon was the same amino acid.

7.  **Non-synonymous Mutations**: Changes that altered the amino acid.

8.  **Insertions**, **Deletions**, **Stop Codons**: The counts of each
    among the isolates at this position.

9.  **Total Mutations**: Sum of all mutation events (syn + non-syn +
    insertion + deletion + stop).

10. **Mutation Frequency**: `Total Mutations / Number_of_Isolates.`

---------------------------
# Protein Analysis Columns

During **ProtAnalysis**, each row represents a **protein alignment
position**:

1.  **Alignment Position**: 1-based index in the protein alignment.

2.  **Reference AA**: The reference amino acid (- if a gap).

3.  **Substitutions**: Which alternate amino acids appear
    (comma-separated).

4.  **Substitution Types**: A
    classification: Substitution, Insertion, Deletion, Stop Codon.

5.  **Substitution Count**: How many times a residue changed among
    isolates.

6.  **Insertions**, **Deletions**, **Stop Codons**: The counts at this
    position.

7.  **Total Mutations**: Summation of all changes.

8.  **Mutation Frequency**: `(Total Mutations) / (Number of Isolates that still have not encountered a stop).`

9.  **Remaining sequences before encountering stop codon**: The number of sequences that have not yet hit a stop codon in the alignment and thus “survive” to this position.

-----------------------------------------
# How Codon-Based Nucleotide Logic Works


When you run `-t n` or `-t both`, the code uses a **reading frame** (1, 2,or 3) to group nucleotides into codons. For a position pos:

1.  Compute `codon_index = (pos - (frame - 1)) // 3`.

2.  Extract the **3-nucleotide** chunk from the reference that begins at 

`ref_codon_start = (frame - 1) + codon_index * 3`.

3.  Translate the reference codon to an amino acid.

4.  For each isolate, extract the same codon region, translate to an
    amino acid.

5.  If the reference and sample amino acid differ,
    it’s **Non-synonymous**; if they match, it’s **Synonymous**.
    Insertions, deletions, or partial codons can lead to “Insertion,”
    “Deletion,” “Unknown,” etc.

**Important**: If your alignment includes leading/trailing UTR or
partial codons, check your --frame or ensure the aligned region starts
exactly at the coding region.

-----------
# Examples

### Example: Nucleotide Only

```bash
python mutation_analysis.py \

-a aligned_nuc.fasta \

-t n \

--reference REF_ISOLATE
```

Generates:

-   **NucAnalysis** sheet

-   **NucMatrix** sheet

-   **Summary Statistics**

### Example: Protein + Grouping

```bash
python mutation_analysis.py \

-a aligned_protein.fasta \

-t p \

--groups group.csv \

-o grouped_protein_analysis.xlsx
```

Produces separate **ProtAnalysis(groupX)** / **ProtMatrix(groupX)** sheets for each group, plus a summary.

### Example: Both

```bash
python mutation_analysis.py \

-a combined_nuc.fasta \

-t both \

--temp
```

1.  Nucleotide analysis.

2.  Translate + re-align with MAFFT.

3.  Protein analysis.

4.  Because `--temp` is set, it keeps intermediate FASTA files.

------------------------------
# Frequently Asked Questions

1.  **Q**: *Why do I see “Non-synonymous” but no difference in the
    actual amino acids?*  
    **A**: Possibly the reading frame is wrong, or partial codons are
    introduced by alignment gaps. Check `--frame` and ensure your
    alignment is correct.

2.  **Q**: *Why do I see `KeyError Worksheet not found` in older
    versions?*  
    **A**: The updated code ensures each new sheet is properly stored.
    You need at least `XlsxWriter ≥ 1.2`.

3.  **Q**: *Can I remove the “Codon Number” or “Mutations” columns in
    the NucAnalysis?*  
    **A**: Yes, but the code is structured so that these columns remain
    for clarity about codon-based changes.

4.  **Q**: *If I have only 1 isolate in a group, do I get meaningful
    stats?*  
    **A**: With just 1 isolate, mutation frequencies will all be 0 or 1
    for that isolate. We usually recommend 2+ isolates per group to see
    meaningful comparisons.

--------------------
# Contact / Issues

For questions, suggestions, or bug reports, please contact:

-   **Name**: Priyanshu Singh Raikwar, Seungwon Ko

-   **Email**: <priyanshu.raikwar@biology.ox.ac.uk>,
    <seungwon.ko@biology.ox.ac.uk>

Or open an issue in the repository .

**End of Document**
