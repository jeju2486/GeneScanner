# caro_project

## Introduction

The purpose of this program is **One-click** user-friendly tool to calculate mutation locations and rates from an alignment file. Of course, there are a lot of more sophisticated tools, this tool targets the beginner bioinformaticians know minimum of programming.

The script is written in Python and supports analysis of nucleotide, protein, or both types of alignments.

## Dependencies

- **Biopython** v1.84 or higher  
- **xlsxwriter** v3.2.0 or higher  
- **Python** 3.6 or higher
- **snp-sites** v2.5.1 or higher  (Optional; for vcf file generating)
- **MAFFT** v7.490 or higher  (Optional; for protein analysis using nucleotide)

## How to Generate the Input File

1. Visit [PubMLST](https://pubmlst.org).
2. Search for your organism of interest (e.g. *Staphylococcus aureus*).
3. Select **Genome Collection**.
4. Choose genomes by filtering (e.g. IDs ≤ 6202, total length ≥ 2 Mbp).
5. Scroll to the bottom and select `Third party` → `SNPsites`.
6. Choose the loci (e.g. *rpsD*).
7. Submit the query, wait for the results, and save the **alignment file**.

### Input File

The input file should be a FASTA **alignment file**.

## How to Run the Code

The script is named `mutation_analysis.py` and accepts several command-line parameters. Below are some common examples.

### For Nucleotide Alignment Analysis

```bash
python mutation_analysis.py -a "nucleotide_alignment.fasta" -t n
```

### For Protein Alignment Analysis

```bash
python mutation_analysis.py -a "protein_alignment.fasta" -t p
```

### For Both Analyses (Nucleotide and Re-aligned Protein)

In this mode, the program first performs nucleotide analysis on the provided alignment. Then it extracts the original (ungapped) nucleotide sequences, pads them (if needed) to ensure their length is a multiple of three, translates them into proteins, realigns them using MAFFT, and performs protein analysis.

```bash
python mutation_analysis.py -a "nucleotide_alignment.fasta" -t both
```

### Optional Parameters

- **Specify the Reading Frame** (default is 1):
  
  ```bash
  python mutation_analysis.py -a "alignment.fasta" -t p -f 2
  ```

- **Specify the Output File** (default names are automatically generated based on the alignment type):
  
  ```bash
  python mutation_analysis.py -a "alignment.fasta" -t p -o "custom_output.xlsx"
  ```

- **Enable Debug Mode**  
  (When enabled, debug messages are printed to the console and saved to `output.log`):
  
  ```bash
  python mutation_analysis.py -a "alignment.fasta" -t p -d
  ```

- **Generate VCF File**
  (Only valid for nucleotide and both modes):
  
  ```bash
  python mutation_analysis.py -a "alignment.fasta" -t n --vcf
  ```

please use help command to see more details

```bash
python mutation_analysis.py --help
```

## Output

The script generates an Excel workbook containing multiple sheets:

- **Nucleotide Analysis**:  
  Shows the reference sequence (the first sequence in the alignment file), alignment position, codon number, codon position, mutation types, mutation frequency, and a new column **Mutated AA(s)** that lists the resulting amino acids for mutations (positioned right after the **Mutations** column).

- **Nucleotide Mutation Matrix / Protein Mutation Matrix**:  
  Displays which isolate has a mutation at which position and, for protein analysis, what amino acid the mutation results in.

- **Summary Statistics**:  
  Provides an overview of mutation rates, including the total positions with mutations, percentages, and additional summary metrics.

## Additional Notes

- In the **translation function**, before translating, each nucleotide sequence is padded with "N" (if necessary) so its length becomes a multiple of three. This ensures that any incomplete codon at the end of the sequence is handled correctly (typically resulting in an unknown amino acid, "X").
- In **debug mode** (`-d`), log messages are both printed to the console and saved to an `output.log` file.
- The program supports realigning protein sequences using MAFFT when running in "both" mode.
- The VCF file generation via snp-sites is optional and only applicable when analysing nucleotide sequences.

## TODO

1. Implement additional visualization tools.
2. Optimize the code for performance.
3. Compare results with alternative methods (e.g., Billy’s results) for validation.
4. Expand and refine the protein alignment analysis (currently, it supports both nucleotide and protein analyses with translation).
5. Investigate integration with tools like Jalview, which offers similar functionalities.
6. Make the reference file selectable
