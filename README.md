# caro_project

## Introduction

Purpose of the program is calculating the mutation location and rate. It is written in Python

## how to run the program

To run the program there is two inputs

### how to generate the inputfile
1. Enter PubMLST
2. Search the Organism Name (e.g. *Staphylococcus aureus*)
3. Select Genome Collection
4. Select the ID (e.g. id <= 6202, total length >= 2mbp)
5. From the bottom, slect `Third party` -> `SNPsites`
6. Select the Loci (e.g. rpsD)
7. Select the submit and wait for the result. Save the **alignment file** 

* input file **`alingment file`**

### how to run the code
**For Nucleotide Alignment:**
 
```ruby
python mutation_analysis.py -a "nucleotide_alignment.fasta" -t "n"
```

**For Protein Alignment:**
 
```ruby
python mutation_analysis.py -a "protein_alignment.fasta" -t "p"
```

**Specify Output File (Optional):**
 
```ruby
python mutation_analysis.py -a "alignment.fasta" -t "p" -o "custom_output.tsv"
```

**Enable Debug Mode (Optional):**
 
```ruby
python mutation_analysis.py -a "alignment.fasta" -t "p" -d
```

## Output
* The script generates a tab-delimited file with mutation analysis results.
Default output file names:
* **Protein Alignment:** `protein_mutation_analysis.tsv`
* **Nucleotide Alignment:** `nucleotide_mutation_analysis.tsv`