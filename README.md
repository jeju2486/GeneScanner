# caro_project

## Introduction

Purpose of the program is calculating the mutation location and rate. It is written in Python

## how to run the program

To run the program there is two inputs

## Dependencies

Biopython v.1.84 or higher
xlsxwriter v.3.2.0 or higher

### how to generate the inputfile
1. Enter PubMLST
2. Search the Organism Name (e.g. *Staphylococcus aureus*)
3. Select Genome Collection
4. Select the ID (e.g. id <= 6202, total length >= 2mbp)
5. From the bottom, slect `Third party` -> `SNPsites`
6. Select the Loci (e.g. rpsD)
7. Select the submit and wait for the result. Save the **alignment file** 

### input file 
**`alingment file`**

### how to run the code
**For Nucleotide Alignment:**
 
```ruby
python mutation_analysis.py -a "nucleotide_alignment.fasta" -t "n"
```

**For Protein Alignment:**
 
```ruby
python mutation_analysis.py -a "protein_alignment.fasta" -t "p"
```

**Specify reading frame (Optional):**
 
```ruby
python mutation_analysis.py -a "alignment.fasta" -t "p" -f 2 -o "custom_output.tsv"
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
* The script generates a three page excel sheets
* **Nucleotide Analysis** : Shows the reference file (Here it is first sequence of alignment file), Alignment position, Codon number and its mutation type and frequency.
* **Nucleotide Muataion Matrix**: Shows which isolate and where in that isolate has mutation and which amino acide it is changed to (Alreday exists in Snpsite plug-in).
* **Summary Statistics**: Summary of mutation rates

## TODO 
1. visualisation (or not)
2. Code optimization
3. Compare with the Billy's result to validate if it works fine
4. Add the function of protein alingment (which would be different with nuclotide alignment and translate). Make user to choose depends on what they want.
5. There is a progrma called 'Jalview' which already published and can do almost what we want :(