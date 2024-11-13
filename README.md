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
7. Select the submit and wait for the result. Save the **alignment file** and **vcf file**

* input 1 **`vcf file`**
* input 2 **`alingment file`**

how to run the code

```ruby
python mutation_analysis.py -v "$vcf_file" -a "alignment_file"
```