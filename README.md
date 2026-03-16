# GeneScanner

![GeneScanner Logo](logo.png)

GeneScanner is a high-throughput mutation analysis tool for aligned nucleotide/protein FASTA data. This pipeline analyzes either **nucleotide** or **protein** aligned FASTA files (or both), identifies and classifies mutations, produces summary statistics, and writes a **multi-sheet Excel** report. 

[![DOI](https://zenodo.org/badge/997396273.svg)](https://doi.org/10.5281/zenodo.17495646)

---

## Table of Contents
1. [GeneScanner](#genescanner)
2. [Description](#description)
3. [System Requirements and Dependencies](#system-requirements-and-dependencies)
4. [Installation](#installation)
5. [How to Use](#how-to-use)

---

## GeneScanner

### Description

**GeneScanner** is a mutation analysis pipeline for aligned sequences.  
It reads **nucleotide**, **protein**, or **both** types of FASTA alignments and produces:

- Per-position mutation summaries with synonymous vs. non-synonymous classification 
- Detection of protein mutations types and location
- Group-wise comparisons 
  
---

## System Requirements and Dependencies
- **Python ≥ 3.6**
    - biopython ≥ 1.81 
    - numpy ≥ 1.23
    - pandas ≥ 1.5
    - XlsxWriter ≥ 3.0
- [**MAFFT**](https://github.com/GSLBiotech/mafft) ≥ 7.526  
- [**snp-sites**](https://github.com/sanger-pathogens/snp-sites) ≥  2.4.0

## Installation

There two main options

1. Directly from source 

Clone this repo
   ```
   git clone https://github.com/jeju2486/GeneScanner.git
   ```

Install
   ```
   pip install GeneScanner/
   ```

This will install all the Python dependencies.

2. Install conda environment (Recommended)

You can also install Python dependencies via conda with the `environment.yml` file.
This will avoid any potential version dependencies.

   ```
   conda env create -f environment.yml
   ```

This will create the environmenet named `GeneScanner-env` with all necessary dependency conflicts.

   ```
   conda activate GeneScanner-env
   ```

Then, you can install GeneScanner from source via pip as done above.

```
   pip install GeneScanner/
```

## How to Use

Full documentation is maintained in the **[Project Wiki](https://github.com/jeju2486/GeneScanner/wiki/GeneScanner)**

For questions or issues, please open an Issue in this repository or contact:

* Priyanshu Singh Raikwar ([priyanshu.raikwar@biology.ox.ac.uk](mailto:priyanshu.raikwar@biology.ox.ac.uk))
* Seungwon Ko ([seungwon.ko@biology.ox.ac.uk](mailto:seungwon.ko@biology.ox.ac.uk))

For installation or viewer issues, please contact to:

* Broncio Aguilar-Sanjuan ([broncio.aguilarsanjuan@biology.ox.ac.uk](mailto:broncio.aguilarsanjuan@biology.ox.ac.uk))
