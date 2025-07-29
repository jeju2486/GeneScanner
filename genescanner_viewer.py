import os
import sys
import numpy as np
import argparse
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

##########################
# Command line arguments and parsing
##########################

parser = argparse.ArgumentParser(description='Visualize GeneScanner output')
parser.add_argument('-i', '--input', type=str, help='Path to the GeneScanner output Excel file')
parser.add_argument('-s', '--gene_starting_position', type=int, default=1, help='Starting position in the aligned sequences')
parser.add_argument('-nc', '--noncoding', action='store_true', help='Indicate if the nucleotide sequence is non-coding')
parser.add_argument('-r', '--reverse', action='store_true', help='Indicate if the nucleotide sequence is on the reverse strand')
parser.add_argument('-t', '--threshold', type=float, default=10.0, help='Mutation frequency threshold (percentage)')
parser.add_argument('-o', '--output', type=str, default=None, help='Path for output image file, otherwise, shown interactively)')
args = parser.parse_args()

input_filepath = args.input
gene_starting_position = args.gene_starting_position
noncoding = args.noncoding
reverse = args.reverse
y_threshold = args.threshold
if args.output:
    output_filepath = args.output

def main():
    #########################################
    # Pre-process Excel file to get a clean DataFrame
    #########################################
    df_raw = pd.read_excel(input_filepath, sheet_name='NucAnalysis', header=None)
    number_of_isolates = int(df_raw.iloc[1][0].split(':')[-1].strip())
    df_raw = df_raw.drop(index=[0,1,2,3]).fillna(0).reset_index(drop=True) # Remove header lines 
    df_raw.columns = df_raw.iloc[0]
    df = df_raw.drop(index=0).reset_index(drop=True) # Final dataframe

    # Add Non-coding Mutations column if gene is not coding
    if noncoding:
        print("Hello :)")
        df['Non-coding Mutations'] = 0

        # Lump coding mutations into non-coding mutations
        A = df['Non-synonymous Mutations']
        B = df['Synonymous Mutations']
        C = df['Stop Codons']
        df['Non-coding Mutations'] = A + B + C

        # Set coding mutations to zero
        df['Non-synonymous Mutations'] = 0
        df['Synonymous Mutations'] = 0
        df['Stop Codons'] = 0
    ######################################
    # Plot settings
    ######################################
    n_mutation_types = 5
    if noncoding:
        n_mutation_types = 3

    fig, (ax1, ax2) = plt.subplots(
        2, 1, 
        figsize=(20, 2*n_mutation_types),
        dpi=200, 
        sharex=True,
        gridspec_kw={'height_ratios': [1, 1]}
    )
    mutation_types = [
        'Synonymous Mutations',
        'Non-synonymous Mutations',
        'Insertions',
        'Deletions',
        'Stop Codons'
    ]
    marker_per_mutation_type = {
        'Synonymous Mutations': 's',
        'Non-synonymous Mutations': 'D',
        'Insertions': 'P',
        'Deletions': 'o',
        'Stop Codons': '*'
    }
    marker_color_per_mutation_type = {
        'Synonymous Mutations': 'black',
        'Non-synonymous Mutations': 'black',
        'Insertions': 'blue',
        'Deletions': 'red',
        'Stop Codons': 'yellow'
    }
    marker_edgecolor_per_mutation_type = {
        'Synonymous Mutations': 'black',
        'Non-synonymous Mutations': 'white',
        'Insertions': 'blue',
        'Deletions': 'black',
        'Stop Codons': 'black'
    }
    zorder_per_mutation_type = {
        'Synonymous Mutations': 0,
        'Non-synonymous Mutations': 10,
        'Insertions': 20,
        'Deletions': 30,
        'Stop Codons': 40
    }

    if noncoding:
        mutation_types.append('Non-coding Mutations')
        mutation_types.remove('Synonymous Mutations')
        mutation_types.remove('Non-synonymous Mutations')
        mutation_types.remove('Stop Codons')

        marker_per_mutation_type['Non-coding Mutations'] = 'X'
        marker_color_per_mutation_type['Non-coding Mutations'] = 'green'
        marker_edgecolor_per_mutation_type['Non-coding Mutations'] = 'black'
        zorder_per_mutation_type['Non-coding Mutations'] = 50

    # Marker sizes
    if noncoding:
        marker_size_per_mutation_type = {
            'Synonymous Mutations': 0,
            'Non-synonymous Mutations': 0,
            'Insertions': 100,
            'Deletions': 100,
            'Stop Codons': 0,
            'Non-coding Mutations': 0
        }
    else:
        marker_size_per_mutation_type = {
            'Synonymous Mutations': 0,
            'Non-synonymous Mutations': 75,
            'Insertions': 100,
            'Deletions': 100,
            'Stop Codons': 250
        }
    ######################################
    # Plot high-frequency mutations
    ######################################
    for mutation in mutation_types:
        X = df['Alignment Position'].astype(int) + gene_starting_position
        Y = df[mutation].astype(int)/number_of_isolates * 100

        if reverse:
            Y = Y.values[::-1]

        for i, y in enumerate(Y):
            if y > y_threshold:
                # Plot scatter points
                ax1.scatter(X[i], y,
                            facecolor=marker_color_per_mutation_type[mutation],
                            marker=marker_per_mutation_type[mutation],
                            s=marker_size_per_mutation_type[mutation],
                            edgecolor=marker_edgecolor_per_mutation_type[mutation],
                            label="",
                            alpha=0.75,
                            zorder=zorder_per_mutation_type[mutation])
                
                # Plot vertical lines
                stem = ax1.stem(X[i], y, markerfmt=' ', linefmt='k-', basefmt=' ')
                stem[1].set_linewidth(1)

    # Plot horizontal line for frequency threshold
    ax1.axhline(y=y_threshold, color='grey', linestyle='--', lw=1)
    ax1.text(len(df) + 5, y_threshold + 2, f'< {y_threshold}%', fontsize=15, rotation=90, color='grey')
    
    # Custom plot settings
    ax1.set_xlim(-10, len(df)+10)
    ax1.set_ylim(-1, 110)
    ax1.set_ylabel('Mutation frequency (%)', fontsize=17)
    ax1.tick_params(axis='y', labelsize=15)
    ax1.set_title(f'GeneScanner Mutational Analysis', fontsize=20)
    ax1.grid(axis='y')

    # Legend box
    mutation_types_no_syn = set(mutation_types) - set(['Synonymous Mutations', 'Non-coding Mutations'])
    legend_elements = [
        plt.Line2D([0], [0], marker=marker_per_mutation_type[mutation], color='w',
                    markerfacecolor=marker_color_per_mutation_type[mutation], markersize=15, markeredgecolor='black',label=mutation)
        for mutation in mutation_types_no_syn
    ]
    ax1.legend(handles=legend_elements, loc='upper center', fontsize=15, title='', ncol=4)

    ########################################
    # Plot low-frequency mutations
    ########################################
    # Plot heatmap
    df_gene = df[['Alignment Position'] + mutation_types].copy()
    df_gene['Gene'] = 'gene'
    df_gene['Alignment Position'] = df_gene['Alignment Position'].astype(int) + gene_starting_position
    melted = df_gene.melt(id_vars=['Gene', 'Alignment Position'], var_name='Mutation Type', value_name='Count')

    heatmap_data = melted.pivot_table(
        index=['Mutation Type'],
        columns='Alignment Position',
        values='Count',
        fill_value=0
    )
    # Reorder and rename mutation types for heatmap
    mutation_order = [
        'Insertions',
        'Deletions',
        'Synonymous Mutations',
        'Non-synonymous Mutations',
        'Stop Codons'
    ]
    mutation_rename = {
        'Synonymous Mutations': 'Synonymous',
        'Non-synonymous Mutations': 'Non-synonymous',
        'Insertions': 'Insertion',
        'Deletions': 'Deletion',
        'Stop Codons': 'Stop Codon'
    }
    if noncoding:
        mutation_order.append('Non-coding Mutations')
        mutation_order.remove('Synonymous Mutations')
        mutation_order.remove('Non-synonymous Mutations')
        mutation_order.remove('Stop Codons')

        mutation_rename['Non-coding Mutations'] = 'Non-coding'

    heatmap_data = heatmap_data.reindex(mutation_order)
    heatmap_data.rename(index=mutation_rename, inplace=True)

    # Plot binary heatmap
    # heatmap_data_bool = (heatmap_data.astype(float) > 0)
    # heatmap_data_bool = (heatmap_data.astype(float) > 0) & ((heatmap_data.astype(float)/number_of_isolates * 100) < y_threshold)
    # sns.heatmap(heatmap_data_bool, cmap='binary', ax=ax2, cbar=False, alpha=1)


    heatmap_data_percentage = (heatmap_data.astype(float) / number_of_isolates * 100)
    # Create the custom mapped data
    heatmap_data_mapped = heatmap_data.copy().astype(float)
    # Apply conditions
    heatmap_data_mapped = np.where(heatmap_data == 0, 0,  # Zero values stay 0
                                  np.where(heatmap_data_percentage < y_threshold, 1,  # Above 0 but below threshold = 1
                                          2))  # Above threshold = 2
    # Create custom colormap - you can choose different color combinations:
    custom_cmap = ListedColormap(['white', 'hotpink', 'black'])
    sns.heatmap(heatmap_data_mapped, cmap=custom_cmap, ax=ax2, cbar=False, alpha=1)

    # Custom heatmap settings
    ax2.set_xlabel('alignment position', fontsize=14)
    ax2.set_ylabel('', fontsize=17)
    ax2.tick_params(axis='both', which='both', length=0)
    
    ## Set x-ticks
    x0 = gene_starting_position
    xticks_spacing = 100
    ax2.set_xticks(range(x0, df['Alignment Position'].max() + x0 + 1, xticks_spacing))
    ax2.set_xticklabels(range(x0, df['Alignment Position'].max() + x0 + 1, xticks_spacing), fontsize=12, rotation=90)

    ## Set y-ticks
    ax2.set_yticks([i + 0.5 for i in range(len(mutation_order))])
    ax2.set_yticklabels([mutation_rename[m] for m in mutation_order], fontsize=14, rotation=0)

    ## Add grid lines for better readability
    ax2.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.75)

    fig.tight_layout()
    if args.output:
        plt.savefig(args.output)
    else:
        # Show the plot interactively
        plt.show()

if __name__ == "__main__":
    main()