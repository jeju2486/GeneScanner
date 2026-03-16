import os
import sys
import numpy as np
import argparse
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

##########################
# Command line arguments and parsing
##########################

parser = argparse.ArgumentParser(description='Visualize GeneScanner output')
parser.add_argument('-i', '--input', type=str, help='Path to the GeneScanner output Excel file')
parser.add_argument('-w', '--worksheet_name', type=str, default=None, help='Name of the worksheet in the Excel file that contains the analysis data (optional)')
parser.add_argument('-s', '--gene_starting_position', type=int, default=1, help='Starting position in the aligned sequences')
parser.add_argument('-nc', '--noncoding', action='store_true', help='Indicate if the nucleotide sequence is non-coding')
parser.add_argument('-r', '--reverse', action='store_true', help='Indicate if the nucleotide sequence is on the reverse strand')
parser.add_argument('-t', '--threshold', type=float, default=10.0, help='Mutation frequency threshold (percentage)')
parser.add_argument('-S', '--style', type=int, default=1, help='Predefined plotting styles for the output image file. Choose 1 or 2. If not specified, defaults to 1.')
parser.add_argument('-o', '--output', type=str, default=None, help='Directory path for output images')
args = parser.parse_args()

# Parse command line arguments
input_filepath = args.input
gene_starting_position = args.gene_starting_position
noncoding = args.noncoding
reverse = args.reverse
y_threshold = args.threshold
style = args.style

def save_figure(fig, output_dir, filename):
    """Save figure to specified directory"""
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, filename)
        fig.savefig(filepath)
        print(f"Saved: {filepath}")
    else:
        fig.show()

def process_worksheet(worksheet_name, input_filepath):
    """Process a single worksheet and generate plots"""
    # Identify data type
    data_type = 'nucleotide'
    if 'ProtAnalysis' in worksheet_name:
        data_type = 'protein'

    print(f"Processing worksheet: {worksheet_name} (Data type: {data_type})")

    df_raw = pd.read_excel(input_filepath, sheet_name=worksheet_name, header=None)
    number_of_isolates = int(df_raw.iloc[1][0].split(':')[-1].strip())
    df_raw = df_raw.drop(index=[0,1,2,3]).fillna(0).reset_index(drop=True) # Remove header lines 
    df_raw.columns = df_raw.iloc[0]
    df = df_raw.drop(index=0).reset_index(drop=True) # Final dataframe

    if data_type == 'nucleotide':
        # Add Non-coding Mutations column if gene is not coding
        if noncoding:
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
            dpi=100, 
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

        if style == 1:
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
        elif style == 2:
            line_color_per_mutation_type = {
                'Synonymous Mutations': 'k',
                'Non-synonymous Mutations': 'm',
                'Insertions': 'b',
                'Deletions': 'r',
                'Stop Codons': 'y'
            }
            if noncoding:
                mutation_types.append('Non-coding Mutations')
                mutation_types.remove('Synonymous Mutations')
                mutation_types.remove('Non-synonymous Mutations')
                mutation_types.remove('Stop Codons')

                line_color_per_mutation_type['Non-coding Mutations'] = 'k'
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
                    if style == 1:
                        ax1.scatter(X[i], y,
                                    facecolor=marker_color_per_mutation_type[mutation],
                                    marker=marker_per_mutation_type[mutation],
                                    s=marker_size_per_mutation_type[mutation],
                                    edgecolor=marker_edgecolor_per_mutation_type[mutation],
                                    label="",
                                    alpha=0.75,
                                    zorder=zorder_per_mutation_type[mutation])
                        
                        stem = ax1.stem(X[i], y, markerfmt=' ', linefmt='k-', basefmt=' ')
                        stem[1].set_linewidth(1)

                    elif style == 2:
                        stem = ax1.stem(X[i], y, markerfmt=' ', linefmt=f'{line_color_per_mutation_type[mutation]}-', basefmt=' ')
                        stem[1].set_linewidth(1.5)

        ax1.axhline(y=y_threshold, color='grey', linestyle='--', lw=1)
        ax1.text(len(df) + 5, y_threshold + 2, f'< {y_threshold}%', fontsize=15, rotation=90, color='grey')
        
        ax1.set_ylim(-1, 110)
        ax1.set_ylabel('Mutation frequency (%)', fontsize=17)
        ax1.tick_params(axis='y', labelsize=15)
        ax1.set_title(f'GeneScanner Mutational Analysis: {worksheet_name}', fontsize=20)
        ax1.grid(axis='y')

        mutation_types_no_syn = set(mutation_types) - set(['Synonymous Mutations', 'Non-coding Mutations'])
        if style == 1:
            legend_elements = [
                plt.Line2D([0], [0], marker=marker_per_mutation_type[mutation], color='w',
                            markerfacecolor=marker_color_per_mutation_type[mutation], markersize=15, markeredgecolor='black',label=mutation)
                for mutation in mutation_types_no_syn
            ]
            ax1.legend(handles=legend_elements, loc='upper center', fontsize=15, title='', ncol=4)

        elif style == 2:
            legend_elements = [
                plt.Line2D([0], [0], color=line_color_per_mutation_type[mutation], lw=2, label=mutation)
                for mutation in mutation_types_no_syn
            ]
            ax1.legend(handles=legend_elements, loc='upper center', fontsize=15, title='', ncol=4)

        ########################################
        # Plot low-frequency mutations
        ########################################
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

        if reverse:
            heatmap_data = heatmap_data.iloc[:, ::-1]

        if style == 1:
            heatmap_data_percentage = (heatmap_data.astype(float) / number_of_isolates * 100)
            heatmap_data_mapped = heatmap_data.copy().astype(float)
            heatmap_data_mapped = np.where(heatmap_data == 0, 0,
                                        np.where(heatmap_data_percentage < y_threshold, 1,
                                                2))
            custom_cmap = ListedColormap(['white', 'hotpink', 'black'])
            sns.heatmap(heatmap_data_mapped, cmap=custom_cmap, ax=ax2, cbar=False, alpha=1)

        elif style == 2:
            heatmap_data_bool = (heatmap_data.astype(float) > 0)
            sns.heatmap(heatmap_data_bool, cmap='binary', ax=ax2, cbar=False, alpha=0.5)

        ax2.set_xlabel('alignment position', fontsize=14)
        ax2.set_ylabel('', fontsize=17)
        ax2.tick_params(axis='both', which='both', length=0)
        
        x0 = gene_starting_position
        xticks_spacing = 100
        ax2.set_xticks(range(x0, df['Alignment Position'].max() + x0 + 1, xticks_spacing))
        ax2.set_xticklabels(range(x0, df['Alignment Position'].max() + x0 + 1, xticks_spacing), fontsize=12, rotation=90)

        ax2.set_yticks([i + 0.5 for i in range(len(mutation_order))])
        ax2.set_yticklabels([mutation_rename[m] for m in mutation_order], fontsize=14, rotation=0)

        ax2.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.75)
        ax2.set_xlim(-1, len(df)+1)

        ax2.plot([], [], color='hotpink', marker='s', linestyle='', markersize=15, label=f'Low-frequency mutations (< {y_threshold}%)')
        ax2.plot([], [], color='black', marker='s', linestyle='', markersize=15, label=f'High-frequency mutations (≥ {y_threshold}%)')

        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), fontsize=12, title='', title_fontsize=14)

        fig.tight_layout(h_pad=2.0)
        save_figure(fig, args.output, f'{worksheet_name}_nucleotide_analysis.png')
        plt.close(fig)

    elif data_type == 'protein':
        ######################################
        # Plot settings
        ######################################
        n_mutation_types = 6

        fig, (ax1, ax2) = plt.subplots(
            2, 1, 
            figsize=(20, 2*n_mutation_types),
            dpi=100, 
            sharex=True,
            gridspec_kw={'height_ratios': [1, 1]}
        )
        mutation_types = [
            'Substitution(Frameshift) Count',
            'Substitution Count',
            'Insertions',
            'Deletions',
            'Stop Codons',
            'Stop Codons(Frameshift)'
        ]

        if style == 1:
            marker_per_mutation_type = {
                'Substitution(Frameshift) Count': 's',
                'Substitution Count': 'D',
                'Insertions': 'P',
                'Deletions': 'o',
                'Stop Codons': '*',
                'Stop Codons(Frameshift)': '^'
            }
            marker_color_per_mutation_type = {
                'Substitution(Frameshift) Count': 'black',
                'Substitution Count': 'black',
                'Insertions': 'blue',
                'Deletions': 'red',
                'Stop Codons': 'yellow',
                'Stop Codons(Frameshift)': 'orange'
            }
            marker_edgecolor_per_mutation_type = {
                'Substitution(Frameshift) Count': 'black',
                'Substitution Count': 'white',
                'Insertions': 'blue',
                'Deletions': 'black',
                'Stop Codons': 'black',
                'Stop Codons(Frameshift)': 'black'
            }
            zorder_per_mutation_type = {
                'Substitution(Frameshift) Count': 0,
                'Substitution Count': 10,
                'Insertions': 20,
                'Deletions': 30,
                'Stop Codons': 40,
                'Stop Codons(Frameshift)': 50
            }
            marker_size_per_mutation_type = {
                'Substitution(Frameshift) Count': 0,
                'Substitution Count': 75,
                'Insertions': 100,
                'Deletions': 100,
                'Stop Codons': 250,
                'Stop Codons(Frameshift)': 150
            }

        elif style == 2:
            line_color_per_mutation_type = {
                'Substitution(Frameshift) Count': 'k',
                'Substitution Count': 'm',
                'Insertions': 'b',
                'Deletions': 'r',
                'Stop Codons': 'y',
                'Stop Codons(Frameshift)': 'g'
            }
        ######################################
        # Plot high-frequency mutations
        ######################################
        for mutation in mutation_types:
            X = df['Alignment Position'].astype(int) + gene_starting_position
            Y = df[mutation].astype(int)/number_of_isolates * 100

            for i, y in enumerate(Y):
                if y > y_threshold:
                    if style == 1:
                        ax1.scatter(X[i], y,
                                    facecolor=marker_color_per_mutation_type[mutation],
                                    marker=marker_per_mutation_type[mutation],
                                    s=marker_size_per_mutation_type[mutation],
                                    edgecolor=marker_edgecolor_per_mutation_type[mutation],
                                    label="",
                                    alpha=0.75,
                                    zorder=zorder_per_mutation_type[mutation])
                        
                        stem = ax1.stem(X[i], y, markerfmt=' ', linefmt='k-', basefmt=' ')
                        stem[1].set_linewidth(1)

                    elif style == 2:
                        stem = ax1.stem(X[i], y, markerfmt=' ', linefmt=f'{line_color_per_mutation_type[mutation]}-', basefmt=' ')
                        stem[1].set_linewidth(1.5)

        ax1.axhline(y=y_threshold, color='grey', linestyle='--', lw=1)
        ax1.text(len(df) + 5, y_threshold + 2, f'< {y_threshold}%', fontsize=15, rotation=90, color='grey')
        
        ax1.set_ylim(-1, 110)
        ax1.set_ylabel('Mutation frequency (%)', fontsize=17)
        ax1.tick_params(axis='y', labelsize=15)
        ax1.set_title(f'GeneScanner Mutational Analysis: {worksheet_name}', fontsize=20)
        ax1.grid(axis='y')

        mutation_types_no_syn = set(mutation_types)
        if style == 1:
            legend_elements = [
                plt.Line2D([0], [0], marker=marker_per_mutation_type[mutation], color='w',
                            markerfacecolor=marker_color_per_mutation_type[mutation], markersize=15, markeredgecolor='black',label=mutation)
                for mutation in mutation_types_no_syn
            ]
            ax1.legend(handles=legend_elements, loc='upper center', fontsize=15, title='', ncol=4)

        elif style == 2:
            legend_elements = [
                plt.Line2D([0], [0], color=line_color_per_mutation_type[mutation], lw=2, label=mutation)
                for mutation in mutation_types_no_syn
            ]
            ax1.legend(handles=legend_elements, loc='upper center', fontsize=15, title='', ncol=4)

        ########################################
        # Plot low-frequency mutations
        ########################################
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
        mutation_order = [
            'Insertions',
            'Deletions',
            'Substitution(Frameshift) Count',
            'Substitution Count',
            'Stop Codons',
            'Stop Codons(Frameshift)'
        ]
        mutation_rename = {
            'Substitution(Frameshift) Count': 'Substitution(FS)',
            'Substitution Count': 'Substitution',
            'Insertions': 'Insertion',
            'Deletions': 'Deletion',
            'Stop Codons': 'Stop Codon',
            'Stop Codons(Frameshift)': 'Stop Codon(FS)'
        }

        heatmap_data = heatmap_data.reindex(mutation_order)
        heatmap_data.rename(index=mutation_rename, inplace=True)

        if reverse:
            heatmap_data = heatmap_data.iloc[:, ::-1]

        if style == 1:
            heatmap_data_percentage = (heatmap_data.astype(float) / number_of_isolates * 100)
            heatmap_data_mapped = heatmap_data.copy().astype(float)
            heatmap_data_mapped = np.where(heatmap_data == 0, 0,
                                        np.where(heatmap_data_percentage < y_threshold, 1,
                                                2))
            if heatmap_data_mapped.max() == 1:
                custom_cmap = ListedColormap(['white', 'hotpink'])
            else:
                custom_cmap = ListedColormap(['white', 'hotpink', 'black'])
            sns.heatmap(heatmap_data_mapped, cmap=custom_cmap, ax=ax2, cbar=False, alpha=1)

        elif style == 2:
            heatmap_data_bool = (heatmap_data.astype(float) > 0)
            sns.heatmap(heatmap_data_bool, cmap='binary', ax=ax2, cbar=False, alpha=0.5)

        ax2.set_xlabel('alignment position', fontsize=14)
        ax2.set_ylabel('', fontsize=17)
        ax2.tick_params(axis='both', which='both', length=0)
        
        x0 = gene_starting_position
        xticks_spacing = 100
        ax2.set_xticks(range(x0, df['Alignment Position'].max() + x0 + 1, xticks_spacing))
        ax2.set_xticklabels(range(x0, df['Alignment Position'].max() + x0 + 1, xticks_spacing), fontsize=12, rotation=90)

        ax2.set_yticks([i + 0.5 for i in range(len(mutation_order))])
        ax2.set_yticklabels([mutation_rename[m] for m in mutation_order], fontsize=14, rotation=0)

        ax2.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.75)
        ax2.set_xlim(-1, len(df)+1)

        ax2.plot([], [], color='hotpink', marker='s', linestyle='', markersize=15, label=f'Low-frequency mutations (< {y_threshold}%)')
        ax2.plot([], [], color='black', marker='s', linestyle='', markersize=15, label=f'High-frequency mutations (≥ {y_threshold}%)')

        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), fontsize=12, title='', title_fontsize=14)

        fig.tight_layout(h_pad=2.0)
        save_figure(fig, args.output, f'{worksheet_name}_protein_analysis.png')
        plt.close(fig)

def main():
    #########################################
    # Pre-process Excel file to get a clean DataFrame
    #########################################
    # Read available sheet names
    available_sheets = pd.ExcelFile(input_filepath).sheet_names
    
    # Determine which worksheets to process
    if args.worksheet_name:
        worksheets_to_process = [args.worksheet_name]
        if args.worksheet_name not in available_sheets:
            print(f"ERROR: Worksheet '{args.worksheet_name}' not found in the Excel file.")
            print(f'Available sheets: {", ".join(available_sheets)}')
            sys.exit(1)
    else:
        # Process all NucAnalysis and ProtAnalysis worksheets
        worksheets_to_process = [ws for ws in available_sheets if 'NucAnalysis' in ws or 'ProtAnalysis' in ws]
        if not worksheets_to_process:
            print("ERROR: No worksheets matching 'NucAnalysis' or 'ProtAnalysis' found.")
            print(f'Available sheets: {", ".join(available_sheets)}')
            sys.exit(1)
    
    # Process each worksheet
    for worksheet in worksheets_to_process:
        try:
            process_worksheet(worksheet, input_filepath)
        except Exception as e:
            print(f"ERROR processing worksheet '{worksheet}': {e}")
            continue

if __name__ == "__main__":
    main()
