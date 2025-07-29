# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 13:12:53 2025

@author: cariefrantz
"""


####################
# IMPORTS
####################

import os
import pandas as pd
import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'Arial'

import RottenIceModules
import RottenIceVars

####################
# VARIABLES
####################
out_dir = os.getcwd() # Output directory for saving the figure and csv

datasets = ['16S','18S','18S-PrimProd']
color_map = RottenIceVars.plotColorsByMonth
month_order = ['M','JN','JY10','JY11']

####################
# FUNCTIONS
####################

def parse_sample_id(sample_id):
    try:
        sample, gene_template = sample_id.split('.')
        # Get sample metadata
        sample_parts = sample.split('-')
        month = sample_parts[0]
        replicate = sample_parts[-1]
        material = sample_parts[-2]
        # Get sequencing info
        template_parts = gene_template.split('-')
        gene = template_parts[0]
        template = template_parts[1]
        return {
            'month': month,
            'material': material,
            'replicate': replicate,
            'gene': gene,
            'template': template
        }
    except Exception:
        return None



def get_group_key(parsed_id):
    """Returns the grouping key (e.g., M-SW)"""
    return f"{parsed_id['month']}-{parsed_id['material']}"

def get_replicate_key(parsed_id):
    """Returns the replicate key (e.g., M-CS-SW-1.16S)"""
    return f"{parsed_id['month']}-{parsed_id['material']}-{parsed_id['replicate']}.{parsed_id['gene']}"

def extract_distances(dist_df):
    parsed_ids = {sid: parse_sample_id(sid) for sid in dist_df.index}
    valid_ids = {sid: pid for sid, pid in parsed_ids.items() if pid is not None}

    # Build dictionary by replicate key
    replicate_pairs = {}
    for sid, pid in valid_ids.items():
        rep_key = get_replicate_key(pid)
        replicate_pairs.setdefault(rep_key, {})[pid['template']] = sid

    results = []
    for rep_key, templates in replicate_pairs.items():
        if 'DNA' in templates and 'cDNA' in templates:
            dna = templates['DNA']
            cdna = templates['cDNA']
            parsed = parse_sample_id(dna)
            group_key = get_group_key(parsed)
            distance = dist_df.loc[dna, cdna]
            results.append(
                {'group': group_key, 'replicate': rep_key,
                 'distance': distance, 'month': parsed['month']}
                )
    
    return pd.DataFrame(results)

def plot_boxplots(df, title, directory):
    # Compute average distance per group (e.g., M-CS-SW.16S)
    group_means = df.groupby('group')['distance'].mean().sort_values(ascending=False)
    df['group'] = pd.Categorical(df['group'], categories=group_means.index, ordered=True)

    plt.figure(figsize=(12, 6))
    sns.boxplot(
        data=df, x='group', y='distance', width = 0.8,
        hue='month', palette = color_map, dodge=False)
    plt.xticks(rotation=90)
    plt.xlabel('Sample set')
    plt.ylabel('DNA-cDNA Weighted UniFrac Distance')
    plt.title(title)
    
    # Sort the legend entries
    handles, labels = plt.gca().get_legend_handles_labels()
    label_handle_map = dict(zip(labels, handles))
    sorted_handles = [label_handle_map[m] for m in month_order if m in label_handle_map]
    sorted_labels = [m for m in month_order if m in label_handle_map]
    plt.legend(sorted_handles, sorted_labels, frameon=False, title='Month')
    
    # Display the plot
    plt.tight_layout()
    plt.show()
    
    plt.savefig(directory + '\\'+'DNA-cDNA-distplot_' + title + '.svg')


#%%

####################
# MAIN FUNCTION
####################

for ds in datasets:
    # Import metadata
    filename, directory, distmat = RottenIceModules.fileGet(
        'Select ' + ds + ' distance matrix', file_type = 'tsv')
    df = extract_distances(distmat)
    df.to_csv(out_dir + '\\' + 'DNA-cDNA-distmat_' + ds + '.csv', index=False)
    plot_boxplots(df, ds, out_dir)
