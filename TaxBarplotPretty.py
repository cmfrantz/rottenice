# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 16:38:21 2025

@author: cariefrantz
@project: RottenIce

BUILDS TAXONOMIC ABUNDANCE BARPLOTS FORMATTED FOR PUBLICATION
This script builds stacked horizontal barplots of relative taxonomic
abundance for each sample. The plots are "prettified" for publication.

This script was created as part of the Rotten Ice Project.
ChatGPT was used to figure out how to adjust the formatting.


Arguments:  None

Requirements:
    Annotated, taxonomy-collapsed taxonomy table in csv format
        where rows = taxonomy calls
        columns = taxonomy, color, then a list of sample names
        Column 'color' defines the colors to use for each taxonomy value.
        Columns ranking, max, and avg are ignored
            (the full list of columns to ignore should be added to the
             script variables)
        The sample columns *must* be ordered how they should be ordered
        in the plot, sorted by months with the months in the correct order.

Example in command line:
    python TaxBarplotPretty.py

Dependencies Install:
    pip install numpy
    pip install pandas
    pip install matplotlib

You will also need to have the following files
    in the same directory as this script.
They contain modules and variables that this script calls.
    RottenIceModules.py
If you get an error indicating that one of these modules is not found,
    change the working directory to the directory containing these files.

Copyright (C) 2025  Carie M. Frantz

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
"""
####################
# IMPORTS
####################
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import RottenIceModules


####################
# VARIABLES
####################

# Data to plot
#datasets = ['16S noUEMC QC', '18S noUBA QC'] # List of the datasets to generate plots for
datasets = ['16S noUEMC QC']
#datasets = ['18S noUBA QC']

# Define the sample list to plot. If False, plots all samples.
samplelist = False
#samplelist = RottenIceModules.genSampleList(['M','JN','JY10','JY11'],['HT','HM','HB'])
# Set whether or not to average replicates of the samples
# If False, all replicates are shown (no averaging)
average_replicates = False
# If True, samples are averaged for a single displayed result
#average_replicates = True

# Define the template to include in the plot (can be 'DNA', 'cDNA', or 'both')
template = 'both'
#template = 'DNA'
#template = 'cDNA'

# Set whether or not to average the samp

# Spreadsheet variables
colorcol = 'color'  # The column header in the data table that defines the color hex value
ignorecols = ['Phylum','phylum','Supergroup','supergroup',
              'Big','PrimProd','Primary Producer',
              'color','ranking','rank','max','avg']     # Non-data columns that could be in the data table

# Plot details
fontsize = 11 # Font size for the plot
title = 'Relative abundance of major taxa' # Title header for the plot
monthmap = {
    'M'     : 'May',
    'JN'    : 'June',
    'JY10'  : 'July 10 dirty floe',
    'JY11'  : 'July 11 clean floe'
    }



####################
# FUNCTIONS
####################

def filter_sample_ids(sample_ids, samplelist, template='both'):
    """
    Filters sample IDs by matching sample names and template.

    Parameters
    ----------
    sample_ids : list of str
        List of full sample ID strings (e.g., 'M-Drain-1.16S-cDNA').
    samplelist : list of str
        List of sample names (e.g., ['M-CS-HT','M-CS-HM','M-Drain']) to match
        at the beginning of each sample ID.
    template : str
        Template to match (either 'DNA','cDNA', or 'both').

    Returns
    -------
    list of str
        Filtered sample IDs that match both the template and the specified site list.
    """
    result = []
    for sid in sample_ids:
        try:
            ssplit = sid.split('.')
            sample = ssplit[0].rsplit('-', 1)[0]  # Get site name (e.g., 'M-Drain')
            s_template = ssplit[1].split('-')[1]  # Get the template
            if sample in samplelist:
                if template == 'both':
                    result.append(sid)
                elif template in ['DNA','cDNA'] and s_template == template:
                    result.append(sid)
        except IndexError:
            continue  # Skip malformed IDs
    return result


def averageSamples(datatable, samples_parsed):
    '''
    Averages all replicates of each sample to be plotted.

    Parameters
    ----------
    datatable : pd.DataFrame
        Table with taxa as rows and sample names as columns
    samples_parsed : pd.DataFrame
        DataFrame with index = sample name and columns = 'month', 'label', 'template'

    Returns
    -------
    datatable_avg : pd.DataFrame
        Table with averaged replicates (columns are unique sample bases)
    samples_parsed_avg : pd.DataFrame
        Corresponding metadata, one row per averaged sample (replicate removed from label)
    '''

    # Make a copy to avoid modifying the original
    samples_parsed = samples_parsed.copy()

    # Remove replicate from 'label'
    samples_parsed['label'] = samples_parsed['label'].astype(str).str.split('-').str[:-1].str.join('-')

    # Combine parts into a unique sample base name
    samples_parsed['sample_base'] = (
        samples_parsed['month'].astype(str) + '-' +
        samples_parsed['label'] + '.' +
        samples_parsed['template'].astype(str)
    )

    # Transpose the data so samples are rows
    datatable_T = datatable.T.copy()
    datatable_T['sample_base'] = samples_parsed['sample_base']

    # Average across replicates
    averaged = datatable_T.groupby('sample_base').mean()

    # Transpose back to original layout
    datatable_avg = averaged.T

    # Build metadata table, matching the new averaged samples
    samples_parsed_avg = (
        samples_parsed
        .drop_duplicates('sample_base')
        .set_index('sample_base')
        .loc[datatable_avg.columns]
    )

    return datatable_avg, samples_parsed_avg


#%%
####################
# MAIN FUNCTION
####################
directory = os.getcwd()

for d in datasets:
    
    ########
    # Import and format the data table
    
    # Import formatted collapsed taxonomy file
    filename, directory, table = RottenIceModules.fileGet(
        'Select taxonomy-collapsed file for ' + d,
        file_type='csv', header_row=0, index_col=0
    )
    table = table[table.index.notnull()]
    
    # Get the list of colors
    colormap = table[colorcol]
    
    # Trim the table to just the data
    cols = [c for c in table.columns if c not in ignorecols]
    datatable = table[cols]
    
    # Trim further to just the samples of interest
    if samplelist != False:
        include_samples = filter_sample_ids(
            cols, samplelist, template = template)
        datatable = table[include_samples]
                
    # Parse sample names into components for labeling the plot
    samples_parsed = pd.DataFrame(
        index=datatable.columns, columns=['month', 'label', 'template'])
    for s in samples_parsed.index:
        m, f, r, g, t = RottenIceModules.samplenameToMonthFractionReplicateGeneTemplate(s)
        samples_parsed.loc[s, 'month'] = m
        samples_parsed.loc[s, 'label'] = f + '-' + r
        samples_parsed.loc[s, 'template'] = t
    
    # Figure out the month labels
    # Make sure 'month' is categorical with order preserved as it appears
    samples_parsed['month'] = pd.Categorical(
        samples_parsed['month'], categories=pd.unique(samples_parsed['month']),
        ordered=True)
    
    # If averaging is selected, average the data
    if average_replicates:
        datatable, samples_parsed = averageSamples(datatable, samples_parsed)
        
    # Keep original order of samples
    samples_parsed = samples_parsed.loc[datatable.columns]
    
    # Identify month groups in original order
    unique_months = samples_parsed['month'].cat.categories
    
    new_order = []
    month_separators = []
    y_pos = 0
    
    for month in unique_months:
        # samples for this month in original order
        month_samples = samples_parsed.index[samples_parsed['month'] == month]
        if new_order:
            y_pos += 1  # add a spacer row between months
        month_start = y_pos
        for col in month_samples:
            new_order.append(col)
            y_pos += 1
        month_end = y_pos - 1
        month_center = (month_start + month_end) / 2
        month_separators.append((month, month_center))
    
    # Reorder datatable preserving original order within months
    datatable = datatable[new_order]
    
    # Add spacer rows between months to match new bar positions (with the spacers)
    spacer = pd.Series(0, index=datatable.index)
    spaced_data = []
    final_order = []
    for i, col in enumerate(datatable.columns):
        spaced_data.append(datatable[col])
        final_order.append(col)
        if i < len(datatable.columns) - 1:
            curr_month = samples_parsed.loc[col, 'month']
            next_month = samples_parsed.loc[datatable.columns[i + 1], 'month']
            if curr_month != next_month:
                spacer_col = spacer.rename(f'space_{i}')
                spaced_data.append(spacer_col)
                final_order.append(spacer_col.name)
    
    datatable_spaced = pd.concat(spaced_data, axis=1)
    
    
    ########
    # Determine if template column should be shown
    unique_templates = samples_parsed['template'].unique()
    include_template_column = len(unique_templates) > 1
    
    # Adjust title to include template if only one
    plot_title = title + ' ' + d
    if not include_template_column:
        plot_title += f' ({unique_templates[0]})'
    
    # Determine height scaling
    num_bars = len([col for col in datatable_spaced.columns if not col.startswith('space_')])
    bar_height = 0.2 if num_bars >= 20 else 0.5
    
    # Build the barplot
    fig, ax = plt.subplots(figsize=(24, num_bars * bar_height))
    datatable_spaced.T.plot(
        kind='barh',
        stacked=True,
        color=colormap.tolist(),
        ax=ax,
        width=0.9
    )
    
    # Flip y-axis to have first sample on top
    ax.invert_yaxis()
    
    # Manually position the bars
    bar_positions = np.arange(len(datatable_spaced.columns))
    ax.set_yticks(bar_positions)
    
    # Hide default y tick labels
    ax.set_yticklabels([])
    
    # Get axis limits for positioning
    xmin, xmax = ax.get_xlim()
    
    # Set base x-position for the label(s) left of the bar plot
    x_label_pos = -0.08 * xmax
    label_x = x_label_pos
    template_x = x_label_pos + (0.045 * xmax if include_template_column else 0.0)
    
    # Add label(s) left of each bar
    for y_pos, col in enumerate(datatable_spaced.columns):
        if col.startswith('space_'):
            continue
        # Sample label (column 1)
        label = samples_parsed.loc[col, 'label']
        ax.text(
            label_x, y_pos, label, ha='left', va='center',
            fontsize=fontsize, fontfamily='monospace'
        )
        # Template (column 2), if applicable
        if include_template_column:
            template = samples_parsed.loc[col, 'template']
            ax.text(
                template_x, y_pos, template, ha='left', va='center',
                fontsize=fontsize, fontfamily='monospace'
            )
    
    # Remove spines and y-axis ticks
    for spine in ['left', 'bottom', 'top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.set_yticklabels([''] * len(bar_positions))
    ax.tick_params(axis='y', which='both', length=0)
    
    # Add month labels (vertical, aligned to left of label column)
    for month, center in month_separators:
        ax.text(
            x=x_label_pos - 0.02, y=center,
            s=monthmap[month], ha='center', va='center', rotation='vertical',
            fontsize=fontsize, weight='bold', transform=ax.transData
        )
    
    # Adjust x-limits to fit labels and bars
    ax.set_xlim(x_label_pos - 0.03 * xmax, xmax)
    
    # Axis labels and title
    ax.set_title(plot_title, fontweight='bold', fontsize=fontsize + 2)
    ax.set_xlabel('Relative Abundance')
    ax.set_ylabel('Samples')
    
    # Add legend (taxa names)
    ax.legend(
        loc='center left',
        bbox_to_anchor=(0.98, 0.9),
        title='Taxa',
        fontsize=fontsize,
        frameon=False
    )
    
    plt.tight_layout()
    plt.show()
    
    # Save figure
    fig.savefig(filename[:-4] + '.pdf', transparent=True)
    fig.savefig(filename[:-4] + '.svg', transparent=True)
    
