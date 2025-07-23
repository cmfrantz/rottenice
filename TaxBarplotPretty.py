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
        columns = taxonomy, color, ranking, max, avg, then a list of sample names
        Column 'color' defines the colors to use for each taxonomy value.
        Columns ranking, max, and avg are ignored
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

datasets = ['16S noUEMC QC', '18S noUBA QC'] # List of the datasets to generate plots for
colorcol = 'color'  # The column header in the data table that defines the color hex value
ignorecols = ['supergroup','Big','Primary Producer','color',
              'ranking','max','avg']     # Non-data columns in the data table
fontsize = 11 # Font size for the plot
title = 'Relative abundance of major taxa' # Title header for the plot
monthmap = {
    'M'     : 'May',
    'JN'    : 'June',
    'JY10'  : 'July 10 dirty floe',
    'JY11'  : 'July 11 clean floe'
    }


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
    # Build the barplot
    fig, ax = plt.subplots(figsize=(24, len(datatable_spaced.columns) * 0.2))  # less height per bar
    datatable_spaced.T.plot(
        kind='barh',
        stacked=True,
        color=colormap.tolist(),
        ax=ax,
        width=0.9  # To decrease spacing between bars, increase this number
    )
    
    # Flip y-axis to have first sample on top
    ax.invert_yaxis()
    
    # Manually position the bars
    bar_positions = np.arange(len(datatable_spaced.columns))
    ax.set_yticks(bar_positions)
    
    # Hide default y tick labels:
    ax.set_yticklabels([])
    
    # Get axis limits for positioning
    xmin, xmax = ax.get_xlim()
    
    # Set base x-position for the two-column labels left of the bar plot
    # Column 1: the sample type label (material+horizon-replicate)
    x_label_pos = -0.08 * xmax  # changes the spacing between the labels and
                                # the plot edge
                                # decreasing the value decreases the spacing
                                # increasing it increases the spacing (and overlap)
    label_x = x_label_pos
    # Column 2: the template (cDNA or DNA)
    template_x = x_label_pos + 0.045 * xmax # changes the spacing between
                                # the two label "columns"
                                # increasing the value increases the spacing
                                # decreasing the value tightens it
    
    # Add the two-column labels left of bar plot
    for y_pos, col in enumerate(datatable_spaced.columns):
        if col.startswith('space_'):
            continue
        # Sample label (column 1)
        label = samples_parsed.loc[col, 'label']
        ax.text(
            label_x, y_pos, label, ha='left', va='center',
            fontsize=fontsize, fontfamily='monospace')
        # Sample template (column 2)
        template = samples_parsed.loc[col, 'template']
        ax.text(
            template_x, y_pos, template, ha='left', va='center',
            fontsize=fontsize, fontfamily='monospace')
    
    # Remove left spine and ticks (no lines or tick marks)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticklabels([''] * len(bar_positions))
    ax.tick_params(axis='y', which='both', length=0)
    
    # Add month labels centered next to each month's group of samples
    for month, center in month_separators:
        ax.text(
            x=x_label_pos-0.02, y=center,
            s=monthmap[month], ha='center', va='center', rotation = 'vertical',
            fontsize=fontsize, weight='bold', transform=ax.transData,
        )
    
    # Adjust x limits to fit labels and bars
    # Decrease the number to shift left, increase to shift right
    ax.set_xlim(x_label_pos - 0.03 * xmax, xmax)
    
    # Axis labels and title
    ax.set_title(title + ' ' + d, fontweight='bold', fontsize = fontsize+2)
    # Title is the common title plus the name of the dataset
    ax.set_xlabel('Relative Abundance')
    ax.set_ylabel('Samples')
    
    # Add the legend (taxa names)
    ax.legend(
        loc='center left',
        bbox_to_anchor=(0.98, 0.9), # horizontal from left, vertical from bottom
        title='Taxa',
        fontsize=fontsize,
        frameon=False
    )
    
    plt.tight_layout()
    plt.show()
    
    # Save the figure
    fig.savefig(filename[:-4] + '.pdf', transparent=True)
    fig.savefig(filename[:-4] + '.svg', transparent=True)
