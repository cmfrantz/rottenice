# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 13:13:03 2025

@author: cariefrantz
"""

####################
# IMPORTS
####################
import os
import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'Arial'

import RottenIceModules

####################
# VARIABLES
####################
outdir = os.getcwd()

otu_datasets = ['16S noUEMC QC', '18S noUBA QC']
algae_dataset = ['AlgaeID']

# Spreadsheet variables
colorcol = 'color'  # The column header in the data table that defines the color hex value
ignorecols = ['Phylum','phylum','Supergroup','supergroup',
              'Big','PrimProd','Primary Producer','SylvieID',
              'color','ranking','rank','max','avg','total','Total']     # Non-data columns that could be in the data table

# Plot details
fontsize = 11 # Font size for the plot
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


def normalizeTable(table):
    '''Normalizes the table to produce relative abundance from absolute counts'''
    col_sums = table.sum(axis=0)
    normalized = table.divide(col_sums.where(col_sums != 0, np.nan), axis=1)
    normalized = normalized.fillna(0)
    return normalized


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

    samples_parsed = samples_parsed.copy()

    # Remove replicate from label
    samples_parsed['samplename'] = samples_parsed['samplename'].astype(str).str.split('-').str[:-1].str.join('-')

    # Combine parts into a unique sample base name
    samples_parsed['sample_base'] = (
        samples_parsed['month'].astype(str) + '-' +
        samples_parsed['samplename'] + '.' +
        samples_parsed['template'].astype(str)
    )

    # Preserve order of first occurrence
    sample_base_order = samples_parsed.drop_duplicates('sample_base')['sample_base']
    samples_parsed['sample_base'] = pd.Categorical(samples_parsed['sample_base'], categories=sample_base_order, ordered=True)

    # Transpose data so samples are rows
    datatable_T = datatable.T.copy()
    datatable_T['sample_base'] = samples_parsed['sample_base']

    # Average replicates in original order
    averaged = datatable_T.groupby('sample_base', sort=False).mean()

    # Transpose back
    datatable_avg = averaged.T

    # Rebuild metadata table in matching order
    samples_parsed_avg = (
        samples_parsed
        .drop_duplicates('sample_base')
        .set_index('sample_base')
        .loc[datatable_avg.columns]
    )

    return datatable_avg, samples_parsed_avg



def parse_sample_names(sample_ids, style = 'short'):
    """
    Parses sample names into a DataFrame with columns:
    ['month', 'label', 'template']
    
    style :
        'short' : labels are sample codes
        'horizon' : labels are long-format horizon names
    """
    parsed = []
    for sid in sample_ids:
        # If the sample name format is the full sample ID with template,
        # e.g., 'M-CS-HT-1.16S-cDNA'
        if '.' in sid:
            month, fraction, replicate, gene, template = RottenIceModules.samplenameToMonthFractionReplicateGeneTemplate(
                sid)
            samplename = '-'.join([month,fraction,replicate])
        # Otherwise, the sample name format is assumed to be short-format,
        # e.g., 'M-CS-HT-1' or 'M-CS-HT'
        else:
            samplename = sid
            template = ''
            splits = sid.split('-')
            month = splits[0]
            if splits[-1].isdigit():
                replicate = splits[-1]
                fraction = splits[-2]
            else:
                replicate = ''
                fraction = splits[-1]
            
        # Determine the label to use
        if style == 'short':
            if replicate == '':
                label = fraction
            else:
                label = '-'.join([fraction,replicate])
        if style == 'horizon':
            f_name, h_name = RottenIceModules.parseFraction(fraction)
            label = h_name
        parsed.append((sid, samplename, month, label, template))
        
    df = pd.DataFrame(parsed,
                      columns=['sample', 'samplename', 'month', 'label', 'template']
                      ).set_index('sample')
    return df



def plot_taxa_barplot(datatable, samples_parsed, colormap, title, monthmap,
                      xlabel = 'Relative Abundance',
                      output_file=None, fontsize=11,
                      label_columns=['label', 'template']):
    """
    Generate a stacked horizontal barplot for relative abundance data.

    Parameters
    ----------
    datatable : pd.DataFrame
        Taxa as rows, sample names as columns (in desired order).
    samples_parsed : pd.DataFrame
        Metadata for samples. Index = column names of datatable. Columns include 'month'.
    colormap : pd.Series
        Color values corresponding to taxa (index matches datatable).
    title : str
        Title for the plot.
    monthmap : dict
        Mapping from month codes to month labels.
    output_file : str or None
        If given, saves to PDF and SVG at this path (excluding extension).
    fontsize : int
        Font size for plot text.
    label_columns : list of str
        Which columns in samples_parsed to use for sample labels.
        If empty, no labels will be drawn. If a two-item list, the two items
        (e.g., 'label' and 'template') will appear in seperate columns.
    """

    samples_parsed['month'] = pd.Categorical(
        samples_parsed['month'], categories=pd.unique(samples_parsed['month']),
        ordered=True)

    unique_months = samples_parsed['month'].cat.categories
    new_order = []
    month_separators = []
    y_pos = 0
    for month in unique_months:
        month_samples = samples_parsed.index[samples_parsed['month'] == month]
        if new_order:
            y_pos += 1
        month_start = y_pos
        for col in month_samples:
            new_order.append(col)
            y_pos += 1
        month_end = y_pos - 1
        center = (month_start + month_end) / 2
        month_separators.append((month, center))

    datatable = datatable[new_order]
    spacer = pd.Series(0, index=datatable.index)
    spaced_data = []
    for i, col in enumerate(datatable.columns):
        spaced_data.append(datatable[col])
        if i < len(datatable.columns) - 1:
            curr_month = samples_parsed.loc[col, 'month']
            next_month = samples_parsed.loc[datatable.columns[i + 1], 'month']
            if curr_month != next_month:
                spacer_col = spacer.rename(f'space_{i}')
                spaced_data.append(spacer_col)
    datatable_spaced = pd.concat(spaced_data, axis=1)

    num_bars = len([c for c in datatable_spaced.columns if not c.startswith('space_')])
    
    # Dynamic bar height based on number of bars
    bar_height = min(0.5, max(0.2, 10 / num_bars))

    fig, ax = plt.subplots(figsize=(24, num_bars * bar_height))
    datatable_spaced.T.plot(kind='barh', stacked=True, color=colormap.tolist(), ax=ax, width=0.9)
    ax.invert_yaxis()
    bar_positions = np.arange(len(datatable_spaced.columns))
    ax.set_yticks(bar_positions)
    ax.set_yticklabels([])

    xmin, xmax = ax.get_xlim()
    x_offsets = [-0.08 * xmax + 0.045 * xmax * i for i in range(len(label_columns))]

    for y_pos, col in enumerate(datatable_spaced.columns):
        if col.startswith('space_'):
            continue
        for i, label_col in enumerate(label_columns):
            label = samples_parsed.loc[col, label_col]
            ax.text(x_offsets[i], y_pos, label, ha='left', va='center',
                    fontsize=fontsize, fontfamily='monospace')

    for spine in ['left', 'bottom', 'top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)

    for month, center in month_separators:
        month_label_x = x_offsets[0] - 0.04 * (xmax - xmin)
        ax.text(month_label_x, center, monthmap[month], ha='center', va='center',
                rotation='vertical', fontsize=fontsize, weight='bold', transform=ax.transData)

    ax.set_xlim(x_offsets[0] - 0.03 * xmax, xmax)

    ax.set_title(title, fontweight='bold', fontsize=fontsize + 2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Samples')

    ax.legend(loc='center left', bbox_to_anchor=(0.98, 0.9),
              title='Taxa', fontsize=fontsize, frameon=False)

    plt.tight_layout()
    plt.show()

    if output_file:
        fig.savefig(output_file + '.pdf', transparent=True)
        fig.savefig(output_file + '.svg', transparent=True)



#%%
####################
# MAIN FUNCTION
####################
    
#####################
# Plot Set 1: Full dataset plots

for d in otu_datasets:

    template = 'both'
    
    # Import formatted collapsed taxonomy file
    filename, directory, table = RottenIceModules.fileGet(
        'Select taxonomy-collapsed Top 20 taxa file for ' + d,
        file_type='csv', header_row=0, index_col=0, directory = outdir
    )
    table = table[table.index.notnull()]
    
    # Get the list of colors
    colormap = table[colorcol]
    
    # Trim the table to just the data
    cols = [c for c in table.columns if c not in ignorecols]
    datatable = table[cols]
    
    # Parse the sample names into their metadata to add to the plot labels
    samples_parsed = parse_sample_names(cols)
    
    # Build the plot
    plot_taxa_barplot(
        datatable, samples_parsed, colormap,
        'Relative abundance of major taxa', monthmap,
        output_file = outdir + '\\otu_barplot_all ' + d, 
        fontsize=fontsize, label_columns=['label', 'template'])
   
    
#####################
# Plot Set 2: Whole-horizon melt averages  

for d in otu_datasets:    
    # Import formatted collapsed taxonomy file
    filename, directory, table = RottenIceModules.fileGet(
        'Select taxonomy-collapsed Top 10 taxa file for ' + d,
        file_type='csv', header_row=0, index_col=0, directory = outdir
    )
    table = table[table.index.notnull()]
    
    # Get the list of colors
    colormap = table[colorcol]
    
    # Trim the table to just the data
    cols = [c for c in table.columns if c not in ignorecols]
    datatable = table[cols]
        
        
    for template in ['DNA','cDNA']:
        
        # Trim further to just the samples of interest
        samplelist = RottenIceModules.genSampleList(['M','JN','JY10','JY11'],['HT','HM','HB'])
        include_samples = filter_sample_ids(
                cols, samplelist, template = template)
        datatable = table[include_samples]
        
        # Parse the sample names into their metadata to add to the plot labels
        samples_parsed = parse_sample_names(include_samples, style = 'horizon')
        
        # Average the samples
        datatable, samples_parsed = averageSamples(datatable, samples_parsed)
        
        # Build the plot
        plot_taxa_barplot(
            datatable, samples_parsed, colormap,
            'Average relative abundance of major taxa (' + d + ' ' + template + ')',
            monthmap,
            output_file = outdir + '\\otu_barplot_H_avg ' + d + ' ' + template, 
            fontsize=fontsize, label_columns=['label'])
    
    
#####################
# Plot Set 3: Algae IDs 
# Get Algae ID file
filename, directory, table = RottenIceModules.fileGet(
    'Select Algae ID taxonomy table',
    file_type='csv', header_row=0, index_col=0, directory = outdir
)

table = table[table.index.notnull()]

# Get the list of colors
colormap = table[colorcol]

# Trim the table to just the data
cols = [c for c in table.columns if c not in ignorecols]

# Get a sorted sample list
possibleSamples = RottenIceModules.genSampleListCS()
samplelist = [s for s in possibleSamples if s in cols]

datatable = table[samplelist]
normtable = normalizeTable(datatable)

# Parse labels
samples_parsed = parse_sample_names(samplelist)

# Build the plots
plot_taxa_barplot(
    datatable, samples_parsed, colormap,
    'Absolute abundance of algal groups', monthmap,
    output_file = outdir + '\\algae_barplot_abs', fontsize = fontsize,
    xlabel = 'Abundance ($10^{6}\\,\\mathrm{cells} \\cdot L^{-1}$)',
    label_columns=['label'])

plot_taxa_barplot(
    normtable, samples_parsed, colormap,
    'Relative abundance of algal groups', monthmap,
    output_file = outdir + '\\algae_barplot_rel', fontsize = fontsize,
    label_columns=['label'])