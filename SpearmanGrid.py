#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""Rotten Ice Project: SpearmanGrid

Created on Sun May 31 15:24:03 2020
@author: cariefrantz
@project: RottenIce

CALCULATES SPEARMAN CORRELATION FOR METADATA VS. TAXONOMIC GROUPS
This script:
    - Calculates the Spearman correlation coefficients relating selected
        metadata parameters to the abundance of sequences retrieved for the
        20 most abundant taxonomic groups at different taxonomic levels.
    - Generates heatmap grids for significant correlations between each
        metadata variable / taxonomic group pair.
    - Produces HTML files for each sequence template type with all of the
        heatmaps (at each taxonomic level).

This script was created as part of the Rotten Ice Project.
It was substantially overhauled in 2025 to use simpler files.

Arguments:  None

Requirements:   
    Metadata table (tsv)
        where rows = sample IDs, columns = metadata characteristics
        This is the same metadata table used in QIIME2 which contains
        metadata for every sample - gene - template combination
        
    ASV tables (csv) for every sequence dataset analyzed
        where rows = taxonomy calls, columns = samples,
        values = taxonomy counts for each sample
        These files are produced by exporting Level 7 taxonomic data from
        QIIME2 barplots in Qiime2View, transposing them,
        and seperating them by template (cDNA vs. DNA)

Example in command line:
    python SpearmanGrid.py

Dependencies Install:
    sudo apt-get install python3-pip python3-dev
    pip install tkinter
    pip install numpy
    pip install pandas
    pip install matplotlib
    pip install math
    pip install bokeh

You will also need to have the following files
    in the same directory as this script.
They contain modules and variables that this script calls.
    RottenIceModules.py
    RottenIceVars.py
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
import pandas as pd

import numpy as np
from scipy.stats import spearmanr

import matplotlib
from matplotlib import pyplot as plt
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import RottenIceModules
import RottenIceVars


####################
# VARIABLES
####################

# Types of sequence data
genes = ['16S','18S','Eukaryotic Primary Producers (18S)']
max_level = 7
templates = RottenIceVars.templates

# Samples to exclude
otherlist = RottenIceVars.other_samples

# Variables to plot
varlist = ['month_num','horizon_num','replicate','lat','lon','depth_in_ice',
           'pH','salinity','SPM','DOC','chl','phaeo','FoFa','PAM','SedLoad',
           'cell_ct','CTC','cells_act','pEPS','POC','nitrogen','CN',
           'temperature','bulk_ice_density','diatom_ct','phyto_ct_select',
           'phyto_ct_all','phyto_ct_other','gels_total','gels_sm','gels_md',
           'gels_lg','gels_xl'
]

# Define taxonomic levels
tax_levels = RottenIceVars.tax_levels

# Colormap to use - pick a diverging colormap
cmap = RottenIceModules.genDivergingCmap()

# Significance cutoff for plot:
# Correlation coefficients with p values greater than this are excluded from
# the plot so that only "significant" correlations are shown
p_cutoff = 0.01
# Number of taxonomic groups to consider in each plot
n_groups = 20

# Info about file naming
file_info = RottenIceVars.file_sets['spearman_taxonomy']
out_filename = file_info['pfx']
title = file_info['title']

# HTML code for the for the HTML files
# Subheading text (Part II after the level number)
subtitle_text = '''
Correlation between biogeochemical variables and relative abundance of major
taxonomic groups. Values shown are Spearman correlation coefficients
calculated using the
<a href="https://docs.scipy.org/doc/scipy/reference/generated/
scipy.stats.spearmanr.html">scipy.stats.spearmanr package</a>
(1.4.1) for python. 
Only significant correlations (p <=''' + str(p_cutoff) + ''') are shown. 
Positive numbers indicate a positive correlation, 
negative numbers indicate a negative correlation; 
larger numbers indicate stronger correlations. 
The top ''' + str(n_groups) + ''' most abundant taxonomic groups for each 
level were analyzed (across all samples, with abundance normalized to total 
 amplicon recovery). Created using the script 
<a href="https://github.com/cmfrantz/rottenice">
SpearmanGrid.py</a>. Analysis done by C. Frantz, June 2025.
'''


####################
# FUNCTIONS
####################

def buildCorrTable(metadata, abundance_table, varlist, samples,
                   significance = 0.05):
    '''
    Calculates Spearman correlation coefficients for metadata vs. relative
    taxonomic abundance

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata table, where columns are metadata parameters and index is
        the sample list.
    abundance_table : pd.DataFrame
        Taxonomic abundance table at some taxonomic depth, where columns are
        samples and rows are the taxa, values are the relative abundance for
        each taxon.
    varlist : list of str
        List of metadata parameters (must be numerical) to correlate.
    samples : list of str
        List of samples to analyze.
    significance : float
        Cutoff value for assessing significance. Default = 0.05.

    Returns
    -------
    corrtable : pd.DataFrame
        Table of Spearman correlation values, if values meet the significance
        threshhold.

    '''
    corrtable = pd.DataFrame(data = None, index = abundance_table.index,
                             columns = varlist)
    # Remove replicate number from sample names
    msamples = samples.copy()
    # msamples = [sample[0:-2] for sample in samples]
    
    # Loop through each taxonomic group
    for group in abundance_table.index:
        # Loop through each metadata variable
        for var in varlist:
            # Grab the list of values for each sample
            varvals = metadata.loc[msamples, var].to_numpy(dtype=float)
            counts = abundance_table.loc[group, samples].to_numpy(dtype=float)
            
            #Calculate Spearman coefficient
            corr, p = spearmanr(varvals, counts, nan_policy='omit')
            
            # If p value <= 0.05, assign the correlation a color on the
            # diverging color scale and add it to the heatmap matrix
            if p <= significance:
                corrtable.loc[group, var] = corr
                
    return corrtable.astype(float)
    

def buildHeatmap(ax, cmap, corrtable, title, vmin = -1, vmax = 1):
    '''
    Builds heatmap figure from a table of correlation values

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Handle for the subplot Axes object.
    cmap : matplotlib.colors.ListedColormap
        Colormap object for mapping colors to values
    corrtable : pandas.DataFrame
        DataFrame containing a table of correlation coefficients between the
        index and columns.
    title : str
        Text used to label the heatmap
    vmin: int or float
        Value setting the minimum possible value for the correlation coeffs.
        Default is -1.
    vmax: int or float
        Value setting the maximum posible value for the correlation coeffs.
        Default is 1.

    Returns
    -------
    im : matplotlib.axes.Axes.imshow
        Displays image, returns handle.

    '''
    # Build the figure
    im = ax.imshow(corrtable, cmap, vmin = -1, vmax = 1)
    
    # Place x & y ticks and labels
    ax.set_xticks(np.arange(corrtable.shape[1]))
    ax.set_yticks(np.arange(corrtable.shape[0]))
    ax.set_xticklabels(
        corrtable.columns, ha = 'right', rotation = 45,
        rotation_mode = 'anchor')
    ax.set_yticklabels(corrtable.index)
    ax.tick_params(labelsize = 9)
    
    # Add grid to plot
    ax.set_xticks(np.arange(corrtable.shape[1]+1)-0.5, minor=True)
    ax.set_yticks(np.arange(corrtable.shape[0]+1)-0.5, minor=True)
    ax.grid(which = 'minor', color = 'k', linestyle = '-', linewidth=1)
    ax.tick_params(which = 'minor', bottom = False, left = False)
    
    # Add title to plot
    ax.set_title(title)
    
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax = ax)
    cbar.ax.set_ylabel('Spearman Correlation', va = 'bottom')
    
    return im
        
    
# Annotate heatmap
def annotateHeatmap(im, data=None, valfmt = "{x:.2f}", **textkw):
    '''Annotate the heatmap with text displaying correlation coefficients'''
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()
        
    # Set default alignment to center, but allow it to be 
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create text for each plotted value.
    # Change the text color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color='w')
            val = data[i, j]
            if np.ma.is_masked(val) or np.isnan(val):
                continue
            text = im.axes.text(j, i, valfmt(val, None), **kw)
            texts.append(text)

    return texts
#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':

    # Import metadata
    filename, directory, metadata = RottenIceModules.fileGet(
        'Select metadata table', tabletype='metadata-qiime')
    metadata = metadata.dropna(how='all').replace('na', np.nan)
    metadata = metadata[varlist]

    # Import taxonomy tables
    tax_tables = {}
    for gene in genes:
        for template in templates:
            dset = gene + '-' + template
            filename, directory, data = RottenIceModules.fileGet(
                'Select ' + dset + ' ASV table', tabletype='OTU-table', directory=directory)
            tax_tables[dset] = data

    # Set up file navigation HTML
    filenames = [(dset, out_filename + '_' + dset) for dset in tax_tables]
    nav_html = RottenIceVars.nav_html_start
    for dset in filenames:
        nav_html += f' <a href="{dset[1]}.html">{dset[0]}</a> /'
    nav_html = nav_html[:-3]

    # Correlation heatmaps per dataset
    spearman_tables = {}

    for dset in tax_tables:
        print(f'\n******** Processing {dset} ********')
        gene = dset.split('-')[0]
        img_filenames = []

        # Format and condense
        ASV_table, samples = RottenIceModules.formatASVTable(
            tax_tables[dset], tax_levels, max_level)
        level_tables = RottenIceModules.condenseASVTable_by_TaxLevel(
            ASV_table, tax_levels[0:max_level], samples)

        # Loop through each taxonomic level
        for l in range(max_level):
            level = tax_levels[l]
            print(f'Taxonomic level {l+1} ({level})...')
            ds = level_tables[level].copy()

            # Normalize and select top groups
            ds = ds / ds.sum(axis=0)
            ds['sums'] = ds.sum(axis=1)
            ds = ds.sort_values(by='sums', ascending=False)
            ds = ds.iloc[:n_groups] if ds.shape[0] >= n_groups else ds

            # Spearman correlation
            print('  Calculating Spearman correlation coefficients...')
            corrtable = buildCorrTable(metadata, ds, varlist, samples, significance=p_cutoff)
            spearman_tables[f'{dset}_L{l+1}'] = corrtable

            # Build heatmap
            fig, ax = plt.subplots(figsize=(20, n_groups * 0.5))
            im = buildHeatmap(ax, cmap, corrtable,
                              f'Spearman correlation of {dset} at taxonomic level {l+1} {level}')
            annotateHeatmap(im, fontsize=6)
            fig.tight_layout()

            # Save figure
            img_filename = f"{out_filename}_{dset}_L{l+1}.svg"
            full_img_path = f"{directory}\\{img_filename}"
            fig.savefig(full_img_path, transparent=True)
            img_filenames.append(img_filename)
            plt.close(fig)

        # Generate HTML per dataset
        html_path = f"{directory}\\{out_filename}_{dset}.html"
        html_title = file_info['title']
        html_body = subtitle_text

        for img in img_filenames:
            html_body += f'<br><img src="{img}" width="1000" alt="Heatmap: {img}"><br>'

        RottenIceModules.genHTMLfile(
            html_path,
            html_title,
            subtitle_text,
            image_filepaths=img_filenames,
            alt_text=["Heatmap: " + img for img in img_filenames],
            page_nav_html=nav_html
        )

        print(f'Saved HTML and images for {dset} to {directory}')
