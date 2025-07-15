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
    Metadata table (csv)
        where rows = samples, columns = metadata characteristics
        header row and index column are specified in the variables below
        
    ASV tables (csv) for every sequence dataset analyzed
        where rows = ASVs, columns = samples
        header row and index column are specified in the variables below

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
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    
    # Import Metadata table
    filename, directory, metadata, = RottenIceModules.fileGet(
        'Select metadata table', tabletype = 'metadata')
    metadata = metadata.dropna(how = 'all')
    metadata = metadata.replace('na', np.nan)
    # Build matrix of metadata values
    metadata = metadata[varlist]
    
    # Import ASV tables
    asv_tables = {}
    for gene in genes:
        for template in templates:
            dset = gene + '-' + template
            filename, directory, data = RottenIceModules.fileGet(
                'Select ' + dset + ' ASV table', tabletype = 'OTU-table',
                directory = directory)
            asv_tables[dset] = data
            
    # Set up file navigation html
    filenames = [(dset, out_filename + '_' + dset) 
                 for dset in asv_tables]
    nav_html = RottenIceVars.nav_html_start
    for dset in filenames:
        nav_html = (nav_html + ' <a href = "' + dset[1] + '.html">'
                    + dset[0] + '</a>  /')
    nav_html = nav_html[:-3]
    
    # Build heatmaps
    spearman_tables = {}
    # Loop through each dataset
    for dset in asv_tables:
        print('*******************************************************\n'
              'Determining Spearman correlations for ' + dset + ' dataset\n'
              '*******************************************************')
              
        gene = dset.split('-')[0]
        
        # Format & condense data by level
        ASV_table, samples = RottenIceModules.formatASVTable(
            asv_tables[dset], tax_levels, max_level)
        level_tables = RottenIceModules.condenseASVTable_by_TaxLevel(
            ASV_table, tax_levels[0:max_level], samples)

        # Set up figure grid
        print('Preparing figure grid...')
        fig, axes = plt.subplots(max_level, 1,
                                 figsize = (20, max_level * n_groups * 0.5))
        
        # Perform correlation calculations
        # and build heatmaps for each taxonomic level
        # Condense datasets
        # Loop through each taxonomic level
        for l in np.arange(0,max_level,1):
            level = tax_levels[l]            
            print('Taxonomic level ' + str(l+1) + ' (' + level + ')...')
            ds = level_tables[level].copy()
            # normalize the dataset
            ds = ds/ds.sum(axis = 0)
            # find the most abundant taxa
            ds['sums'] = ds.sum(axis=1)
            ds = ds.sort_values(by=['sums'], ascending=False)
            if ds.shape[0] >= n_groups:
                ds = ds.iloc[0:n_groups]
                
            # Calculate Spearman correlation and return significant values
            print('  Calculating Spearman correlation coefficients for taxa')
            corrtable = buildCorrTable(
                metadata, ds, varlist, samples, significance = p_cutoff)
            
            # Save correlation table
            spearman_tables[dset+'_L'+str(l+1)] = corrtable
        
            # Build heatmap table to display results of Spearman correlations
            ax = axes[l]
            im = buildHeatmap(
                ax, cmap, corrtable,
                ('Spearman correlation of ' + dset +
                 ' data at taxonomic level ' + str(l+1) + ' ' + level))
            texts = annotateHeatmap(im, fontsize = 6)
            
        # Save plots
        print('Saving ' + dset + ' figure files...')
        fig.tight_layout()
        filename = (out_filename + '_' + dset)
        img_filename = filename + '.svg'
        fig.savefig(directory + '\\' + img_filename,
                    transparent = True)

        # Generate HTML
        RottenIceModules.genHTMLfile(
            directory + '\\' + filename + '.html',
            file_info['title'], subtitle_text,
            img_filename,
            page_nav_html = nav_html,
            alt_text = ('Heatmaps showing Spearman correlation coefficients '
                        + 'for metadata and abundant taxa'))
