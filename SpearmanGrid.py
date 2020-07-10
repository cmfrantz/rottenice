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

This script was created as part of the Rotten Ice Project

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
    pip install progress
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


Copyright (C) 2020  Carie M. Frantz

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
from progress.bar import Bar

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
genes = RottenIceVars.genes
templates = RottenIceVars.templates

# Samples to exclude
otherlist = RottenIceVars.other_samples

# Variables to plot
varlist = ['day_num', 'horizon_code', 'sample_mat_code',
           'coord_lat', 'coord_lon', 'temperature', 'pH', 'bulk_density',
           'salinity_insitu', 'salinity', 'SPM', 'DOC', 'Chl', 'Phaeo',
           'FoFa', 'PAM', 'SedLoad', 'bact_cell_ct', 'CTC',
           'bact_active_cell_ct', 'diatom_ct', 'phyto_ct_all',
           'phyto_ct_other', 'pEPS', 'POC', 'nitrogen', 'CN']

# Colormap to use - pick a diverging colormap
cmap = RottenIceModules.genDivergingCmap()

# Significance cutoff for plot:
# Correlation coefficients with p values greater than this are excluded from
# the plot so that only "significant" correlations are shown
p_cutoff = 0.01
# Number of taxonomic groups to consider in each plot
n_groups = 20

# Info about file naming
file_info = RottenIceVars.file_sets['spearman']
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
 amplicon recovery). Analysis done by C. Frantz, June 2020.
'''

####################
# FUNCTIONS
####################

def buildCorrTable(metadata, abundance_table, samples):
    '''Build the metadata vs. taxonomic abundance table
        containing calculated Spearman correlation coefficients'''
    corrtable = pd.DataFrame(data = None, index = abundance_table.index,
                             columns = varlist)
    # Remove replicate number from sample names
    msamples = [sample[0:-2] for sample in samples]
    
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
            if p <= 0.05:
                corrtable.loc[group, var] = corr
                
    return corrtable.astype(float)
    

def buildHeatmap(ax, corrtable, dset):
    '''Build the heatmap of correlation coefficients'''
    # Build the figure
    im = ax.imshow(corrtable, cmap = cmap, vmin = -1, vmax = 1)
    
    # Place x and y ticks and labels
    ax.set_xticks(np.arange(corrtable.shape[1]))
    ax.set_yticks(np.arange(corrtable.shape[0]))
    ax.set_xticklabels(corrtable.columns, ha = 'right',
                       rotation = 45, rotation_mode = 'anchor')
    ax.set_yticklabels(corrtable.index)
    ax.tick_params(labelsize = 9)
    
    # Add grid to plot
    ax.set_xticks(np.arange(corrtable.shape[1]+1)-0.5, minor = True)
    ax.set_yticks(np.arange(corrtable.shape[0]+1)-0.5, minor = True)
    ax.grid(which = 'minor', color = 'k', linestyle = '-', linewidth = 1)
    ax.tick_params(which = 'minor', bottom = False, left = False)
    
    # Add title to the plot
    ax.set_title(dset + ' L' + str(level))
    
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax = ax)
    cbar.ax.set_ylabel("Spearman Correlation", va = 'bottom')
    
    return im
        
    
# Annotate heatmap
def annotateHeatmap(im, data=None, valfmt = "{x:.2f}", **textkw):
    '''Annotate the heatmap with text displaying correlation coefficients'''
    if not isinstance(data, (list, np.ndarray)):                                # What does t his do?
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
    # Loop through each dataset
    for dset in asv_tables:
        print('*******************************************************\n'
              'Determining Spearman correlations for ' + dset + ' dataset\n'
              '*******************************************************')
        gene = dset.split('-')[0]
        max_level = genes[gene]['max_level']
        
        # Format data
        data, samples = RottenIceModules.formatOTUtableData(
            asv_tables[dset], max_level = max_level,
            tax_reassign_list = genes[gene]['tax_reassign_list'])
        # Get rid of any samples with no data and update the sample list ############
        data = data.dropna(axis = 1)
        samples = data.columns
        non_sample_cols = ['taxonomy'] + RottenIceModules.levelCols(max_level)
        samples = [sample for sample in samples
                   if sample not in non_sample_cols]
        samples = [sample for sample in samples
                   if not any(other in sample for other in otherlist)]
        data = data[samples + non_sample_cols]
        
        # Set up figure grid
        print('Preparing figure grid...')
        fig, axes = plt.subplots(max_level, 1,
                                 figsize = (20, max_level * n_groups * 0.5))
        
        # Perform correlation calculations
        # and build heatmaps for each taxonomic level
        bar = Bar('Calculating at each taxonomic level...',
                  max = max_level) # Progress bar
        # Loop through each taxonomic level
        for level in np.arange(1,max_level+1,1):
            # Condense the dataset to the level
            ds = RottenIceModules.condenseDataset(data, level, samples)
            # Normalize the dataset
            ds = ds / ds.sum(axis = 0)
            # Find the most abundant ASVs
            ds['sums'] = ds.sum(axis = 1)
            ds = ds.sort_values(by = ['sums'], ascending = False)
            if ds.shape[0] >= n_groups:
                ds = ds.iloc[0:n_groups-1]
            # Calculate Spearman correlation and return significant values
            corrtable = buildCorrTable(metadata, ds, samples)
            # Build heatmap
            ax = axes[level-1]
            im = buildHeatmap(ax, corrtable, dset)
            texts = annotateHeatmap(im, fontsize = 6)
            bar.next()
        bar.finish()
        
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