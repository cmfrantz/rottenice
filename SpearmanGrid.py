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
    pip install os
    pip install tkinter
    pip install progress
    pip install pandas
    pip install numpy
    pip install scipy
    pip install matplotlib


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
import os
from tkinter import filedialog
from tkinter import *
from progress.bar import Bar

import pandas as pd

import numpy as np
from scipy.stats import spearmanr

import matplotlib
from matplotlib import pyplot as plt
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


####################
# VARIABLES
####################

# Input file variables
metadata_row = 2    # Row in metadata file containing unique metadata names
sample_col = 0      # Column in metadata file containing unique sample names
sample_row = 0      # Row in ASV table files containing unique sample names
ASV_col = 0         # Column in metadata file containing unique ASV names

# Types of sequence data
# Genes sequenced and their corresponding max meaningful taxonomic level
genes = {
    '16S'   : 7,    # Tax level 7 is species level for 16S
    '18S'   : 11    # Tax level 11 is species level for 18S
    }
# Types of sequences (templates)
templates = ['DNA', 'cDNA']

# Variables to plot
varlist = ['day_num', 'horizon_code', 'sample_mat_code',
           'coord_lat', 'coord_lon', 'temperature', 'pH', 'bulk_density',
           'salinity_insitu', 'salinity', 'SPM', 'DOC', 'Chl', 'Phaeo',
           'FoFa', 'PAM', 'SedLoad', 'bact_cell_ct', 'CTC',
           'bact_active_cell_ct', 'diatom_ct', 'phyto_ct_all',
           'phyto_ct_other', 'pEPS', 'POC', 'nitrogen', 'CN']

# Samples to exclude
otherlist = ['space', 'Blank', 'EL']

# Colormap to use
cmap = 'PRGn'

# Significance cutoff for plot:
# Correlation coefficients with p values greater than this are excluded from
# the plot so that only "significant" correlations are shown
p_cutoff = 0.05

# HTML code for the for the HTML files
# Header links
linkhead = '''        
<h1>Rotten Ice Bio Data Analysis Navigation</h1>
<b>ASV bar plots: </b>
<a href = "http://faculty.weber.edu/cariefrantz/RottenIce/ASV-barplots/
16S_barplots_L1_viridis.html">Link</a><br />
<b>Algae ID bar plots: </b>
<a href = "http://faculty.weber.edu/cariefrantz/RottenIce/
AlgaeID_barplots_L15_viridis.html">Link</a><br />
<b>Beta diversity clustered heatmaps: </b>
<a href = "http://faculty.weber.edu/cariefrantz/RottenIce/B-div-heatmaps/
heatmaps_16S_all.html">Link</a><br />
<b>Metadata vs. taxonomic group correlation analysis: </b>
<a href = "spearman_16S-DNA.html">16S-DNA</a>  |  
<a href = "spearman_16S-cDNA.html">16S-cDNA</a>  |  
<a href="spearman_18S-DNA.html">18S-DNA</a>  |  
<a href="spearman_18S-cDNA.html">18S-cDNA</a><br />
'''
# Subheading text (Part II after the level number)
subhead2 = '''
showing correlation between biogeochemical variables and 
relative abundance of major taxonomic groups.
Values shown are Spearman correlation values calculated using the 
<a href="https://docs.scipy.org/doc/scipy/reference/generated/
scipy.stats.spearmanr.html">scipy.stats.spearmanr package</a>
(1.4.1) for python. 
Only significant correlations (p <=0.05) are shown. 
Positive numbers indicate a positive correlation, 
negative numbers indicate a negative correlation; 
larger numbers indicate stronger correlations. 
The top 20 most abundant taxonomic groups (across all samples)
for each level were analyzed. 
Analysis done by C. Frantz, June 2020.</p>
'''

####################
# FUNCTIONS
####################

def loadFile(text, sep=',', header_lines=0, index_col=0, 
             initialdir = os.getcwd()):
    '''Load file containing a table'''
    # User input dialog to locate the metadata file
    root = Tk()                     # Opens user input window
    root.filename = filedialog.askopenfilename(initialdir = initialdir,
                                               title = text)
    filename = root.filename
    root.destroy()                  # Closes the user input window
    dirPath = os.path.dirname(filename)     # Directory
    print ('Loading ' + filename + '...')
    # Read in and return data
    data = pd.read_csv(filename, sep = sep, header = header_lines,
                       index_col = index_col)
    return data, dirPath


def genTaxlist(tax_vals):
    '''Returns the unique taxonomic assignments in a taxonomy list'''
    taxlist = list(set(tax_vals))
    taxlist.sort()
    return taxlist


def levelCols(lmax):
    '''Generates the set of column names for all levels up to a given level'''
    cols = ['L' + str(i) for i in np.arange(1,lmax+1)]
    return cols


def formatASVTable(data):
    '''Prepares ASV table for analysis'''
    # Delete samples to ignore
    samples = list(data.columns)[0:-1]  # Pull list of samples from header
    samples = [sample for sample in samples
               if not any(other in sample for other in otherlist)]
    data = data.loc[:, samples + ['taxonomy']]
    
    # Normalize the data to reflect fractional abundance
    sums = data[samples].sum(axis = 0)
    data.loc[:,samples] = data[samples]/sums
    
    # Format taxonomy list to read better
    for i in np.arange(15):
        delstr = 'D_'+str(i)+'__'
        data.loc[:,'taxonomy'] = data['taxonomy'].str.replace(delstr, '')    
    
    # Break listed taxonomy into taxonomy at each taxonomic level
    taxlist = genTaxlist(data['taxonomy'])
    bar = Bar('Parsing taxonomy...', max = len(taxlist)) # Progress bar
    for value in taxlist:
        splitlist=[value]
        # Get list of levels
        if '; __' in value:
            splitlist = value.split('; __')
        if '; ' in value:
            splitlist = value.split('; ')
        # Fix last level if needed
        if splitlist[-1]:
            if splitlist[-1][-1] == ';':
                splitlist[-1] = splitlist[-1][0:-1]
        else: splitlist = splitlist[0:-1]
        # Fill in taxonomy value at each level    
        for i in np.arange(len(splitlist)):
            data.loc[data['taxonomy']==value, 'L'+str(i+1)] = splitlist[i]
        bar.next()  # Advance progress bar
    bar.finish()    # Stop progress bar
    
    # Get rid of nans in taxonomy levels
    end_level = 14
    cols = levelCols(end_level)
    data.loc[:, cols] = data[cols].replace(np.nan, '')
    
    return data, samples


def condenseDataset(data, level, samples):
    '''Condenses the dataset at a given level to combine ASVs with the same
        taxonomy call '''
    # Generate list of unique taxonomic values at the selected level
    cols = levelCols(level)
    taxnames = genTaxNames(data[cols])
    data['taxonomy'] = taxnames
    taxlist = genTaxlist(taxnames)
    
    # Set up a new dataframe to hold the condensed data
    condData = pd.DataFrame(index = taxlist, columns = samples)
    
    # Loop through taxonomy list and find combined ASV counts for each group
    for value in taxlist:
        # List all ASV ids that were assigned the given taxonomy value
        asv_ids = list(data.loc[data['taxonomy']==value].index)
        if len(asv_ids)==1:
            condData.loc[value, samples] = data.loc[asv_ids, samples].values
        else:
            sums = data.loc[asv_ids, samples].sum(axis = 0)
            condData.loc[value, samples] = sums
                
    return condData


def genTaxNames(taxcols):
    '''Generates combined taxonomy names from individual levels'''
    newtax = ['>'.join(list(taxcols.loc[i])) for i in list(taxcols.index)]
    return newtax


def topASVs(data, level, samples):
    '''Finds the 20 most abundant ASVs at the given taxonomic level'''
    # Add up abundance for all ASVs assigned the same taxonomy at given level
    ds = condenseDataset(data, level, samples)
    ds['sums'] = ds[samples].sum(axis=1)
    # Return the 20 most abundant ASVs
    if len(ds.index) >= 20:
        ds = ds.nlargest(20, 'sums')
    else: ds = ds.nlargest(len(ds.index), 'sums')
    ds = ds[samples].sort_index()
    return ds


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
    

def buildHeatmap(ax, corrtable):
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
    ax.set_title(gene + '-' + template + ' L' + str(level))
    
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
    metadata, directory = loadFile('Select metadata table (csv)',
                                   header_lines = [metadata_row],
                                   index_col = sample_col)
    metadata = metadata.dropna(how = 'all')
    metadata = metadata.replace('na', np.nan)
    # Build matrix of metadata values
    metadata = metadata[varlist]
    
    # Import ASV tables
    asv_tables = {}
    for gene in genes:
        for template in templates:
            dset = gene + '-' + template
            data, directory = loadFile('Select ' + dset + ' ASV table (csv)',
                                       initialdir = directory,
                                       header_lines = [sample_row],
                                       index_col = ASV_col)
            asv_tables[dset] = data
            
    # Build heatmaps
    # Loop through each dataset
    for dset in asv_tables:
        print('*******************************************************/n'
              'Determining Spearman correlations for ' + dset + ' dataset/n'
              '*******************************************************')
        
        # Format data
        data, samples = formatASVTable(asv_tables[dset])
        
        # Set up figure grid
        print('Preparing figure grid...')
        maxlevel = genes[dset[:3]]
        fig, axes = plt.subplots(maxlevel, 1, figsize = (20, maxlevel*8))
        
        # Perform correlation calculations
        # and build heatmaps for each taxonomic level
        bar = Bar('Calculating at each taxonomic level...',
                  max = maxlevel) # Progress bar
        # Loop through each taxonomic level
        for level in np.arange(1,maxlevel+1,1):
            # Condense the dataset to the level and return top ASVs
            ds = topASVs(data, level, samples)
            # Calculate Spearman correlation and return significant values
            corrtable = buildCorrTable(metadata, ds, samples)
            # Build heatmap
            ax = axes[level-1]
            im = buildHeatmap(ax, corrtable)
            texts = annotateHeatmap(im, fontsize = 6)
            bar.next()
        bar.finish()
        
        # Save plots
        print('Saving ' + dset + ' figure files...')
        fig.tight_layout()
        fig.savefig('spearman-heatmap_' + gene + '-' + template + '.svg',
                    transparent = True)

        # Generate HTML for display page
        head = ("<h1>Factors related to community changes: "
                + gene + " " + template + " sequences</h1>")
        subhead1 = '<p>Heatmaps at taxonomic levels 1-' + str(level)
        imghtml = ('<p><img src = "spearman-heatmap_' + gene + '-'
                   + template + '.svg" alt = "Spearman Heatmap ' + gene
                   + ' ' + template + '"></p>')
        html = linkhead + head + subhead1 + subhead2 + imghtml
        hfile = open('spearman_' + gene + '-' + template + '.html', 'w')
        hfile.write(html)
        hfile.close()
