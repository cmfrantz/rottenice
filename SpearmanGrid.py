# -*- coding: utf-8 -*-
"""
Created on Sun May 31 15:24:03 2020

@author: cariefrantz


Script calculates Spearman correlation coefficients relating metadata to taxonomic group abundance
It was created as part of the Rotten Ice Project
It uses the following files:
-Sample metadata table (csv)
-ASV tables (csv) for every sequence dataset analyzed

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

import os
import pandas as pd
import numpy as np
from tkinter import filedialog
from tkinter import *
from progress.bar import Bar
from scipy.stats import spearmanr
import matplotlib
from matplotlib import pyplot as plt


##########
#VARIABLES
templates = {
    '16S'   : {
        'levels'    : 7
        },
    '18S'   : {
        'levels'    : 14
        }
    }
seqtypes = ['DNA', 'cDNA']
varlist = ['day_num', 'horizon_code', 'sample_mat_code', 'coord_lat', 'coord_lon', 'temperature', 'pH', 'bulk_density', 'salinity_insitu', 'salinity', 'SPM', 'DOC', 'Chl', 'Phaeo', 'FoFa', 'PAM', 'SedLoad', 'bact_cell_ct', 'CTC', 'bact_active_cell_ct', 'diatom_ct', 'phyto_ct_all', 'phyto_ct_other', 'pEPS', 'POC', 'nitrogen', 'CN']
otherlist = ['space', 'Blank', 'EL']    # Samples to exclude
cmap = 'PRGn' # Colormap to use

linkhead = '''        
<h1>Rotten Ice Bio Data Analysis Navigation</h1>
<b>ASV bar plots: </b><a href = "http://faculty.weber.edu/cariefrantz/RottenIce/ASV-barplots/16S_barplots_L1_viridis.html">Link</a><br />
<b>Algae ID bar plots: </b><a href = "http://faculty.weber.edu/cariefrantz/RottenIce/AlgaeID_barplots_L15_viridis.html">Link</a><br />
<b>Beta diversity clustered heatmaps: </b><a href = "http://faculty.weber.edu/cariefrantz/RottenIce/B-div-heatmaps/heatmaps_16S_all.html">Link</a><br />
<b>Metadata vs. taxonomic group correlation analysis: </b>
<a href = "spearman_16S-DNA.html">16S-DNA</a>  |  <a href = "spearman_16S-cDNA.html">16S-cDNA</a>  |  <a href="spearman_18S-DNA.html">18S-DNA</a>  |  <a href="spearman_18S-cDNA.html">18S-cDNA</a><br />
'''

#%%

# Load file with a table
def loadFile(text, sep=',', header_lines=0, index_col=0):
    root = Tk()
    root.filename = filedialog.askopenfilename(initialdir = "/", title = text)
    filename = root.filename
    print ('Loading ' + filename + '...')   # Prints the filename
    root.destroy()                          # Closes the user input window
    data = pd.read_csv(filename, sep = sep, header = header_lines, index_col = index_col)  # Reads in all of the matrix data
    return data, filename


# Find unique taxonomic assignments in a taxonomy list
def genTaxlist(tax_vals):
    taxlist = list(set(tax_vals))
    taxlist.sort()
    return taxlist


# Generates the column names up to a given level
def levelCols(lmax):
    cols = ['L' + str(i) for i in np.arange(1,lmax+1)]
    return cols


# Read the distance matrix in from tab-seperated file and format
def formatASVTable(dset):
    [data, filename] = loadFile('Select ' + dset + ' ASV table (csv)...')
    
    # Delete samples to ignore
    samples = list(data.columns)[0:-1]
    samples = [sample for sample in samples if not any(other in sample for other in otherlist)]
    data = data[samples + ['taxonomy']]
    
    # Normalize the data to reflect fractional abundance
    sums = data[samples].sum(axis = 0)
    data[samples] = data[samples]/sums
    
    # Format taxonomy list to read better
    print('Formatting taxonomy...')
    for i in np.arange(15):
        delstr = 'D_'+str(i)+'__'
        data['taxonomy'] = data['taxonomy'].str.replace(delstr, '')    
    
    # Break taxmap into levels
    taxlist = genTaxlist(data['taxonomy'])
    bar = Bar('Formatting taxonomy...', max = len(taxlist))
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
            
        for i in np.arange(len(splitlist)):
            data.loc[data['taxonomy']==value, 'L'+str(i+1)] = splitlist[i]
            
        bar.next()
    bar.finish()
    
    # Get rid of nans in taxonomy levels
    end_level = 14
    cols = levelCols(end_level)
    data[cols]=data[cols].replace(np.nan, '')
    
    # Condense dataset
    data = condenseDataset(data, end_level, samples)
    
    return data, samples


# Generate the data to use for the plots
def condenseDataset(data, level, samples):

    print('Condensing dataset at L' + str(level) + '...')
    condData = pd.DataFrame()
    # Reset taxonomic ID based on the selected level
    cols = levelCols(level)
    data['taxonomy'] = genTaxNames(data[cols])
    
    # Generate list of unique taxonomic values at the selected level
    taxlist = genTaxlist(data['taxonomy'].values)
    
    # Combine values with same taxonomy value
    for value in taxlist:
        rows = data.loc[data['taxonomy']==value].copy()
        if rows.shape[0]==1:
            condData = condData.append(rows)
        else:
            mat = rows[samples].to_numpy()
            rows.loc[list(rows.index)[0], samples] = list(np.sum(mat, axis = 0))
            condData = condData.append(rows.iloc[0])
    
    # Get the new dataframe sorted
    condData = condData.reindex(columns = data.columns)
        
    return condData


# Generate taxonomy names
def genTaxNames(taxcols):
        # Taxcols are the taxonomy columns for all levels to join
        newtax = ['>'.join(list(taxcols.loc[i])) for i in list(taxcols.index)]
        return newtax


# Find top ASVs
def topASVs(data, level, samples):
    ds = condenseDataset(data, level, samples)
    ds['sums'] = ds[samples].sum(axis=1)
    if len(ds.index) >= 20:
        ds = ds.nlargest(20, 'sums')
    else: ds = ds.nlargest(len(ds.index), 'sums')
    ds = ds.set_index('taxonomy')
    ds = ds[samples].sort_index()
         
    return ds


# Build metadata vs. abundance table
def buildCorrTable(metadata, abundance_table, samples):
    corrtable = pd.DataFrame(data = None, index = abundance_table.index, columns = varlist)
    #corrtable = pd.DataFrame(data = None, index = varlist, columns = abundance_table.index)
    msamples = [sample[0:-2] for sample in samples]
    
    for group in abundance_table.index:
        for var in varlist:
            # Grab the list of values for each sample
            varvals = metadata.loc[msamples, var].to_numpy(dtype=float)
            counts = abundance_table.loc[group, samples].to_numpy(dtype=float)
            
            #Calculate Spearman coefficient
            corr, p = spearmanr(varvals, counts, nan_policy='omit')
            
            # If p value <= 0.05, assign the correlation a color on the diverging color scale and add it to the heatmap matrix
            if p <= 0.05:
                corrtable.loc[group, var] = corr
                
    corrtable = corrtable.astype(float)
                
    return corrtable
    
#%%
# Build the heatmap
def buildHeatmap(ax, corrtable):
    # Build the figure
    im = ax.imshow(corrtable, cmap = cmap, vmin = -1, vmax = 1)
    #im = ax.pcolor(corrtable, edgecolors = 'k', linewidths = 1, cmap = cmap, vmin = -1, vmax = 1)
    ax.set_xticks(np.arange(corrtable.shape[1]))
    ax.set_yticks(np.arange(corrtable.shape[0]))
    ax.set_xticklabels(corrtable.columns, rotation = 45, ha = 'right', rotation_mode = 'anchor')
    ax.set_yticklabels(corrtable.index)
    ax.set_title(template + '-' + seqset + ' L' + str(level))
    ax.tick_params(labelsize = 9)
    
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax = ax)
    cbar.ax.set_ylabel("Spearman's' Correlation", va = 'bottom')
    
    # Add grid
    ax.set_xticks(np.arange(corrtable.shape[1]+1)-0.5, minor = True)
    ax.set_yticks(np.arange(corrtable.shape[0]+1)-0.5, minor = True)
    ax.grid(which = 'minor', color = 'k', linestyle = '-', linewidth = 1)
    ax.tick_params(which = 'minor', bottom = False, left = False)
    
    return im
        
    
# Annotate heatmap
def annotateHeatmap(im, data=None, valfmt = "{x:.2f}", **textkw):
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

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color='w')
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts    

#%%
####
# MAIN FUNCTION
if __name__ == '__main__':
    
    # Import Metadata table
    metadata, filename = loadFile('Select metadata table (csv)...', sep = ',', header_lines = 2, index_col = 0)
    metadata = metadata.dropna(how = 'all')
    metadata = metadata.replace('na', np.nan)
    # Build matrix of metadata values
    metadata = metadata[varlist]
        
    # Build heatmaps
    for template in templates:
        
        # Import data
        for seqset in seqtypes:
            print('Importing ' + template + '-' + seqset + ' table...')
            templates[template][seqset], templates[template][seqset + '-samples'] = formatASVTable(template + '-' + seqset)
                            
            # Set up figure grid
            maxlevel = templates[template]['levels']
            fig, axes = plt.subplots(maxlevel, 1, figsize = (20, maxlevel*8))
            
            # Wrangle data and make heatmaps
            for level in np.arange(1,templates[template]['levels']+1,1):
                print('Building level ' + str(level) + ' figure...')
                # Condense the dataset to the level and return top ASVs
                ds = topASVs(templates[template][seqset], level, templates[template][seqset+'-samples'])
                # Calculate Spearman's correlation and return significant values
                corrtable = buildCorrTable(metadata, ds,  templates[template][seqset+'-samples'])
                # Build heatmap
                ax = axes[level-1]
                im = buildHeatmap(ax, corrtable)
                texts = annotateHeatmap(im, fontsize = 6)

            # Draw and save plots   
            fig.tight_layout()
            fig.savefig('spearman-heatmap_' + template + '-' + seqset + '.svg', transparent = True)

            # Create html
            head = "<h1>Factors related to community changes: " + template + " " + seqset + " sequences</h1>"
            subhead = '<p>Heatmaps at taxonomic levels 1-' + str(level) + ' showing correlation between biogeochemical variables and relative abundance of major taxonomic groups. Values shown are Spearman&#39;s correlation values calculated using the <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html">scipy.stats.spearmanr package</a> (1.4.1) for python. Only significant correlations (p <=0.05) are shown. Positive numbers indicate a positive correlation, negative numbers indicate a negative correlation; larger numbers indicate stronger correlations. The top 20 most abundant taxonomic groups (across all samples) for each level were analyzed. Analysis done by C. Frantz, June 2020.</p>'
            imghtml = '<p><img src = "spearman-heatmap_' + template + '-' + seqset + '.svg" alt = "Spearman Heatmap ' + template + ' ' + seqset + '"></p>'
            html = linkhead + head + subhead + imghtml
            
            hfile = open('spearman_' + template + '-' + seqset + '.html', 'w')
            hfile.write(html)
            hfile.close()
