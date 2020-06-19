#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""Rotten Ice Project: WebTaxBarplotsStacked

Created on Sat May 16 13:30:06 2020
@author: cariefrantz
@project: RottenIce

BUILDS INTERACTIVE (BOKEH) TAXONOMIC BARPLOTS: SEQS VERSION
This script builds an HTML webpage with interactive taxonomic bar plots.
It builds stacked plots displaying absolute and relative sequence abundances
for both DNA and cDNA, enabling direct comparisons of DNA vs. cDNA results.

Each taxonomic level is saved as a seperate webpage.

This version of the WebTaxBarplots script was used for sequencing data.
By modifying variables at the top of the code, it can easily be adapted for
any paired sets of OTU tables.

There is another version (WebTaxBarplotsSimple) that was used for
visualizing phytoplankton taxonomy data from microscopy counts.
It produces side-by-side plots and is best for a low number of samples.

This script was created as part of the Rotten Ice Project


Arguments:  None

Requirements:      
    ASV tables (csv) for every sequence dataset analyzed
        where rows = ASVs, columns = samples
        header row and index column are specified in the variables below

Example in command line:
    python WebTaxBarplotsStacked.py

Dependencies Install:
    sudo apt-get install python3-pip python3-dev
    pip install tkinter
    pip install progress
    pip install pandas
    pip install numpy
    pip install matplotlib
    pip install math


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
from tkinter import *
from tkinter import filedialog
from progress.bar import Bar

import pandas as pd
import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt

from bokeh.io import output_file, show, export_svgs
from bokeh.layouts import row, column
from bokeh.plotting import figure
from bokeh.models import FactorRange, Div


####################
# VARIABLES
####################

# DATASET-SPECIFIC VARIABLES
# (Modify these to match the dataset variables)

# Default directory
directory = '/'

# Input file variables
sample_row = 0      # Row in ASV table files containing unique sample names
OTU_col = 0         # Column in metadata file containing unique ASV names

# List of major sample groups and 'others' for sample grouping by prefix
# Sample names with the prefixes on the left
#   will be given the category label on the right
groups = {
    'M'     : 'May',
    'JN'    : 'June',
    'JY10'  : 'July 10',
    'JY11'  : 'July 11',
    'Other' : 'Other'
    }
# Delimiter splitting the sample group prefix from the sample name
grpdelim = '-'   
# Sample names including any of these strings will be grouped in the 'Other'
# category and will be excluded from the focus taxa analysis
otherlist = ['space', 'Blank', 'EL']

# List of taxonomic values to reassign,
# with their taxonomic level and new classification (after ':')
taxReassignList = {
    ('Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
     + 'Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
     + 'Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
     + 'Ambiguous_taxa; Ambiguous_taxa; D_14__')
            : 'Unassigned'}

# List of genes sequenced with maximum meaningful taxonomic level
# Each of these will be on its own page
genes       = {
    '16S'   : {
        'max_level'     : 7
        },
    '18S'   : {
        'max_level'     : 11
        }
    }

# PLOT PROPERTIES VARIABLES
# (Modify these to tweak how the plots look)
plttitletxt = ('Rotten ice project Illumina sequencing data, '
               'processed by C. Frantz May 2020 using pipeline-finished, '
               'clean ASV tables')

# Colors
cmap = 'viridis'    # Change this to change the default colormap
cmaptext = ''''Choose a colormap, e.g.,
viridis, plasma, prism, nipy_spectral, jet, gist_earth, etc.
See https://matplotlib.org/tutorials/colors/colormaps.html for a full list. 
> '''

# Number of taxa plotted
maxtaxa = 100   # Max number of taxa to show in a plot
fcutoff = 0.1   # Taxa with relative abundance greater than this value
                # in any sample in the dataset will be displayed
# Plot dimensions
pltht = 200     # Height of the plot box
pltw = 1500     # Width given to the actual plot
minpadding = 5  # Padding around plot
minystep = 5    # Number of y axis ticks (~5)
pltbottompad = 100 # Padding given to x axis label
keymaxwidth = 925 # Max width needed for key
keyvalht = 24   # Height taken up by each key entry
keytoppad = 40  # Height taken up by key title

# Plot details
# The types of plots to generate for each dataset
plttype = {
    'abs' : {
        'title'     : 'absolute',
        'data'      : 'pltdata',
        'xaxis'     : True,
        'ymax'      : 'auto',
        'widgets'   : 'hover,box_zoom,reset',
        'legend'    : False
        },
    'rel' : {
        'title'     : 'relative',
        'data'      : 'normpltdata',
        'xaxis'     : False,
        'ymax'      : 1,
        'widgets'   : 'hover,box_zoom,wheel_zoom,pan,reset',
        'legend'    : False
        }
    }

# All of the plots to draw in order and what goes in them
pltlist = {
    1 : {
        'template'  : 'DNA',
        'type'      : 'abs',
        'xlabelloc' : 'above',
        'ylabel'    : 'DNA ASV count \u27F6',
        'invert'    : False
        },
    2 : {
        'template'  : 'DNA',
        'type'      : 'rel',
        'xlabelloc' : 'below',
        'ylabel'    : 'DNA Relative abundance \u27F6',
        'invert'    : False
        },
    3 : {
        'template'  : 'cDNA',
        'type'      : 'rel',
        'xlabelloc' : 'below',
        'ylabel'    : '\u27F5 cDNA Relative abundance',
        'invert'    : True
        },
    4 : {
        'template'  : 'cDNA',
        'type'      : 'abs',
        'xlabelloc' : 'below',
        'ylabel'    : '\u27F5 cDNA ASV count',
        'invert'    : True
        }
    }


# PROGRAM VARIABLES
# (Modifying these may cause the program to break... proceed with caution)
# Data types
templates   = ['DNA', 'cDNA']


####################
# FUNCTIONS
####################
# Listed alphabetically by name

def buildLegend(taxlist, colors):
    '''Draws a legend to add to the figure'''
    # Generate dummy plot to hold legend
    p0 = figure(height = len(taxlist) * keyvalht + keytoppad)
    items = []
    for i in np.arange(len(taxlist)):
        items.append(p0.square(1,1, size = 15, color = colors[i],
                               legend_label = taxlist[i]))
    
    # Add legend
    p0.legend.location = 'top_left'
    p0.legend.title = 'Key'
    p0.legend.label_text_font_size = '9px'
    p0.legend.spacing = 0
    p0.legend.margin = 0
    
    # Make everything else invisible
    for item in items: item.visible=False
    p0.xaxis.visible = False
    p0.yaxis.visible = False
    p0.xgrid.visible = False
    p0.ygrid.visible = False
    p0.outline_line_width = 0
        
    return p0


def buildNavDiv(cmap):
    '''Generate the HTML header information that includes links to all of the
        plots generated in this script'''
    text = '<p>Plot Navigation  ---   '
    for gene in genes:
        max_level = genes[gene]['max_level']
        text = text + '<b>' + gene + ':</b>  '
        for level in np.arange(1,max_level+1):
            text = (text + ' <a href = "' + gene + '_barplots_L' + str(level)
                    + '_' + cmap + '.html' + '">L' + str(level) + '</a>  /')
        text = text[0:-1] + '  ---  '
    return Div(text = text[0:-7])


def buildPlot(ds, pltnum, xrange, colorlist):
    '''Builds inidividual plots based on predefined plot formatting variables
        and read-in plot data and formatting information'''
    
    # Get information about the plot being built from the initial variables
    pt = pltlist[pltnum]['type']
    title = ds['title'] + ' ' + plttype[pt]['title'] + ' level ' + str(level)
    # print(title + '...')
    
    # Build figure
    figoptions = {
        'x_range'             : xrange, 
        'x_axis_location'     : pltlist[pltnum]['xlabelloc'],
        'y_axis_label'        : pltlist[pltnum]['ylabel'],
        'plot_height'         : pltht + minpadding*2, 
        'plot_width'          : pltw + minpadding*2,
        'min_border_top'      : minpadding,
        'min_border_left'     : minpadding,
        'min_border_right'    : minpadding,
        'min_border_bottom'   : minpadding,
        'toolbar_location'    : 'right',
        'tools'               : plttype[pt]['widgets'], 
        'tooltips'            : '@samples $name: @$name'
        }
    p = figure(**figoptions)
    
    # Add the colorbars
    pltdata = ds[plttype[pt]['data']]
    taxlist = list(pltdata)[1:]
    p.vbar_stack(taxlist, x='samples', width = 0.9, color = colorlist, 
                         source=pltdata)
    
    # Format y axis
    # Determine the max value for the plot
    if isinstance(plttype[pt]['ymax'], (int, float)):
        maxval = plttype[pt]['ymax']
    else: 
        maxval = roundMaxVal(ds['maxval'])
    # Set the y axis range
    if pltlist[pltnum]['invert'] == False :
        p.y_range.start = 0
        p.y_range.end = maxval
    else:
        p.y_range.start = maxval
        p.y_range.end = 0
    
    # Format x axis
    p.axis.minor_tick_line_color = None
    p.x_range.range_padding = 0.02
    p.xgrid.visible = False
    if plttype[pt]['xaxis'] == True:
        p.xaxis.major_label_orientation = math.pi/2
    else: p.xaxis.visible = False
    
    # Add legend
    if plttype[pt]['legend'] == True:
        legend = buildLegend(taxlist, colorlist)
        p.add_layout(legend, 'right')
        
    # Save figure
    p.output_backend = 'svg'
    export_svgs(p, filename = 'barplot_' + title + '.svg')
    return p


def colorsFromCmap(n, cmap):
    '''Generate the set of colors to use from the chosen colormap'''
    colormap = plt.get_cmap(cmap)
    colorset = colormap(np.linspace(0,1,n))
    colors=[]
    for color in colorset: colors.append(matplotlib.colors.rgb2hex(color))
    return colors         

    
def colorMap2Tax(taxset, max_level, cmap):
    '''Map colors to taxonomy values for each level
        (this keeps colors consistent across levels for taxonomic groups)
        taxset is the full taxonomy breakdown
        cmap is the string for the color map to use'''
    
    cols_all = levelCols(max_level)
    
    colnames = ['color-' + col for col in cols_all]
    for c in colnames:
        taxset[c] = ''
    
    # Find the greatest level with <= 256 unique values
    print('Finding level to use as colormap index...')
    uniques=[]
    for level in np.arange(1, max_level+1):
        cols_sub = levelCols(level)
        uniquetax = list(set(genTaxNames(taxset[cols_sub])))
        uniques.append(uniquetax)
        if len(uniquetax) > 256: break
    cmaplevel = level-1
    
    # Generate colors for each OTU at the highest taxonomic level
    # Generate tax names for each group at the level
    cols_clevel = levelCols(cmaplevel)
    uniquetax = list(set(genTaxNames(taxset[cols_clevel])))
    print('Using Level ' + str(cmaplevel) + ' for colormap index ('
          + str(len(uniquetax)) + ' unique values)')
    uniquetax.sort()
    colorlist = colorsFromCmap(len(uniquetax), cmap)
    colordict = dict(zip(uniquetax, colorlist))
    
    # Fill in the values for the max level and all lower levels
    taxset['taxonomy'] = genTaxNames(taxset[cols_clevel])
    cols_low = np.setdiff1d(levelCols(max_level), levelCols(cmaplevel-1))
    print('Assigning colors to L' + str(cmaplevel) + '-L'
          + str(max_level) + '...')
    cols_low_clr = ['color-' + col for col in cols_low]
    for i in list(taxset.index):
        taxset.loc[i][cols_low_clr] = colordict[taxset.loc[i]['taxonomy']]
    
    # Fill in the values for the higher levels
    print('Assigning colors to L1-L' + str(cmaplevel-1))
    # Loop through remaining levels.
    # At each level, pick the mid value for each set.
    for level in range(cmaplevel-1, 0, -1):
        cols_sub = levelCols(level)
        taxlist = genTaxNames(taxset[cols_sub])
        # Get unique list of tax values at this level
        # Loop through each unique value
        for taxval in uniques[level-1]:
            # Find all color values in the next level belonging to this tax set
            rows = [i for i in np.arange(len(taxlist))
                    if taxlist[i] == taxval]
            colorlist = list(
                taxset.iloc[rows]['color-L' + str(level+1)].values)
            # Determine the middle color
            midclr = colorlist[int(len(colorlist)/2)]
            # Assign this color to all with this tax value
            taxset.iloc[rows, max_level + level-1] = midclr
            
    # Return the colormap
    cols_all_clr = ['color-L' + str(i) for i in np.arange(1,max_level+1)]        
    colormap = taxset[cols_all_clr]
   
    return colormap   


def dataset2Dict(data, taxlist, samples, sgroups):
    '''Generate the dictionary containing plot data and information'''
    # Build the dict & color list
    pltdata = {'samples' : sgroups}
    colorlist = []
    # Read in the data and save colors
    for value in taxlist:
        rows = data.loc[data['taxonomy']==value]
        if rows.shape[0] == 1:
            pltdata[value] = list(rows[samples].to_numpy()[0])
        else:
            mat = rows[samples].to_numpy()
            pltdata[value] = list(np.sum(mat, axis=0))
        colorlist.append(rows['colors'][0])
    return pltdata, colorlist


def dataset2Plot(dataset, level, colormap):
    '''Prepare the data to use for each plots'''
    data = dataset['groupeddata']
    data['colors'] = colormap['color-L' + str(level)]
    # Reset taxonomic ID based on the selected level
    cols = levelCols(level)
    data['taxonomy'] = genTaxNames(data[cols])
    
    # Extract data, samples, samplegroups
    samples = dataset['samples']
    sgroups = groupSamples(samples)
    dataset['sampleGroups'] = sgroups
    
    # Create normalized data
    normdata = data.copy().fillna(0)
    maxsum = 0
    for sample in samples:
        sumvals = np.sum(normdata[sample])
        if sumvals > 0 :    normdata[sample] = normdata[sample].values/sumvals
        else:               normdata[sample].values[:] = 0
        if sumvals > maxsum: maxsum = sumvals
    
    # Generate list of unique taxonomic values at the selected level
    taxlist = genTaxList(data['taxonomy'].values)
    # Create corresponding datasets for plotting
    [dataset['pltdata'], colorlist] = dataset2Dict(
        data.fillna(0), taxlist, samples, sgroups)
    [dataset['normpltdata'], colorlist] = dataset2Dict(
        normdata, taxlist, samples, sgroups)
    dataset['maxval'] = maxsum
        
    return dataset, colorlist


def condenseDataset(data, level, samples):
    '''Condenses the dataset at a given level to combine ASVs with the same
        taxonomy call '''
    print('Condensing dataset at L' + str(level) + '...')
    # Generate list of unique taxonomic values at the selected level
    cols = levelCols(level)
    taxnames = genTaxNames(data[cols])
    data['taxonomy'] = taxnames
    taxlist = genTaxList(taxnames)
    
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

def fileGet(title, directory = directory):
    '''Gets ASV table file from user input'''
    # Import csv file (remove extra headers)
    root=Tk()
    filename=filedialog.askopenfilename(initialdir=directory, title = title,
                                        filetypes = [('CSV', '*.csv')])
    dirPath = os.path.dirname(filename)     # Directory
    root.destroy()
    return filename, dirPath


def fileRead(filename, max_level):
    '''Reads in and formats each ASV table'''
    
    # Read in and format the table
    print('Loading ' + filename + '...')
    data = pd.read_csv(filename, sep = ',',
                       header = sample_row, index_col = OTU_col)
    
    # Get sample list
    samples = list(data.columns)[0:-1]
       
    # Reclassify any values with assignments in the taxReassignList
    for val in list(taxReassignList):
        data.loc[data['taxonomy']==val, 'taxonomy'] = taxReassignList[val]
    
    # Format taxonomy list to read better
    print('Formatting taxonomy...')
    for i in np.arange(15):
        delstr = 'D_'+str(i)+'__'
        data['taxonomy'] = data['taxonomy'].str.replace(delstr, '')        

    # Break taxmap into levels
    taxlist = genTaxList(data['taxonomy'])
    bar = Bar('', max = len(taxlist))
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
    end_level = len(data.columns)-len(samples)-1
    cols = levelCols(end_level)
    data[cols] = data[cols].replace(np.nan, '')
    
    # Condense dataset based on the max level
    # data = datasetCondense(data, max_level, samples)

    return data, samples


def findFocusTaxa(dslist):
    'Identify the set taxa to include in the plot'
    focusTaxa = []
    for dataset in dslist:
        # Find the most abundant ASVs
        focusTaxa = focusTaxa + findAbundantTaxa(datasets[dataset]['data'],
                                                 datasets[dataset]['samples'])
    focusTaxa = np.unique(focusTaxa)
    for dataset in dslist:
        datasets[dataset]['focusTaxa'] = focusTaxa
        
    return datasets
        

def findAbundantTaxa(data, samples):
    ''''Find the most abundant ASVs in the dataset
        Data is the set of data for the given level
        Samples is the list of samples
        Only base the list off of values in the non-Other samples'''
    usecols = [s for s in samples
               if not any(other in s for other in otherlist)]
    
    # Find all ASVs with rel abundance > cutoff
    print('Finding any ASVs above min fraction cutoff...')
    normdata = data[usecols] / np.sum(data[usecols], axis = 0)
    rows = []
    for datarow in normdata.index:
        vals = normdata.loc[datarow, usecols]
        if any(val > fcutoff for val in vals):
            rows.append(datarow)
    
    # Find the most abundant ASVs
    print('Finding list of most abundant ASVs..')
    n = maxtaxa - len(rows)
    data['totals'] = np.sum(data[usecols], axis = 1)
    data = data.sort_values('totals', ascending = False)
    rows = rows + list(data.iloc[0:n].index)
    abundantTaxa = data.loc[rows,'taxonomy']
    abundantTaxa = list(np.unique(abundantTaxa))
    
    return abundantTaxa


def genTaxList(tax_vals):
    '''Returns the unique taxonomic assignments in a taxonomy list'''
    taxlist = list(set(tax_vals))
    taxlist.sort()
    return taxlist


def genTaxNames(taxcols):
    '''Generates combined taxonomy names from individual levels'''
    newtax = ['>'.join(list(taxcols.loc[i])) for i in list(taxcols.index)]
    return newtax


def groupTaxa(rawdata, focusTaxa):
    '''Merge taxa into new groups based on the focus taxa and "others"'''
    data = rawdata.copy()
    excludedOthers = []
    
    # Loop through each level
    bar = Bar('', max = 14) # Progress bar
    for L in np.arange(1,15):
        excludedOthers = excludedOthers + ['\nAssigning L' + str(L) + ':\n']
        # Get full taxonomy of any taxonomic groups
        #     with a taxonomic value at this level
        # Find non-empty values in the taxonomy column for this level
        indices = [i for i in data.index if data.loc[i]['L'+str(L)]]
        # Generate taxonomy names at this level
        data['taxonomy'] = genTaxNames(data[levelCols(L)]) 
        # Generate list of unique taxonomic names at this level                 
        taxlist = genTaxList(data.loc[indices]['taxonomy'])
        # Remove any empty values
        taxlist = list(filter(None, taxlist))                               
        # If the list isn't empty...
        if taxlist:
            # Loop through each unique taxonomy name
            for value in taxlist:
                # Get the full taxonomic name for everything in the category
                indices = [i for i in data.index
                           if data.loc[i]['taxonomy']==value]
                taxvals = np.unique(list(rawdata.loc[indices, 'taxonomy']))
                # If no items in the tax list are members of the focus taxa,
                #   call all values 'other' for that level
                if not any(item in taxvals for item in focusTaxa):
                    excludedOthers = (
                        excludedOthers + [str(val) for val in taxvals])
                    data.loc[indices, 'L'+str(L)] = 'Other'
                    data.loc[indices, 'L'+str(L+1):'L'+str(14)] = ''
        bar.next()
    bar.finish()
    return data, excludedOthers
                

def groupSamples(samples):
    '''Group samples as outlined in the initial variables'''
    sampleset=[]
    for sample in samples:
        if any(other in sample for other in otherlist):
            sampleset.append(('Other', sample))
        else:
            ssplit = sample.split('-')
            sampleset.append((groups[ssplit[0]], '-'.join(ssplit[-2:])))  
    return sampleset


def levelCols(lmax):
    '''Generates the set of column names for all levels up to a given level'''
    cols = ['L' + str(i) for i in np.arange(1,lmax+1)]
    return cols


def roundMaxVal(maxval):
    '''Find a maximum value for an axis that is rounded
        to the nearest whole step for plot ticks'''
    minstepsize = maxval/minystep
    logval = math.floor(math.log10(minstepsize))
    roundval = minystep * math.ceil(minstepsize/10**logval)*10**logval
    return roundval
#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    
    datasets = {}
    # Get the ASV table data files (csv) to use
    for gene in list(genes):
        for template in templates:
            dsetID = gene + '_' + template
            datasets[dsetID]={
                'gene'      : gene,
                'template'  : template,
                'title'     : gene + ' ' + template,
                'max_level' : genes[gene]['max_level']
                }
            datasets[dsetID]['filename'], directory  = fileGet(
                'Select ' + datasets[dsetID]['title'] + ' OTU Table CSV',
                directory)
    
    # Get the colors (disable to use default colors)
    # cmap = input(cmaptext)
    
    # Load the raw data and taxonomy
    print('\n***************************\n'
          'Loading and formatting data\n'
          '***************************')
    for ds in list(datasets):
        [datasets[ds]['data'], datasets[ds]['samples']] = fileRead(
            datasets[ds]['filename'], datasets[ds]['max_level'])
    
    
    # Build plot sets for each gene
    # Loop through each gene (seperate set of pages)
    for gene in list(genes):
        print('\n***************************\n'
              'Processing ' + gene + ' data\n'
              '***************************')
        # Identify focus taxa
        dslist = [dataset for dataset in datasets
                  if datasets[dataset]['gene']==gene]
        datasets = findFocusTaxa(dslist)
        
        # Loop through each dataset for each gene
        for ds in dslist:
            # Rename all members of taxonomic groups that are in the same
            # L-1 taxonomic group as a focus taxon 'Other'
            print('Grouping ' + ds + ' ASVs at each taxonomic level...')
            groupeddata, others = groupTaxa(
                datasets[ds]['data'], datasets[ds]['focusTaxa'])
            datasets[ds]['groupeddata'] = groupeddata
            datasets[ds]['excludedOthers'] = others

        # Get the colors for the plot
        max_level = datasets[ds]['max_level']
        colormap = colorMap2Tax(
            datasets[ds]['groupeddata'][levelCols(max_level)].copy(),
            max_level, cmap)
    
        # Build the plots
        for level in np.arange(1,max_level+1):
            print('Generating Level ' + str(level) + ' plots...')
            
            # Set up HTML file
            output_file(
                gene + '_barplots_L' + str(level) + '_' + cmap + '.html',
                title = datasets[ds]['title'] + ' barplots L' + str(level))
            
            # Generate the plot data for each plot
            for ds in dslist:
                [datasets[ds], colorlist] = dataset2Plot(datasets[ds].copy(),
                                                       level, colormap)
            
            # Build the plots
            plots = []
            # Loop through each plot in the plot list
            bar = Bar('', max = len(pltlist))
            for pltnum in pltlist:
                # Retrieve the corresponding dataset
                dataset = datasets[gene + '_' + pltlist[pltnum]['template']]
                
                # link x axes in the plot so all plots are on same x scale
                if not plots:
                    sampleset = dataset['sampleGroups']
                    xrange = FactorRange(*sampleset)
                else: xrange = plots[0].x_range
                
                # Build the plot
                p = buildPlot(dataset, pltnum, xrange, colorlist)  
                plots.append(p)
                bar.next()
            bar.finish()
            
            # Build the legend
            legend = buildLegend(list(dataset['pltdata'])[1:], colorlist)
    
            # Build the HTML file
            # HTML page navigation bar
            navdiv = buildNavDiv(cmap)
            titlediv = Div(text = '<p><b>' + template +
                           ' Data, taxonomic level ' + str(level) + ':</b> '
                           + plttitletxt + '<br />Top set is '
                           + datasets[list(datasets)[0]]['title']
                           + ', bottom set is '
                           +  datasets[list(datasets)[1]]['title'] )
            grid = column([navdiv, titlediv, row([column(plots), legend])])
            # Save and show HTML file
            show(grid)
       
