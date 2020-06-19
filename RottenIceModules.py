#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""Rotten Ice Project: Modules

Created on Thu Jun 18 21:49:22 2020
@author: cariefrantz
@project: RottenIce

MODULES SHARED BETWEEN DIFFERENT SCRIPTS IN THE ROTTEN ICE PROJECT
This script contains a collection of functions that are used in multiple \
    scripts in the Rotten Ice project.

Arguments:  None

Requirements of each function are described in the documentation.

Example in external script (must be either in the PATH or same directory):
    import RottenIceModules

Dependencies Install:
    sudo apt-get install python3-pip python3-dev
    pip install tkinter
    pip install progress
    pip install numpy
    pip install pandas
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

import numpy as np
import pandas as pd
import matplotlib as plt
import math

from bokeh.io import output_file, show, export_svgs
from bokeh.plotting import figure

####################
# FUNCTIONS
####################
# Listed alphabetically by name

def buildBarPlot(plt_data, colorlist, y_max, title='', y_label = 'Abundance',
                 invert = False, show_x_ax = True, legend = False,
                 figoptions={}):
    '''Builds a bokeh interactive bar plot of taxonomy data

    Parameters
    ----------
    plt_data : dict
        Dictionary containing data for each bar to plot.
        First dict key should be 'samples' : [list of sample name strings]
        Other dict keys should be assigned taxonomy strings.
        Dict values should be lists of OTU counts for each sample (in order).
    colorlist : list of tuples
        List of colors to use for each bar. Length of colorlist should be \
            the same as the number of taxonomic values in the passed dict.
    y_label : str, optional)
        String to use for the y-axis label. The default is 'Abundance'.
    y_max : float or int
        Y axis maximum value. Use 1 for relative abundance plots. \
            determines it automatically from the data.
    title : str, optional
        Plot title. The default is ''.
    invert : Boolean, optional
        If invert is true, the plot will be mirrored (upside-down). \
            The default is False.
    show_x_ax : Boolean, optional
        Whether or not to display the x axis labels. The default is True.
    legend: Boolean, optional
        Whether or not to display a figure legend. The default is False.
    figoptions : dict, optional
        Dictionary containing optional figure formatting options for \
            bokeh.figure. The default is {}.

    Returns
    -------
    p   : bokeh figure
        Bokeh figure object containing the plot.

    '''
    # Create figure
    print('Building ' + title + ' plot...')
    p = figure(**figoptions)
    
    # Add bars corresponding to each taxon
    taxlist = list(plt_data)[1:]
    vbars = p.vbar_stack(taxlist, x = 'samples', width = 0.9, color = colorlist,
                 source = plt_data)
    
    # Set y axis scale
    if invert == False:
        p.y_range.start = 0
        p.y_range.end = y_max
    else:
        p.y_range.start = y_max
        p.y_range.end = 0
    
    # Set x axis properties
    p.axis.minor_tick_line_color = None
    p.x_range.range_padding = 0.05
    p.xgrid.visible = False
    if show_x_ax == True:
        p.xaxis.major_label_orientation = math.pi/2
    else:
        p.xaxis.visible = False
    
    # Add legend
    if legend == True:
        legend = genLegend(taxlist, vbars)
        p.add_layout(legend, 'right')                                           # Figure out genLegend
        
    # Save figure
    p.output_backend = 'svg'
    export_svgs(p, filename)
    return p
    

def buildNavDiv(cmap):
    '''Generate the HTML header information that includes links to all of the   Fill in
        HTML pages generated for this project'''


def colorsFromCmap(n, cmap):
    '''Generate the set of colors to use from the chosen colormap               Fill in

    Parameters
    ----------
    n : int
        DESCRIPTION.
    cmap : str
       The matplotlib colormap name to use for picking colors.

    Returns
    -------
    colors : type
        DESCRIPTION.
    '''
    colormap = plt.get_cmap(cmap)
    colorset = colormap(np.linspace(0,1,n))
    colors=[]
    for color in colorset: colors.append(plt.colors.rgb2hex(color))
    return colors         

    
def colorMap2Tax(taxset, cmap, max_level = 'Auto'):
    '''This function maps colors to taxonomy values for plots. \
        Doing this keeps colors consistent across levels for taxonomic groups.   

    Parameters
    ----------
    taxset : pandas.DataFrame
        Full taxonomy breakdown matrix with levels broken out.
    cmap : str
       The matplotlib colormap name to use for picking colors.
    max_level : int or 'Auto', optional
        The maximum taxonomic level to consider. The default is 'Auto'.

    Returns
    -------
    colormap : list of tuples
        List of colors (as tuples) to use for each taxon.
    '''
   
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
            # Find all color values in the next level belonging to this
            # taxonomic set
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


def condenseDataset(data, level, samples):
    '''This script condenses an ASV table with taxonomy to only include
        a subset of samples. It also combines all ASVs with the same taxonomic
        call at the level selected.

    Parameters
    ----------
    data : DataFrame
        ASV table dataset with taxonomy calls:
            index:      ASV IDs
            headers:    <samples> + 'taxonomy' + any others
    level : int
        Taxonomic level to condense the dataset to.
    samples : list
        Sample IDs matching <data> headers to keep in the condensed dataset.

    Returns
    -------
    condData : DataFrame
        Condensed DataFrame containing:
            index:      unique taxonomic calls at the given level
            header:     <samples>
            data:       sum of all ASVs for each taxonomic group and sample
    '''
    print('Condensing dataset at L' + str(level) + '...')
    # Generate list of unique taxonomic values at the selected level
    cols = levelCols(level)
    taxnames = genTaxNames(data[cols])
    data['taxonomy'] = taxnames
    taxlist = getUnique(taxnames)
    
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


def fileGet(title, directory = os.getcwd(), file_type = 'csv',
            header_row = 0, index_col = 0):
    '''This script imports data from a file chosen with a user input GUI.
    
    Parameters
    ----------
    title : str
        Text to put in the user input window.
    directory : str, optional
        The start directory to search in. The default is os.getcwd().
    file_type : str, optional
        The type of file. Options are currently only:
            'csv' (default)     comma-seperated file
            'tsv'               tab-seperated file
    header_row : int, optional
        The row of the file to capture for the data header. The default is 0.
    index_col : int, optional
        The column of the file to capture for the data index. The default is 0.

    Returns
    -------
    filename : str
        Selected file full filepath.
    dirPath : str
        Filepath for the directory containing the selected file.
    data : pandas.DataFrame
        DataFrame containing the indexed data read in from the file.
    '''
    
    # Define filetypes
    if file_type == 'csv':
        ftype = [('CSV', '*.csv')]
        sep = ','
    elif file_type == 'tsv':
        ftype = [('TSV'), '*.tsv, *.txt']
        sep = '\t'
    
    # Open user input file dialog to pick file
    root=Tk()
    filename=filedialog.askopenfilename(initialdir=directory, title = title,
                                        filetypes = ftype)
    dirPath = os.path.dirname(filename)     # Directory
    root.destroy()
    
    # Read file
    print('Loading ' + filename + '...')
    data = pd.read_csv(filename, sep = sep,
                       header = header_row, index_col = index_col)
    
    return filename, dirPath, data


def findAbundantTaxa(datasets, samples, fcutoff = 0.1, maxtaxa = 100):
    '''This script identifies the limited set of most abundant taxa to \
        include in downstream analysis (e.g., bar plots) from ASV tables of \
            sequence results of the same gene with the same samples \
                (e.g., 16S DNA & 16S cDNA for the same samples).

    Parameters
    ----------
    datasets : list of pandas.DataFrame
        List containing ASV table DataFrames with shared taxonomic values.
        DataFrames must include sample name headers as well as 'taxonomy' \
            header.
            
    samples : list of str
        List containing sample names (as str) to be included in the analysis.   # Make sure scripts using this pass sample names with 'others' excluded
        All sample names should also be headers in all included DataFrames.
            
    fcutoff : float (optional)
        If any taxa is greater than this value in abundance in any sample, \
            it will be included. Default is 0.1 (10%).
            
    maxtaxa : int (optional)
        The maximum number of abundant taxa to return.

    Returns
    -------
    abundantTaxa : list of str
        Alphabetically-sorted list of the most abundant taxa in the selected \
            samples.
    '''
    abundantTaxa = []
    # Search in each dataset
    for dataset in datasets:
        # Find the most abundant ASVs
        # Find all ASVs with rel abundance > cutoff
        print('Finding ASVs above min fraction cutoff...')
        # Normalize data to get relative abundance values
        normdata = dataset[samples] / np.sum(dataset[samples], axis = 0)
        rows = []
        for datarow in normdata.index:
            vals = normdata.loc[datarow, samples]
            if any(val > fcutoff for val in vals):
                rows.append(datarow)
        
        # Find the most abundant ASVs
        print('Finding most abundant ASVs..')
        n = maxtaxa - len(rows)
        dataset['totals'] = np.sum(dataset[samples], axis = 1)
        dataset = dataset.sort_values('totals', ascending = False)
        rows = rows + list(dataset.iloc[0:n].index)
        abundantTaxa.append(dataset.loc[rows,'taxonomy'])
    return newtax
                

def findMaxAxVal(maxval, minticks, minstepsize):
    '''Find a maximum value for an axis that is rounded to the nearest whole \
        step for plot ticks
    
    Parameters
    ----------
    maxval : float
        Maximum value of the data.

    Returns
    -------
    maxtickval : float
        Max value for the axis.
    '''
    minstepsize = maxval/minticks
    logval = math.floor(math.log10(minstepsize))
    maxaxval = minticks * math.ceil(minstepsize/10**logval)*10**logval
    return maxaxval


def formatOTUtableData(OTU_table, max_level = 14, taxReassignList = []):
    '''This script reads in and formats an imported raw ASV table by adding \
        taxonomy data.

    Parameters
    ----------
    OTU_table : pandas.DataFrame
        This is the raw imported OTU (or ASV) table with
            index:  OTU IDs
            header: a list of sample names followed by 'taxonomy' at the end
    max_level : int (optional)
        This is the maximum taxonomic level present in the dataset. \
            The default is 14 (i.e., 'D_14__').
    taxReassignList : dict (optional)
        List of taxonomic names in the dataset with the values they should \
            be reassigned. The default is none.

    Returns
    -------
    data : pandas.DataFrame
        Formatted data as a DataFrame.
        New headers include a full taxonomic breakdown
    samples : list
        List of samples in the dataset.

    '''
    # Get sample list
    samples = list(OTU_table.columns)[0:-1]
       
    # Reclassify any values with assignments in the taxReassignList
    if taxReassignList:
        for val in list(taxReassignList):
            OTU_table.loc[OTU_table['taxonomy']==val,
                          'taxonomy'] = taxReassignList[val]
    
    # Format taxonomy list to read better
    print('Formatting taxonomy...')
    for i in np.arange(max_level + 1):
        delstr = 'D_'+str(i)+'__'
        OTU_table['taxonomy'] = OTU_table['taxonomy'].str.replace(delstr, '')        

    # Break taxmap into levels
    taxlist = getUnique(OTU_table['taxonomy'])
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
            OTU_table.loc[OTU_table['taxonomy']==value,
                          'L'+str(i+1)] = splitlist[i]
            
        bar.next()
    bar.finish()

    # Get rid of nans in taxonomy levels
    end_level = len(OTU_table.columns)-len(samples)-1
    cols = levelCols(end_level)
    OTU_table[cols] = OTU_table[cols].replace(np.nan, '')
    
    # Convert values to float
    OTU_table[samples] = OTU_table[samples].astype(float)

    return OTU_table, samples


def getUnique(value_list):
    '''Return alphebetized list of unique values from an original list
    
    Parameters
    ----------
    value_list : list of str
        List to search.

    Returns
    -------
    uniques : list of str
        Each unique value sorted alphabetically.
    '''
    uniques = list(set(value_list))
    uniques.sort()
    return uniques
	


def levelCols(lmax):
    ''''Generates the set of column names for all levels up to a given level

    Parameters
    ----------
    lmax : int
        Maximum taxonomic level for column list.

    Returns
    -------
    cols : list of str
        List of taxonomy column header names corresponding to each level up \
            to the specified max taxonomic level.
    '''
    cols = ['L' + str(i) for i in np.arange(1,lmax+1)]
    return cols

