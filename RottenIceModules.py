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
    pip install bokeh


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
from matplotlib import cm
import math

from bokeh.io import export_svgs
from bokeh.models import Legend, LegendItem, Div
from bokeh.plotting import figure

import RottenIceVars




####################
# FUNCTIONS
####################
'''
Functions are listed alphabetically within categories:
    
    * Sample name [de]coding
    * Data Prep - General
    * Data Prep - OTU Tables
    * Colors
    * Plots - General
    * Plots - Bokeh
    * HTML
    
'''    

#---------------------------------
# SAMPLE NAME [DE]CODING
#---------------------------------

def genSampleList(months, fractions, others = [], replicates = [],
                  locations = ['CS']):
    '''Generates a list of sample names using Rotten Ice Project sample codes

    Parameters
    ----------
    months : list of str
        Month codes to include. Options are 'M', 'JN', 'JY10', and 'JY11'.
    horizons : list of str
        Horizon codes to include, e.g., 'HT', 'HM', 'HB', 'IT', 'SW'.
    others : list of str, optional
        List of additional non-coded samples to include. The default is [].
    replicates : int or str (optional)
        List of replicate codes. The default is [].
    locations : list of str, optional
        Location codes. The default is ['CS'].

    Returns
    -------
    samplelist : list of str
        List of sample names.
    '''
    samplelist = []
    # Loop through locations, months, fractions, and replicates
    for location in locations:
        for month in months:
            # Format month code
            if 'JY' in month:
                monthcode = month
            else:
                monthcode = month + '-' + location
            # Format horizon code
            for fraction in fractions:
                # Add replicates if sample names include them
                if replicates:
                    for replicate in replicates:
                        samplelist.append(monthcode + '-'
                                          + fraction + '-' + str(replicate))
                else:
                    samplelist.append(monthcode + '-' + fraction)
    # Add others to the end
    if others:
        samplelist = samplelist + others
    return samplelist


def genSampleName(month, fraction, location = ['CS'], replicate = []):
    'Generates sample name from sample metadata'
    if 'JY' in month:
        loccode = month
    else:
        loccode = month + '-' + location
    if replicate:
        samplename = loccode + '-' + fraction + '-' + replicate
    else:
        samplename = loccode + '-' + fraction
    return samplename


def metadataFromSamplename(samplename):
	'''Determines the month and horizon from a sample name'''
	breakout = samplename.split('-')
	month = breakout[0]
	if 'JY' in month:
		fraction = breakout[1]
	else:
		fraction = breakout[2]
	return month, fraction



#---------------------------------
# DATA PREP: GENERAL
#---------------------------------

def fileGet(title, tabletype = 'Generic', directory = os.getcwd(),
            file_type = 'csv', header_row = 0, index_col = 0):
    '''This script imports data from a file chosen with a user input GUI.
    
    Parameters
    ----------
    title : str
        Text to put in the user input window.
    tabletype: str
        Type of table to input.
        Options are 'Generic' (default), 'metadata', and 'OTU-table'.
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
    
    # If type of file is specified, use the proper format
    if tabletype == 'metadata':
        header_row = RottenIceVars.metadata_head_row
        index_col = RottenIceVars.metadata_index_col
        file_type = RottenIceVars.metadata_filetype
    elif tabletype == 'OTU-table':
        header_row = RottenIceVars.OTU_table_head_row
        index_col = RottenIceVars.metadata_index_col
        file_type = RottenIceVars.metadata_filetype
    
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



#---------------------------------
# DATA PREP: OTU TABLES
#---------------------------------

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
    taxnames = genTaxName(data[cols])
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


def findAbundantTaxa(datasets, sample_lists, fcutoff = 0.1, max_taxa = 100):
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
            
    samples : list of index
        List containing lists of sample names / indices (as str or index) \
            to be included in the analysis of each passed dataset.
        All sample names should also be headers in all included DataFrames.
            
    fcutoff : float (optional)
        If any taxa is greater than this value in abundance in any sample, \
            it will be included. Default is 0.1 (10%).
            
    max_taxa : int (optional)
        The maximum number of abundant taxa to return.

    Returns
    -------
    abundantTaxa : list of str
        Alphabetically-sorted list of the most abundant taxa in the selected \
            samples.
    '''
    abundantTaxa = []
    # Search in each dataset
    for dataset, samples in zip(datasets, sample_lists):
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
        n = max_taxa - len(rows)
        dataset['totals'] = np.sum(dataset[samples], axis = 1)
        dataset = dataset.sort_values('totals', ascending = False)
        rows = rows + list(dataset.iloc[0:n].index)
        abundantTaxa.append(dataset.loc[rows,'taxonomy'])
        
    # Get unique values in the list
    abundantTaxa = [taxon for taxset in abundantTaxa for taxon in taxset]
    abundantTaxa = getUnique(abundantTaxa)
    
    return abundantTaxa



def formatOTUtableData(OTU_table, max_level = 14, tax_reassign_list = []):
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
    tax_reassign_list : dict (optional)
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
       
    # Reclassify any values with assignments in the tax_reassign_list
    if tax_reassign_list:
        for val in list(tax_reassign_list):
            OTU_table.loc[
                OTU_table['taxonomy']==val,
                'taxonomy'] = tax_reassign_list[val]
    
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
        elif '; ' in value:
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


def genFilenamesByLevel(filename_prefix, max_level, filename_suffix = '',
                       cmap = '', filetype = ''):
    '''Generates a set of filenames up to a given taxonomic level'''
    # Generate the filename suffix
    if len(filename_suffix) > 0:
        filename_suffix = '_' + filename_suffix
    if len(cmap) > 0:
        filename_suffix = cmap + filename_suffix
    if len(filetype) > 0:
        filename_suffix = filename_suffix + '.' + filetype
    # Generate list of filenames
    filenames = []
    for level in range(1, max_level + 1):
        filenames.append(filename_prefix + '_L' + str(level)
                         + '_' + filename_suffix)
    return filenames


def genTaxName(taxcols):
    '''Generates combined taxonomy names from individual levels'''
    taxname = ['>'.join(list(taxcols.loc[i])) for i in list(taxcols.index)]
    return taxname


def groupTaxa(data_table, focusTaxa):
    '''Merge taxa into new groups based on the focus taxa and "others"'''
    data = data_table.copy()
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
        data['taxonomy'] = genTaxName(data[levelCols(L)]) 
        # Generate list of unique taxonomic names at this level                 
        taxlist = getUnique(data.loc[indices]['taxonomy'])
        # Remove any empty values
        taxlist = list(filter(None, taxlist))                               
        # If the list isn't empty...
        if taxlist:
            # Loop through each unique taxonomy name
            for value in taxlist:
                # Get the full taxonomic name for everything in the category
                indices = [i for i in data.index
                           if data.loc[i]['taxonomy']==value]
                taxvals = np.unique(list(data_table.loc[indices, 'taxonomy']))
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



#---------------------------------
# COLORS
#---------------------------------

def colorsFromCmap(n, cmap):
    '''Generate the set of colors to use from the chosen colormap

    Parameters
    ----------
    n : int
        Number of colors needed (number of items to color).
    cmap : str
       The matplotlib colormap name to use for picking colors.

    Returns
    -------
    colors : list of str
        List of hex color codes for each item.
    '''
    colormap = cm.get_cmap(cmap)
    colorset = colormap(np.linspace(0,1,n))
    colors=[]
    for color in colorset: colors.append(cm.colors.rgb2hex(color))
    return colors         

    
def colorMap2Tax(tax_set, max_level, cmap = RottenIceVars.cmap):
    '''This function maps colors to taxonomy values for plots. \
        Doing this keeps colors consistent across levels for taxonomic groups.   

    Parameters
    ----------
    tax_set : pandas.DataFrame
        Full set of taxonomic data including taxonomy breakdown matrix \
            with levels broken out.
    cmap : str
       The matplotlib colormap name to use for picking colors.
    max_level : int
        The maximum taxonomic level to consider.

    Returns
    -------
    colormap : list of tuples
        List of colors (as tuples) to use for each taxon.
    '''
   
    cols_all = levelCols(max_level)
    
    colnames = ['color-' + col for col in cols_all]
    for c in colnames:
        tax_set[c] = ''
    
    # Find the greatest level with <= 256 unique values
    print('Finding level to use as colormap index...')
    uniques=[]
    for level in np.arange(1, max_level+1):
        cols_sub = levelCols(level)
        tax_names = genTaxName(tax_set[cols_sub])
        uniquetax = getUnique(tax_names)
        uniques.append(uniquetax)
        if len(uniquetax) > 256: break
    cmaplevel = level-1
    
    # Generate colors for each OTU at the highest taxonomic level
    # Generate tax names for each group at the level
    cols_clevel = levelCols(cmaplevel)
    uniquetax = getUnique(genTaxName(tax_set[cols_clevel]))
    print('Using Level ' + str(cmaplevel) + ' for colormap index ('
          + str(len(uniquetax)) + ' unique values)')
    uniquetax.sort()
    colorlist = colorsFromCmap(len(uniquetax), cmap)
    colordict = dict(zip(uniquetax, colorlist))
    
    # Fill in the values for the max level and all lower levels
    tax_set.loc[:, 'taxonomy'] = genTaxName(tax_set[cols_clevel])
    cols_low = np.setdiff1d(levelCols(max_level), levelCols(cmaplevel-1))
    print('Assigning colors to L' + str(cmaplevel) + '-L'
          + str(max_level) + '...')
    cols_low_clr = ['color-' + col for col in cols_low]
    for i in list(tax_set.index):
        tax_set.loc[i, cols_low_clr] = colordict[tax_set.loc[i, 'taxonomy']]
    
    # Fill in the values for the higher levels
    print('Assigning colors to L1-L' + str(cmaplevel-1))
    # Loop through remaining levels.
    # At each level, pick the mid value for each set.
    for level in range(cmaplevel-1, 0, -1):
        cols_sub = levelCols(level)
        taxlist = genTaxName(tax_set[cols_sub])
        # Loop through each unique value
        for taxval in uniques[level-1]:
            # Find all color values in the next level belonging to this
            # taxonomic set
            rows = [i for i in np.arange(len(taxlist))
                    if taxlist[i] == taxval]
            colorlist = tax_set.iloc[rows]['color-L' + str(level+1)]
            # Determine the middle color
            if len(colorlist) == 1:
                midclr = colorlist[0]
            else: midclr = colorlist[int(len(colorlist)/2)]
            # Assign this color to all with this tax value
            tax_set.loc[colorlist.index, 'color-L' + str(level)] = midclr
            
    # Return the colormap
    cols_all_clr = ['color-L' + str(i) for i in np.arange(1,max_level+1)]        
    colormap = tax_set[cols_all_clr]
    return colormap   



#---------------------------------
# PLOTS - GENERAL
#---------------------------------

def findMaxAxVal(max_val, min_ticks = 2):
    '''Find a maximum value for an axis that is rounded to the nearest whole \
        step for plot ticks
    
    Parameters
    ----------
    max_val : float
        Maximum value of the data.
    min_ticks : int, optional
        Minimum number of ticks. Default is 2.

    Returns
    -------
    max_ax_val : float
        Max value for the axis.
    '''
    min_step_size = max_val/min_ticks
    log_val = math.floor(math.log10(min_step_size))
    max_ax_val = min_ticks * math.ceil(min_step_size/10**log_val)*10**log_val
    return max_ax_val


def genMarkerKwargs(samplename):
    '''Generates the kwargs for formatting scatterplot markers'''
    month, fraction = metadataFromSamplename(samplename)
    kwargs = {
        'color'     : RottenIceVars.plotColorsByMonth[month],
        'marker'    : RottenIceVars.plotMarkersByFraction[fraction],
        }
    if RottenIceVars.plotMarkerBorderByMonth[month] == True:
        kwargs.update(RottenIceVars.plotMarkerLineProperties)
    return kwargs


def genSamples_w_Markers(months, fractions):
    '''Generates a list of samples and list of marker properties for each \
        sample from a list of months and fractions being plotted.'''
    samplelist = genSampleList(months, fractions)
    markers = []
    for sample in samplelist:
        kwargs = genMarkerKwargs(sample)
        kwargs.update({'alpha': 0.8})
        markers.append(kwargs)
    return samplelist, markers


def plotMarkerPropertiesFromSamplelist(
        samplelist, monthparams = RottenIceVars.plotColorsByMonth,
        fractionparams = RottenIceVars.plotMarkersByFraction):
	'''Determines the plot marker and marker color from the sample name'''
	colors = []
	markers = []
	for sample in samplelist:
		month, fraction = metadataFromSamplename(sample)
		colors.append(RottenIceVars.plotColorsByMonth[month])
		markers.append(RottenIceVars.plotMarkersByFraction[fraction])
	return colors, markers



#---------------------------------
# PLOTS - BOKEH
#---------------------------------

def buildBarPlot(plt_data, colorlist, y_max, save_image = False,
                 filename='bar', invert = False,
                 show_x_ax = True, x_range = [],
                 legend = False, **kwargs):
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
    y_max : float or int
        Y axis maximum value. Use 1 for relative abundance plots.
    save_image: boolean, optional
        If true, a png snapshot of the bar plot is saved.
        If false, no snapshot is saved.
    filename : str, optional
        String filepath + filename prefix to use when saving barplot png image
    invert : Boolean, optional
        If invert is true, the plot will be mirrored (upside-down). \
            The default is False.
    show_x_ax : Boolean, optional
        Whether or not to display the x axis labels. The default is True.
    x_range : bokeh.models.FactorRange, optional
        Bokeh range object that manually sets the range label.
        The default is [], which assumes that the x range should be the same \
            as the data columns (samples).
    legend: Boolean, optional
        Whether or not to display a figure legend. The default is False.
    kwargs : dict, optional
        Dictionary containing optional figure formatting options for \
            bokeh.figure.

    Returns
    -------
    p   : bokeh figure
        Bokeh figure object containing the plot.

    '''
    # Update plot properties kwargs
    title = 'bar'
    if x_range:
        plot_formatting_args = {'x_range' : x_range}
    else:
        plot_formatting_args = {'x_range' : plt_data['samples']}
    if kwargs:
        if 'title' in list(kwargs):
            title = kwargs['title']
        kwargs.update(plot_formatting_args)
    else:
        kwargs = plot_formatting_args
    
    # Create figure

    # print('Building ' + title + ' plot...')
    p = figure(**kwargs)
    
    # Add bars corresponding to each taxon
    taxlist = list(plt_data)
    taxlist.remove('samples')
    vbars = p.vbar_stack(
        taxlist, x = 'samples', width = 0.9, color = colorlist,
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
        legenditems = []
        # Create a legend item for each taxon
        for i in np.arange(len(taxlist)):
            legenditems.append(LegendItem(label = taxlist[i],
                                          renderers = [vbars[i]]))
        # Build the legend from the legend items    
        legend = Legend(items = legenditems, location = 'top_right', 
                        glyph_height = 15, glyph_width = 15,
                        label_text_font_size = '9px',)
        p.add_layout(legend, 'right')
        
    # Save figure
    if save_image:
        p.output_backend = 'svg'
        export_svgs(p, filename = 'barplot_' + title + '.svg')
        
    return p


def buildBokehNavDiv(page_title, subtitle_text, local_nav_html):
    '''Builds header div for bokeh web plots

    Parameters
    ----------
    page_title : str
        Title to use for the page header
    subtitle_text : str
        Text (including any html formatting) to be used for the subtitle \
            paragraph. It will be wrapped in <p>subtitle_text</p>.
    local_nav_html: str
        HTML to be used for navigating between same-type plots \
            (e.g. pages of plots at different taxonomic levels). \
                It will be wrapped in <p>local_nav_html</p>.

    Returns
    -------
    header_div : bokeh.models.Div
        Div class containing the header html for a bokeh page.
    '''
    html_nav_head = genHTMLhead(
        page_title, set_nav_html = local_nav_html,
        subtitle_text = subtitle_text)
    
    header_div = Div(text = html_nav_head)
    return header_div
     
    
def buildPlotDicts(data, level, samples, colormap):
    '''Processes raw data to produce two datasets and a mapped color list:
        - raw counts (absolute data)
        - relative abundances (relative data)

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame of OTU table. Should contain each sample as a column as \
            well as broken out taxa columns. Each row is a unique OTU ID.
    level : int
        Taxonomic level for the plot
    samples: list of str
        Samples to include in the plot
    colormap: pandas.DataFrame
        Dataframe mapping colors to taxonomic calls generated by \
            colorMap2Tax. Indices of colormap and data should match.

    Returns
    -------
    abspltdict : dict
        Dictionary containing each taxonomic ID as keys and lists of counts \
            for each sample as values. Also contains the key 'sample' with \
                a list of sample names.
    relpltdict : dict
        Dictionary containing each taxonomic ID as keys and lists of \
            normalized counts for each sample as values. Also contains the \
                key 'sample' with a list of sample names.
    colorlist : list of str
        List of color values corresponding to each dictionary key.

    '''
        
    # Map taxonomy names to colormap
    colors = colormap.loc[data.index, 'color-L' + str(level)]
    cols = levelCols(level)
    taxnames = genTaxName(data[cols])
    colors.index = taxnames
    
    # Condense dataset to the selected level
    data = condenseDataset(data, level, samples)
    data = data.fillna(0)
    # Create normalized data
    normdata = data.copy()
    for sample in samples:
        normdata[sample] = normdata[sample].values/np.sum(normdata[sample])
    
    # Make color list
    colorlist = []
    for taxgroup in data.index:
        if isinstance(colors.loc[taxgroup], str):
            colorlist = colorlist + [colors.loc[taxgroup]]
        else:
            colorlist = colorlist + [list(colors.loc[taxgroup])[0]]
    
    # Make dictionaries
    abspltdict = df2Dict(data)
    relpltdict = df2Dict(normdata)
    return abspltdict, relpltdict, colorlist


def df2Dict(data):
    '''Stitches together a data dictionary for bokeh plots'''
    pltdata = {'samples' : data.columns.values}
    for value in data.index:
        pltdata[value] = data.loc[value].values
    return pltdata


def genLegendOutside(taxlist, colors, key_val_ht = 24, key_top_pad = 40):
    '''Draws a legend to add to the figure

    Parameters
    ----------
    taxlist : list of str
        List of taxa to put in the legend.
    colors : list of str
        List of colors corresponding to each taxon in taxlist.
    keyvalht : int, optional
        Size of the key color box in px. The default is 10.
    keytoppad : int, optional
        Padding at the top of the key box in px. The default is 40.

    Returns
    -------
    l : TYPE
        DESCRIPTION.

    '''
    # Generate dummy plot to hold legend
    l = figure(height = len(taxlist) * key_val_ht + key_top_pad)
    items = []
    for i in np.arange(len(taxlist)):
        items.append(l.square(1,1, size = 15, color = colors[i],
                               legend_label = taxlist[i]))
    
    # Add legend
    l.legend.location = 'top_left'
    l.legend.title = 'Key'
    l.legend.label_text_font_size = '9px'
    l.legend.spacing = 0
    l.legend.margin = 0
    
    # Make everything else invisible
    for item in items: item.visible=False
    l.xaxis.visible = False
    l.yaxis.visible = False
    l.xgrid.visible = False
    l.ygrid.visible = False
    l.outline_line_width = 0
        
    return l



#---------------------------------
# HTML FILES
#---------------------------------

def addHTMLhead(filename, page_title, page_nav_html = '', subtitle_text = ''):
    '''Adds conserved HTML header code to an existing HTML file'''
    # Get the header
    html_head = genHTMLhead(
        page_title,
        page_nav_html = page_nav_html,
        subtitle_text = subtitle_text)
    
    # Open the existing file
    with open(filename) as file:
        file_text = file.read()
        
    # Insert the header at the appropriate point
    found = False
    insert_points = ['<body>', '<html>']
    for tag in insert_points:
        if tag in file_text:
            file_text = file_text.replace(tag, tag + html_head)
            found = True
            break
    if not found:
        file_text = html_head + file_text
    
    # Overwrite the file
    with open(filename, 'w') as file:
        file.write(file_text)
    print('Added header to ' + filename)



def genHTMLhead(page_title, page_nav_html='', subtitle_text=''):
    '''Generates HTML header code for a page'''
    # Generate conserved header from RottenIceVars
    html_head = RottenIceVars.web_nav_header + '<p>'
    file_sets = RottenIceVars.file_sets
    for fset in file_sets:
        html_head = (html_head + '<b>'
                     + file_sets[fset]['title']
                     + ':  </b><a href="'
                     + RottenIceVars.web_dir
                     + file_sets[fset]['land_page'] + '">Link</a><br />')
    html_head = html_head + RottenIceVars.other_nav + '</p>'
    
    # Add page header
    html_page_head = '<h1>' + page_title + '</h1>'
    if len(page_nav_html)>0: page_nav_html = '<p>' + page_nav_html + '</p>'
    if len(subtitle_text)>0: subtitle_text = '<p>' + subtitle_text + '</p>'
    
    return html_head + html_page_head + page_nav_html + subtitle_text


def genHTMLfile(filename, page_title, subtitle_text, image_filepath,
                page_nav_html = '', alt_text = ''):
    '''Creates a new HTML file for the fileset from a saved plot image'''
    html_head = genHTMLhead(
        page_title, page_nav_html = page_nav_html,
        subtitle_text = subtitle_text)
    html_img = ('<p><img src = "' + image_filepath
                + '" alt = "' + alt_text + '"></p>')
    with open(filename, 'a') as file:
        file.write('<HTML>' + html_head + html_img + '</HTML>')
    print('Saved ' + filename)


def genPlotNavHTML_Sets(filename_sets):
    '''Generates plot navigation HTML from a set of filename prefix lists. \
        Each list has filenames corresponding to each level for which plots \
            were made. If there is only one set of files, \
                (e.g., only algae data or only 16S data), \
                    use genPlotNavHTML instead.'''
    nav_html = RottenIceVars.nav_html_start
    crop = len(nav_html)
    for file_set in filename_sets:
        set_pfx = '  <em>' + file_set + ': </em>'
        set_html = genPlotNavHTML(filename_sets[file_set])
        nav_html = nav_html + set_pfx + set_html[crop:] + '  //  '
    nav_html = nav_html[:-6]
    return nav_html
    
    
def genPlotNavHTML(filenames):
    '''Generates plot navigation HTML from a list of filename prefixes \
        corresponding to each level for which plots were made. \
            If there are multiple sets of files (e.g., 16S & 18S data), \
                use genPlotNavHTML_Sets instead. '''
    nav_html = RottenIceVars.nav_html_start
    for i, filename in enumerate(filenames):
        nav_html = (nav_html + ' <a href = "' + filename + '.html">'
                    + 'L' + str(i+1) + '</a>  /')
    nav_html = nav_html[:-3]
    return nav_html

