#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""Rotten Ice Project: Modules

Created on Thu Jun 18 21:49:22 2020
@author: cariefrantz
@project: RottenIce

MODULES SHARED BETWEEN DIFFERENT SCRIPTS IN THE ROTTEN ICE PROJECT
This script contains a collection of functions that are used in multiple
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
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.table import table

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
    * Plots - PCA
    * HTML
    
'''    

#---------------------------------
# SAMPLE NAME [DE]CODING
#---------------------------------

def genSampleList(months, fractions, others = [], replicates = [],
                  locations = ['CS']):
    '''Generates a list of sample names using Rotten Ice Project sample codes;
        list is formatted for use in tables, e.g., where months are one axis
        and possible fractions are another. The list includes samples that
        may not actually have data. For only real samples,
        use genSampleListCS().

    Parameters
    ----------
    months : list of str or 'all'
        Month codes to include. Options are 'M', 'JN', 'JY10', and 'JY11'.
    horizons : list of str or 'all'
        Horizon codes to include, e.g., 'HT', 'HM', 'HB', 'IT', 'SW'.
    others : list of str, optional
        List of additional non-coded samples to include. The default is [].
    replicates : int or str (optional)
        List of replicate codes. The default is [].
    locations : list of str, optional
        Location codes. The default is ['CS'].
    mode : options are 'table' or 'CS_real'. Default is 'table'
        'table' sets up a list for use in table formats where one axis is
            a list of months and the other is a list of possible fractions.
        'CS_real' is the actual list of CS samples.

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


def genSampleListCS():
    '''Generates a list of all sample names collected in the CS dataset'''
    fraction_sets = RottenIceVars.fraction_sets_CS
    month_names = RottenIceVars.months
    samplelist = []
    for month in fraction_sets:
        if 'July' in month:
            samplelist = (samplelist + [month_names[month] + '-' + fraction
                              for fraction in fraction_sets[month]])
        else:
            samplelist = (samplelist
                              + [month_names[month] + '-CS-' + fraction
                                 for fraction in fraction_sets[month]])
    return samplelist
    


def genSampleName(month, fraction, location = 'CS', replicate = ''):
    'Generates sample name from sample metadata'
    if 'JY' in month:
        loccode = month
        if fraction in ['BT', 'BM', 'BB']:
            fraction = 'B'
        elif fraction in ['P1','P2']: fraction = 'Drain'
    else:
        loccode = month + '-' + location
        if fraction in ['PW', 'Drain']: samplename = 'nan'
    if len(replicate)>0:
        samplename = loccode + '-' + fraction + '-' + replicate
    else:
        samplename = loccode + '-' + fraction
    return samplename


def metadataFromSamplename(samplename):
	'''Determines the month and fraction from a sample name'''
	breakout = samplename.split('-')
	month = breakout[0]
	if 'JY' in month:
		fraction = breakout[1]
	else:
		fraction = breakout[2]
	return month, fraction


def samplenameToMonthFractionReplicateGeneTemplate(samplename):
    '''
    Takes a full sequencing sample name like M-CS-HT-1.16S-DNA
    and returns the month (M), fraction (HT), replicate (1), gene (16S) and 
    template (DNA)
    '''
    split1 = samplename.split('.')
    sample = split1[0].split('-')
    month = sample[0]
    fraction = sample[-2]
    replicate = sample[-1]
    gene_template = split1[1]
    split2 = gene_template.split('-')
    gene = split2[0]
    template = split2[1]
    return month, fraction, replicate, gene, template

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
        Options are 'Generic' (default), 'metadata', 'OTU-table',
            and 'alpha-div'
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
    if tabletype in [
            'metadata', 'OTU-table', 'ASV-table-fromBIOM', 'alpha-div']:
        table_fmt = RottenIceVars.data_table_fmts[tabletype]
        header_row = table_fmt['head_row']
        index_col = table_fmt['index_col']
        file_type = table_fmt['filetype']
        
    # Define filetypes
    if file_type == 'csv':
        ftype = [('CSV', '*.csv')]
        sep = ','
    elif file_type == 'tsv':
        ftype = [('TSV', '*.tsv;*.txt')]
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
            index:      taxonomy
            headers:    <samples> + any others
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
    
    # Split out the taxonomy into levels
    tax_split = data.index.to_series().str.split(';', expand=True)
    tax_clean = tax_split.applymap(
        lambda x: x.split('__')[1] if '__' in x else '')
    tax_clean.columns = [f'L{i+1}' for i in range(tax_clean.shape[1])]
    data[tax_clean.columns] = tax_clean
    
    # print('Condensing dataset at L' + str(level) + '...')
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
            
    sample_lists : list of index
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


def formatASVTable(ASV_table, tax_levels, max_level=7):
    '''
    Formats an ASV table by breaking down taxonomy into individual columns
    and normalizing counts to be relative abundance

    Parameters
    ----------
    ASV_table : pandas.DataFrame
        DataFrame containing the ASV table, where index = taxonomic calls for
        the ASVs, columns = sample names.
    tax_levels: list of str
        list naming the different taxonomic levels, to be used as col headers
    max_level : int
        Number of the maximum taxonomic level to analyze. Default = 7 (species)

    Returns
    -------
    ASVTable : pandas.DataFrame
        Formatted ASV table.
    samples : list of str
        List of samples in the inputted ASV_table

    '''
    ASVTable = ASV_table.copy().astype(float).dropna(axis=0)
    
    # Normalize the data (to get relative vs. absolute abundances)
    for column in ASVTable.columns:
        ASVTable[column] = ASVTable[column].values/ASVTable[column].sum()
    
    # Add taxonomy level columns    
    samples = list(ASVTable.columns)
    ASVTable[tax_levels[0:max_level]] = ''
    
    # Format taxonomy list and break into levels
    print('Formatting taxonomy...')
    bar = Bar("", max = len(ASVTable.index))
    for ASV in ASVTable.index:
        splitlist = ASV.split(';')
        splitlist = splitlist[0:max_level]
        for i,s in enumerate(splitlist):
            if s == '__':
                splitlist[i]=''
            else:
                splitlist[i]=s[3:]
            ASVTable.loc[ASV,tax_levels[i]] = '; '.join(splitlist[:i+1])
        bar.next()
    bar.finish()

    
    return ASVTable, samples


def condenseASVTable_by_TaxLevel(formatted_ASV_table,levels,samples):
    '''
    Generates new ASV tables grouped by taxa at different taxonomic levels

    Parameters
    ----------
    formatted_ASV_table : pandas.DataFrame
        ASV table passed through the formatASVTable script.
    levels : list of str
        List of taxonomic levels (column headers in formatted_ASV_table) to
        condense.
    samples : list of str
        List of samples to include in the condensed ASV table.

    Returns
    -------
    tables : TYPE
        DESCRIPTION.

    '''
    tables = {}
    for level in levels:
        print('Condensing dataset at level ' + level + '...')
        # Generate list of unique taxa at the selected level
        taxlist = np.unique(formatted_ASV_table[level])
        
        # Set up a new condensed table
        table_cond = pd.DataFrame(data=None, index=taxlist, columns=samples)
        
        # Condense
        for tax in taxlist:
            table_cond.loc[tax] = formatted_ASV_table[formatted_ASV_table[level]==tax][samples].sum(axis=0)
            
        # Save
        tables[level] = table_cond
        
    return tables



def splitTaxLevels(ASV_table):
    '''
    Takes a formatted ASV table and splits out the taxonomy values,
    generating columns containing each taxonomy call for each taxonomic level

    Parameters
    ----------
    ASV_table : pandas.DataFrame
        This is the raw imported ASV table with
            index:  taxonomy calls in a format that looks like
                d__Eukaryota;p__Chlorophyta;c__Chlorophyceae;
                o__Chlamydomonadales;f__Chlamydomonadales;g__Microglena;s__
            header: a list of sample names

    Returns
    -------
    ASV_taxonomy : pandas.DataFrame
        ASV_table with new columns corresponding to taxonomic levels,
        with the taxonomic call for each level.
    max_level : int
        The maximum taxonomic level found in the dataset.

    '''
    # Convert the index of ASV table (taxonomy calls) to a series
    taxonomy_series = ASV_table.index.to_series()
    # Split out each taxonomic value using ';'
    taxonomy_split = taxonomy_series.str.split(';', expand=True)
    # Remove prefixes
    taxonomy_clean = taxonomy_split.apply(
        lambda col: col.map(lambda x: x.split('__')[1] if isinstance(x, str)
                            and '__' in x else '')
        )
    # Rename columns L1-L7 (or max level)
    taxonomy_clean.columns = [
        f"L{i+1}" for i in range(taxonomy_clean.shape[1])]
    # Determine the max taxonomic level
    max_level = taxonomy_clean.shape[1]
    # Join the columns back to your original DataFrame
    ASV_taxonomy = ASV_table.copy()  # if needed, to avoid modifying in place
    ASV_taxonomy[[f"L{i+1}" for i in range(max_level)]] = taxonomy_clean
    
    return ASV_taxonomy, max_level
    
    
    
    
    
def formatOTUtableData(OTU_table, max_level = 14, tax_reassign_list = []):
    '''
    This script reads in and formats an imported raw ASV table by adding \
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
    OTU_table = OTU_table.copy()
    OTU_table['taxonomy'] = OTU_table.index
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
            
        for L in range(1, min(len(splitlist) + 1, max_level + 1)):
            OTU_table.loc[OTU_table['taxonomy']==value,
                          'L'+str(L)] = splitlist[L-1]
            
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


def groupTaxa(data_table, focusTaxa, max_level):
    '''Merge taxa into new groups based on the focus taxa'''
    data = data_table.copy()
    excludedOthers = []
    
    # Loop through each level
    bar = Bar('', max = max_level) # Progress bar
    for L in np.arange(1, max_level+1):
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
                    if L < max_level:
                        data.loc[indices, 'L'+str(L+1):'L'+str(max_level)] = ''
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
            # Find indices (positions) where taxlist == taxval
            rows_pos = [i for i in range(len(taxlist)) if taxlist[i] == taxval]
            # Get labels (index values) corresponding to those rows
            rows_labels = tax_set.index[rows_pos]
        
            colorlist = tax_set.loc[rows_labels, 'color-L' + str(level+1)]
            
            # Determine the middle color
            if len(colorlist) == 1:
                midclr = colorlist.iloc[0]
            else: midclr = colorlist.iloc[len(colorlist)//2]
            
            # Assign this color to all with this tax value
            tax_set.loc[rows_labels, 'color-L' + str(level)] = midclr
            
    # Return the colormap
    cols_all_clr = ['color-L' + str(i) for i in np.arange(1,max_level+1)]        
    colormap = tax_set[cols_all_clr]
    return colormap   


def genDivergingCmap():
    '''Generates a diverging colormap from default colormap colors'''
    # Pull extreme values from the default colormap
    cmap = cm.get_cmap(RottenIceVars.cmap)
    start = cmap(0)
    end = cmap(0.5)
    #center = viridis((end-start)/2+start)
    center = (1,1,1)
    
    low_vals = []
    high_vals = []
    
    for i in range(3):
        low_vals.append(np.linspace(start[i], center[i], 128))
        high_vals.append(np.linspace(center[i], end[i], 128))
    
    low_vals = np.vstack(np.transpose(low_vals))
    high_vals = np.vstack(np.transpose(high_vals))
    newcolors = np.vstack([low_vals, high_vals])
    
    new_diverging = ListedColormap(newcolors, name = 'Diverging')
    return new_diverging



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
        normdata[sample] = normdata[sample].values/np.sum(
            normdata[sample].values)
    
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


#%%
#---------------------------------
# PLOTS - PCA
#---------------------------------

def biplot(data_table, samples, variables, title,
           color_col, color_map, marker_col, marker_map,
           n_arrows=None, number_arrows = False, marker_size = 10,
           plot_size = 6, legend_pad = 3):
    '''
    Creates a PCA biplot of metadata

    Parameters
    ----------
    data_table : pandas.DataFrame
        Table containing the data to analyze (which must be numeric),
        with samples as rows, variables as columns, including the columns used
        to map color and marker types (can be non-numeric).
    variables : list of str
        List of the variables (column names) to include in the PCA analysis.
        Variable values must be numeric.
    color_col : str
        Column in data_table that determines the marker colors
        (used in color_map).
    color_map : dict of str
        Dictionary mapping values in color_col to colors.
    marker_col : str
        Column in data_table that determines the marker shape
        (used in marker_map).
    marker_map : dict of str
        Dictionary mapping values in marker_col to marker shapes.
    title : str
        Title for labeling the plot
    n_arrows : int, optional
        The maximum number of variable vector arrows to display.
        If defined, the plot will display only the N most important variables.
        The default is None.
    number_arrows : boolean, optional
        Whether or not to number the arrows and add an additional key instead
        of labeling them with the variable names.
        The default is False.
    arrow_scale : float, optional
        Factor used to scale the variable vector arrows. The default is 1.0.
        Increasing this makes the arrows longer/larger.
    marker_size : int, optional
        Size of the sample markers in the PCA plot. The default is 10.
    plot_size : int, optional
        Base height of the plot. It is increased if a key table is added.
        The height determines the width to keep the plot square.
        The default is 6.
    legend_pad : int, optional
        Portion of the plot to reserve for the legend. The default is 2.

    Returns
    -------
    pca_scores : Array of float with dimension [n samples, m pcoa axes]
        PCA-transformed sample data (from pca.transform()).
    var_loadings : Array of float with dimension [n metadata parameters, m pca axes]
        PCA component loadings, i.e., how much each variable influences each
        pca axis (from pca.components).T).
    fig : matplotlib figure
        The biplot figure handle

    '''
    ### PCA CALCULATIONS
    
    # Build matrix with only the indicated metadata values
    valid_rows = [s for s in samples if s in data_table.index]
    invalid_rows = [s for s in samples if s not in data_table.index]
    valid_cols = [v for v in variables if v in data_table.columns]
    invalid_cols = [v for v in variables if v not in data_table.columns]
    for s in invalid_rows:
        print(f"Warning: Sample '{s}' not found in data table index and will be skipped.")
    for v in invalid_cols:
        print(f"Warning: Variable '{v}' not found in data table columns and will be skipped.")
    df = data_table.loc[valid_rows, valid_cols]
    
    # Identify and remove non-numeric columns
    non_numeric_cols = df.select_dtypes(exclude=['number']).columns
    for col in non_numeric_cols:
        print(f"Column '{col}' is non-numeric and will be removed.")
    df_clean = df.drop(columns=non_numeric_cols)
            
    # Delete any samples/rows with missing data
    rows_with_na = df_clean[df_clean.isna().any(axis=1)]
    for idx in rows_with_na.index:
        print(f"Error: '{idx}' contains missing values and will be skipped.")
    df_clean = df_clean.dropna()
    clean_index = df_clean.index
    
    # Subset and scale the data
    X = df_clean.values
    X_scaled = StandardScaler().fit_transform(X)
    
    # Perform 2-dimensional PCA analysis
    pca = PCA(n_components=2)
    pca_scores = pca.fit_transform(X_scaled)
    var_loadings = pca.components_.T  # shape (n_features, 2)
    
    # Align full metadata table to samples used for PCA
    data_for_plot = data_table.loc[clean_index]
    
    
    ### PLOT BIPLOT
    fig_width = plot_size + legend_pad  # add extra width for legend
    fig_height = plot_size              # fix height for square plot
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Add main PCA plot axes with fixed size square area
    ax = fig.add_axes([0.1, 0.1, plot_size/fig_width, plot_size/fig_height])
    ax.set_aspect('equal')
    
    # Extract groupings for coloring and markers
    color_vals = data_for_plot[color_col]
    marker_vals = data_for_plot[marker_col]
    
    # Plot samples
    for c_val in color_map:
        for m_val in marker_map:
            idx = (color_vals == c_val) & (marker_vals == m_val)
            if idx.any():
                ax.scatter(pca_scores[idx, 0], pca_scores[idx, 1],
                           color=color_map[c_val],
                           marker=marker_map[m_val], s=marker_size,
                           label=f"{c_val} / {m_val}",
                           edgecolor='black', alpha=0.8)
    
    # Compute vector magnitudes
    vector_lengths = np.linalg.norm(var_loadings, axis=1)
    
    # Select top N variables if n_arrows is specified
    if n_arrows is not None and n_arrows < len(var_loadings):
        top_indices = np.argsort(vector_lengths)[-n_arrows:]
    else:
        top_indices = np.arange(var_loadings.shape[0])
    
    # Scale arrows based on plot size
    x_range = pca_scores[:, 0].ptp()
    y_range = pca_scores[:, 1].ptp()
    max_plot_radius = 0.3 * max(x_range, y_range)
    max_arrow_length = vector_lengths.max()
    arrow_scale = max_plot_radius / max_arrow_length
    
    arrow_labels = {}
    
    for count, i in enumerate(top_indices, start=1):
        # Label arrows by number or variable name
        if number_arrows:
            arrow_label = str(count)
            arrow_labels[count] = variables[i]
        else:
            arrow_label = variables[i]
    
        # Draw arrows
        ax.arrow(0, 0,
                 var_loadings[i, 0] * arrow_scale,
                 var_loadings[i, 1] * arrow_scale,
                 color='gray', alpha=0.8,
                 head_width=0.03, head_length=0.05)
    
        # Add arrow labels
        ax.text(var_loadings[i, 0] * arrow_scale * 1.15,
                var_loadings[i, 1] * arrow_scale * 1.15,
                arrow_label, color='gray',
                ha='center', va='center', fontsize=9)
    
    # Axis labels with explained variance
    pca_var = np.var(pca_scores, axis=0) / np.sum(np.var(pca_scores, axis=0))
    ax.set_xlabel(f"PC1 ({pca_var[0]*100:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PC2 ({pca_var[1]*100:.1f}%)", fontsize=10)
    
    # Legend outside the plot
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
              fontsize=10, frameon=False)
    
    ax.set_title(title, fontsize=12)
    ax.set_facecolor('white')
    ax.grid(False)
    
    # Show all plot borders
    for side in ['top', 'right', 'bottom', 'left']:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(1.2)
        ax.spines[side].set_color('black')
    
    ax.tick_params(direction='out', length=6, width=1, colors='k')
    
    plt.show()
    
    # Export arrow label key as a table if requested
    if number_arrows is True:
        label_table = pd.DataFrame.from_dict(
            arrow_labels, orient='index', columns=['Variable'])
        label_table.index.name = 'Arrow'
        label_table = label_table.reset_index()
    else:
        label_table = False
    
    
    
    '''
    # This version of the figure code added a table key of arrow values
    # to the plot, but it is not working properly (the table overlaps the plot)

    

    # --- Fixed layout variables ---
    pca_inches = plot_size          # square PCA plot width and height (inches)
    legend_pad_inches = legend_pad  # extra width for the legend area
    table_height_inches = 0         # vertical space reserved for the arrow label table

    # Estimate table height only if arrow labels are requested
    if number_arrows:
        rows_per_unit = 10
        table_height_inches = math.ceil(n_arrows / rows_per_unit) * 1.0 if n_arrows else 1.0

    # Define full figure size
    fig_width = pca_inches + legend_pad_inches
    fig_height = pca_inches + table_height_inches
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Compute Axes placement in figure fractions
    plot_left = 0.1
    plot_bottom = table_height_inches / fig_height
    plot_width_frac = pca_inches / fig_width
    plot_height_frac = pca_inches / fig_height

    # --- Create square PCA Axes ---
    ax = fig.add_axes([plot_left, plot_bottom, plot_width_frac, plot_height_frac])
    ax.set_aspect('equal')

    ### PLOT PCA POINTS

    color_vals = data_for_plot[color_col]
    marker_vals = data_for_plot[marker_col]

    for c_val in color_map:
        for m_val in marker_map:
            idx = (color_vals == c_val) & (marker_vals == m_val)
            if idx.any():
                ax.scatter(pca_scores[idx, 0], pca_scores[idx, 1],
                           color=color_map[c_val],
                           marker=marker_map[m_val], s=marker_size,
                           label=f"{c_val} / {m_val}",
                           edgecolor='black', alpha=0.8)

    # Compute vector magnitudes and select top variables
    vector_lengths = np.linalg.norm(var_loadings, axis=1)
    if n_arrows is not None and n_arrows < len(var_loadings):
        top_indices = np.argsort(vector_lengths)[-n_arrows:]
    else:
        top_indices = np.arange(var_loadings.shape[0])

    # Scale arrows relative to sample spread
    x_range = pca_scores[:, 0].ptp()
    y_range = pca_scores[:, 1].ptp()
    max_plot_radius = 0.3 * max(x_range, y_range)
    max_arrow_length = vector_lengths.max()
    arrow_scale = max_plot_radius / max_arrow_length

    arrow_labels = {}
    for count, i in enumerate(top_indices, start=1):
        arrow_label = str(count) if number_arrows else variables[i]
        if number_arrows:
            arrow_labels[count] = variables[i]
        ax.arrow(0, 0,
                 var_loadings[i, 0] * arrow_scale,
                 var_loadings[i, 1] * arrow_scale,
                 color='gray', alpha=0.8,
                 head_width=0.03, head_length=0.05)
        ax.text(var_loadings[i, 0] * arrow_scale * 1.15,
                var_loadings[i, 1] * arrow_scale * 1.15,
                arrow_label, color='gray',
                ha='center', va='center', fontsize=9)

    # Axis labels, title, and legend
    pca_var = np.var(pca_scores, axis=0) / np.sum(np.var(pca_scores, axis=0))
    ax.set_xlabel(f"PC1 ({pca_var[0]*100:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PC2 ({pca_var[1]*100:.1f}%)", fontsize=10)
    ax.set_title(title, fontsize=12)

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
              fontsize=10, frameon=False)

    ax.set_facecolor('white')
    ax.grid(False)
    for side in ['top', 'right', 'bottom', 'left']:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(1.2)
        ax.spines[side].set_color('black')
    ax.tick_params(direction='out', length=6, width=1, colors='k')

    ### VARIABLE LABEL TABLE (IF NEEDED)

    if number_arrows:
        label_table = pd.DataFrame.from_dict(
            arrow_labels, orient='index', columns=['Variable'])
        label_table.index.name = 'N'
        label_table = label_table.reset_index()

        cell_text = label_table.values.tolist()
        col_labels = ['Arrow', 'Variable']
        col_widths = [0.08, 0.92]

        # Add table in reserved space at bottom
        table_ax = fig.add_axes([0.05, 0.02, 0.9, table_height_inches / fig_height * 0.8])
        table_ax.axis('off')

        table_obj = table_ax.table(
            cellText=cell_text,
            colLabels=col_labels,
            cellLoc='left',
            colWidths=col_widths,
            loc='center'
        )

        table_obj.auto_set_font_size(False)
        table_obj.set_fontsize(8)

        for row in range(len(cell_text) + 1):  # include header
            table_obj[(row, 1)].set_text_props(ha='left')
            table_obj[(row, 1)].PAD = 0

    plt.show()
    '''
    
    return pca_scores, var_loadings, fig, label_table




#%%
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


def genHTMLfile(filename, page_title, subtitle_text, image_filepaths,
                alt_text = '', page_nav_html = ''):
    '''Creates a new HTML file for the fileset from saved plot image(s)

    Parameters
    ----------
    filename : str
        Name and path of the file (includes *.html).
    page_title: str
        Title text for the page.
    subtitle_text: str
        Subtitle text for the page.
    image_filepaths: str or list of str
        Filepath (in the final web directory) that will link to image(s)
    alt_text : str or list of str
        Alternate text for each image.
    page_nav_html : str
        HTML code for the page navigation if page is part of a set of pages.

    Returns
    -------
    None

    '''
    # Get the generic header
    html_head = genHTMLhead(
        page_title, page_nav_html = page_nav_html,
        subtitle_text = subtitle_text)
    # Determine if single image or list of images
    if type(image_filepaths) == str:
        image_filepaths = [image_filepaths]
        alt_text = [alt_text]
    elif type(alt_text) == str:
        alt_text = [alt_text] * len(image_filepaths)
    # Generate html
    html_img = ''
    for i, image in enumerate(image_filepaths):
        html_img = (html_img + '<p><img src = "' + image
                    + '" alt = "' + alt_text[i] + '"></p>')
    # Save file
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
        set_html = genPlotLevelNavHTML(filename_sets[file_set])
        nav_html = nav_html + set_pfx + set_html[crop:] + '  //  '
    nav_html = nav_html[:-6]
    return nav_html
    
    
def genPlotLevelNavHTML(filenames):
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

