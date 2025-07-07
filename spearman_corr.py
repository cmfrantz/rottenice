# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""DRotten Ice Project: SpearmanGrid

Created on Fri July 4 16:09:03 2025
@author: cariefrantz
@project: RottenIce

CALCULATES SPEARMAN CORRELATION FOR METADATA VS. TAXONOMIC GROUPS
This script:
    - Calculates the Spearman correlation coefficients relating selected
        metadata parameters to the abundance of sequences retrieved for the
        most abundant taxonomic groups at different taxonomic levels.
    - Generates heatmap grids for significant correlations between each
        metadata variable / taxonomic group pair.
    - Produces HTML files for each sequence template type with all of the
        heatmaps (at each taxonomic level).

Arguments:  None

Requirements:   
    Metadata table (tsv)
        where rows = samples, columns = metadata characteristics
        header row (meatadata ids)      = 0
        index column (sample names)     = 0
        
    ASV table (csv) for the sequence dataset to analyze
        where rows = ASVs taxonomy, columns = samples
        header row (list of samples)    = 0
        index column (ASV taxa)         = 0
        To prepare, open qiime2 barchart in qiime2 view, set to highest
        taxonomic level, sort by sample index, download csv.
        Then open csv in Excel, copy->transpose to flip rows & columns,
        delete metadata at the bottom, sort data by index, save as new csv.

Example in command line:
    python SpearmanGrid.py

Dependencies and versions used:
    python v. 3.11.7 packaged by Anaconda, Inc.
    fpdf (1.7.2)
    matplotlib (3.8.0)
    numpy (1.26.4)
    pandas (2.1.4)
    progress (1.6)
    scipy (1.11.4)
    tkinter

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
from scipy.stats import spearmanr
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from fpdf import FPDF


####################
# VARIABLES
####################

# Types of sequence data
genes = ['16S-DNA','16S-cDNA','18S-DNA','18S-cDNA',
         '18S-PrimaryProducers-DNA','18S-PrimaryProducers-cDNA']

# Processing settings
max_n_taxa = 20     # Max number of unique taxa to correlate
max_level = 7       # Max taxonomic level to process
tax_levels = ['domain','phylum','class','order','family','genus','species']
significance = 0.01 # P-value cutoff for significance

# Variables to plot
varlist = ['month_num','horizon_num','replicate','lat','lon','depth_in_ice',
           'pH','salinity','SPM','DOC','chl','phaeo','FoFa','PAM','SedLoad',
           'cell_ct','CTC','cells_act','pEPS','POC','nitrogen','CN',
           'temperature','bulk_ice_density','diatom_ct','phyto_ct_select',
           'phyto_ct_all','phyto_ct_other','gels_total','gels_sm','gels_md',
           'gels_lg','gels_xl'
]
           

# Colormap to use
# See https://matplotlib.org/tutorials/colors/colormaps.html for a full list.
color_set = 'viridis'


####################
# FUNCTIONS
####################
def genDivergingCmap(color_set):
    '''
    Generates a diverging colormap from default colormap colors

    Parameters
    ----------
    color_set : string
        Name of the colormap to be used.
        See https://matplotlib.org/tutorials/colors/colormaps.html

    Returns
    -------
    new_diverging : TYPE
        DESCRIPTION.

    '''
    # Pull extreme values from the default colormap
    cmap = matplotlib.colormaps[color_set]
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

def fileGet(title, directory = os.getcwd(),
            file_type = 'tsv', header_row = 0, index_col = 0):
    '''
    Imports data from a file chosen with user input GUI.

    Parameters
    ----------
    title : str
        Text to put in the user input window.
    directory : sstr, optional
        The start directory to search in. The default is os.getcwd().
    file_type : str, optional
        The type of file. Options are currently:
            'csv'               Comma-seperated file
            'tsv' (default)     Tab-separated file.
    header_row : int, optional
        The row of the file to capture for the data headers. The default is 0.
    index_col : int, optional
        The column of the file to capture for the data index. The default is 0.

    Returns
    -------
    filename : str
        Selected filepath.
    dirPath : str
        Directory path for the selected file.
    data : pandas.DataFrame
        DataFrame containing the indexed data read in from the file.

    '''
    
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


def buildCorrTable(metadata, abundance_table, varlist, samples,
                   significance = 0.05):
    '''
    Calculates Spearmann correlation coefficients for metadata vs. relative
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
    
    #Build empty table
    corrtable = pd.DataFrame(
        data = None, index = abundance_table.index, columns = varlist)
    
    # Loop through each taxon & calculate
    for tax in abundance_table.index:
        # Loop through each metadata variable & calculate
        for var in varlist:
            # Grab the lists of values to correlate
            varvals = metadata.loc[samples, var].to_numpy(dtype=float)
            counts = abundance_table.loc[tax, samples].to_numpy(dtype=float)
            
            #Calculate Spearman coefficient
            corr, p = spearmanr(varvals, counts, nan_policy='omit')
            
            # If p value <= 0.05, add it to the heatmap matrix
            if p <= significance:
                corrtable.loc[tax,var] = corr
                
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
    #
    #IMPORT FILES
    
    # Import metadata table
    filename, directory, metadata = fileGet(
        'Select the metadata table', file_type='tsv')
    metadata = metadata.dropna(how = 'all')
    metadata= metadata.replace('na', np.nan)
    # Add this in if paring down the variable list
    # metadata = metadata[varlist]
    
    # Import ASV tables
    asv_tables = {}
    for gene in genes:
            filename, directory, data = fileGet(
                'Select ' + gene + ' ASV table', directory=directory,
                file_type = 'csv')
            asv_tables[gene] = data
    
    #
    # BUILD METADATA VS TAXONOMIC ABUNDANCE CORRELATION HEATMAPS
    
    # Generate colormap to use
    cmap = genDivergingCmap(color_set)
    
    # Prepare container for tables
    spearman_tables = {}

    # Loop through each sequencing dataset
    for gene in genes:
        
        print('*******************************************************\n'
              'Determining Spearman correlations for ' + gene + ' dataset\n'
              '*******************************************************')
        # Format & condense data by level
        ASV_table, samples = formatASVTable(
            asv_tables[gene], tax_levels, max_level)
        level_tables = condenseASVTable_by_TaxLevel(
            ASV_table, tax_levels[0:max_level], samples)
        
       
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
            if ds.shape[0] >= max_n_taxa:
                ds = ds.iloc[0:max_n_taxa]
                
            # Calculate Spearman correlation and return significant values
            print('  Calculating Spearman correlation coefficients for taxa')
            corrtable = buildCorrTable(
                metadata, ds, varlist, samples, significance = significance)
            
            # Save correlation table
            spearman_tables[gene+'_L'+str(l+1)] = corrtable
            
            # Build heatmap table to display results of Spearman correlations
            fig, ax = plt.subplots(figsize=(30,10))
            im = buildHeatmap(
                ax, cmap, corrtable,
                ('Spearman correlation of ' + gene +
                 ' data at taxonomic level ' + str(l+1) + ' ' + level))
            texts = annotateHeatmap(im, fontsize = 6)
            fig.tight_layout()
            plt.show()
        
            # Save plots
            print('  Saving ' + gene + ' L' + str(l+1) + 'figure...')
            img_filename = 'SpearmanCorrFig_' + gene + '_L' + str(l+1)
            fig.savefig(
                directory + '\\' + img_filename + '.svg', transparent = True)
            fig.savefig(
                directory + '\\' + img_filename + '.pdf', format='pdf')   
            
              

    #
    # DO METADATA CORRELATION ANALYSIS
    
    # CalculateSpearman correlations and return significant values
    print('  Calculating Spearman correlation coefficients for metadata')
    corrtable = pd.DataFrame(data = None, index = varlist, columns = varlist)
    for var1 in varlist:        
        # Loop through each variable & calculate
            for var2 in varlist:
                # Grab the lists of values to correlate
                var1vals = metadata.loc[samples, var1].to_numpy(dtype=float)
                var2vals = metadata.loc[samples, var2].to_numpy(dtype=float)
                
                #Calculate Spearman coefficient
                corr, p = spearmanr(var1vals, var2vals, nan_policy='omit')
                
                # If p value <= 0.05, add it to the heatmap matrix
                if p <= significance:
                    corrtable.loc[var1,var2] = corr
    corrtable = corrtable.astype(float)
    
    # Save correlation table
    spearman_tables['metab_metadata'] = corrtable
                    
    # Build heatmap table to display results of Spearman correlations
    fig, ax = plt.subplots(figsize=(30,30))
    im = buildHeatmap(
        ax, cmap, corrtable, 'Spearman correlation of metadata')
    texts = annotateHeatmap(im, fontsize = 6)
    fig.tight_layout()
    plt.show()
    # Save plots
    print('  Saving metadata figure...')
    fig.savefig(
        directory + '\\' + 'SpearmanCorrFig_Metadata.svg', transparent = True)
    fig.savefig(
        directory + '\\' + 'SpearmanCorrFig_Metadata.pdf', format = 'pdf')
    
    
    # Build analysis document for each variable
    print ('  Saving metadata analysis pdf')
    pdf = FPDF()
    for var in varlist:
        # Get and sort the correlating variable list
        coeffs = corrtable.loc[var].dropna()[1:]
        coeffs.sort_values(ascending = False, inplace = True)
        pos = coeffs[coeffs>0]
        neg = coeffs[coeffs<0].sort_values()
        
        # Create a new page & enter results
        pdf.add_page()
        pdf.set_font('Arial', size = 12, style = 'B')
        pdf.cell(200, 10, txt = var, ln=1, align = 'C')
        pdf.set_font('Arial', size = 12, style = '')
        pdf.cell(
            0, 20, txt = (
                'Spearman correlation coefficients for correlations with ' +
                'other metadata parameters where p <=' + str(significance)),
            ln=1, align = 'C')
        pdf.set_font('Arial', size = 12, style = 'U')
        pdf.cell(0, 20, txt = 'Positive correlations:', ln=1, align = 'L')
        for i,p in enumerate(pos):
            pdf.set_font('Arial', size = 12, style = '')
            pdf.cell(
                0, 10,
                txt = str(round(pos.iloc[i],2)) + '  ' + pos.index[i],
                ln=1, align = 'L')
        pdf.set_font('Arial', size = 12, style = 'U')
        pdf.cell(0,20, txt = 'Negative correlations:', ln=1, align = 'L')
        for i,n in enumerate(neg):
            pdf.set_font('Arial', size = 12, style = '')
            pdf.cell(
                0, 10,
                txt = str(round(neg.iloc[i],2)) + '  ' + neg.index[i],
                ln=1, align = 'L')
    pdf.output(directory + '\\' + 'SpearmanCorr_Metadata.pdf')
    
    # Save and export tables as Excel workbook
    with pd.ExcelWriter(directory + '\\SpearmanCorrTable.xlsx') as writer:
        for n, df in enumerate(list(spearman_tables)):
            spearman_tables[df].to_excel(writer, sheet_name = df)
        