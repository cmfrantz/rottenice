# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 12:36:02 2025

@author: cariefrantz
@project: RottenIce

GENERATES 2D PCOA PLOT COLLECTION

This script generates a collection of PCoA plots split out by month (all
fractions) and horizon (all months) for 16S and 18S community data.

This script was created as part of the Rotten Ice Project

Arguments:  None

Requirements:   
    Metadata table (tsv)
        where rows = samples, columns = metadata characteristics
        unlike the other scrips in this RottenIce collection, the metadata
        file used here is the raw tsv used in QIIME2
        
    PCoA ordinate tables (tsv) for the 16S and 18S data
        The PCoA ordinate table should be exported from QIIME2:
            qiime tools export \
            # --input-path pcoa.qza \
            # --output-path directory/export

Example in command line:
    python pcoa.py

Dependencies Install:
    sudo apt-get install python3-pip python3-dev
    pip install tkinter
    pip install pandas
    pip install matplotlib


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

import os
from tkinter import *
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import RottenIceVars


####################
# VARIABLES

monthmap = {
    'M-CS'  : 'May',
    'JN-CS' : 'June',
    'JY10'  : 'July - dirty floe',
    'JY11'  : 'July - clean floe'
    }

materialmap = {
    'I'     : 'ice-only',
    'H'     : 'whole horizon',
    'B'     : 'brine',
    'P'     : 'percolate',
    'D'     : 'drain',
    'PW'    : 'pondwater',
    'SW'    : 'seawater',
    'BW'    : 'below-ice water'
}

marker_map_horizon = {
    'top'       : RottenIceVars.plotMarkersByHorizon['T'],
    'middle'    : RottenIceVars.plotMarkersByHorizon['M'],
    'bottom'    : RottenIceVars.plotMarkersByHorizon['B'],
    'whole'     : RottenIceVars.plotMarkersByHorizon['A'],
    'water'     : RottenIceVars.plotMarkersByHorizon['W']
    }

''' Old color map
color_map_fraction = {
    'I'         : '#2066a8', # blue
    'H'         : '#f6d6c2', # tan
    'B'         : '#d47264', # orange-red
    'P'         : '#B05F53', # red
    'D'         : '#874940', # dark red
    'PW'        : '#b5d1ae', # light green
    'SW'        : '#326b77', # dark green
    'BW'        : '#A799B7' # lavender
    }
'''
# Pull the standard color map for fractions
color_map_fraction = {
    'I'         : RottenIceVars.plotColorsByMaterial['I'],
    'H'         : RottenIceVars.plotColorsByMaterial['H'],
    'B'         : RottenIceVars.plotColorsByMaterial['B'],
    'P'         : RottenIceVars.plotColorsByMaterial['P'],
    'D'         : RottenIceVars.plotColorsByMaterial['I'],
    'PW'        : RottenIceVars.plotColorsByMaterial['PW'],
    'SW'        : RottenIceVars.plotColorsByMaterial['SW'],
    'BW'        : RottenIceVars.plotColorsByMaterial['BW']
    }

# Pulls hex values for standard colors for each month
color_map_month = RottenIceVars.plotColorsByMonth

size_map_template = {
    'cDNA'      : 80,
    'DNA'       : 30
    }


####################
# FUNCTIONS
####################

def fileGet(title, directory = os.getcwd(),
            file_type = 'tsv', header_row = 0, skiprows = None, index_col = 0):
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
        If no header, header_row = None
    skiprows : int, optional
        Number of lines at the start of the file to skip. The default is None.
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
                       header = header_row, skiprows = skiprows, 
                       index_col = index_col)
    
    return filename, dirPath, data  


def prepMetadataWpcoa(metadata, title = '', directory = os.getcwd()):
    '''
    Select the pcoa file and add to a metadata table. Opens a user interface
    to select files.
    
    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata table with sample IDs as the index, metadata parameters as
        columns.
    title : str, optional
        Adds a title to the user input box to indicate which file to grab.
        The default is ''.
    directory : str, optional
        The working directory. The default is os.getcwd().

    Returns
    -------
    pcoa_md : pandas.DataFrame
        Metadata table with 3 axes of PCoA coordinates added
        (pcoa_1, pcoa_2, pcoa_3).
        Sample IDs are the index, metadata parameters are columns.
    expvals : list of float
        Percent of variance explained by each of the PCoA axes
        (first value is first axis, and so on)
    pcoa_minmax : pandas.DataFrame
        The min and max values of the PCoA coordinates for each of the first
        three PCoA axes. Each row is an axis (first row, pcoa_1).
        This is used to keep the plot axis ranges consistent.

    '''
    # Import PCoA coordinates
    # Open user input file dialog to pick file
    fname, directory, pcoa_coords = fileGet(
        "Select the " + title + " PCoA ordinate file", directory = directory,
        file_type='tsv', header_row=None, skiprows=9, index_col=0)
    # Get file directory info for later saving of plots
    froot, fext = os.path.splitext(fname)
    
    # Get the amount of variance explained
    with open(fname, 'r') as f:
        lines = f.readlines()
        expline = lines[4].strip()
        expvals = expline.split('\t')

    # Add PCoA coordinates to metadata table
    pcoa_coords_3 = pcoa_coords[[1,2,3]]
    pcoa_coords_3.columns=['pcoa_1','pcoa_2','pcoa_3']
    pcoa_md = pd.merge(
        pcoa_coords_3, metadata, left_index = True, right_index = True)
    
    # Determine the axis limits (min/max coordinates)
    pcoa_minmax = pd.DataFrame(
        {'min': pcoa_coords_3.min(), 'max': pcoa_coords_3.max()})
    
    
    # Return compiled metadata file and explained variance
    return pcoa_md, expvals, pcoa_minmax
    
    


def plot2Dpcoa(
        metadataWpcoa, pcoa_expvals, title = '', xlim = None, ylim = None,
        marker_map = False, md_marker = None,
        color_map = False, md_color = None,
        size_map = False, md_size = None):
    '''
    Creates a 2D PCoA plot of data using Matplotlib
    

    Parameters
    ----------
    metadataWpcoa : pandas.dataframe
        Dataframe with sample IDs as the index, metadata parameters as columns,
        with pcoa_1 and pcoa_2 the pcoa coordinates.
    pcoa_expvals : list of float
        List of the fraction of variance explained by each pcoa dimension.
    title : str, optional
        Plot title. Default is none.
    xlim : list of float, optional
        X-axis limits. Default is None, which does not specify a limit.
    ylim : list of float, optional
        Y-axis limits. Default is None, which does not specify a limit.
    marker_map : dict of str, optional
        Dictionary mapping metadata values to marker types. Default is False.
        If a map is provided, md_marker must be defined.
    md_marker : str or None, optional
        Metadata column used to define the markers. Default is None.
    color_map : dict of str, optional
        Dictionary mapping metadata values to colors. Default is False.
        If a map is provided, md_color must be defined.
    md_color : str or None, optional
        Metadata column used to define the colors. Default is None.
    size_map : dict of int, optional
        Dictionary mapping metadata values to marker sizes. Default is False.
        If a map is provided, md_size must be defined.
    md_size : str or None, optional
        Metadata column used to define the marker size. Default is None.

    Returns
    -------
    None

    '''
    
    # Set up the plot
    fig, ax = plt.subplots(figsize = (8,6))

    # Layer on the markers
    if marker_map:
        for mkr in marker_map:
            d_m = metadataWpcoa[metadataWpcoa[md_marker]==mkr]
            
            # Layer on the colors
            if color_map:
                for clr in color_map:
                    d_c = d_m[d_m[md_color]==clr]
                    
                    # Layer on the different marker sizes
                    if size_map:
                        for sz in size_map:
                            d_s = d_c[d_c[md_size]==sz]
                            
                            plt.scatter(
                                d_s['pcoa_1'], d_s['pcoa_2'],
                                c = color_map[clr],
                                marker = marker_map[mkr],
                                s = size_map[sz])
                    else:
                        plt.scatter(
                            d_c['pcoa_1'], d_c['pcoa_2'],
                            c = color_map[clr],
                            marker = marker_map[mkr])
            else:
                plt.scatter(
                    d_m['pcoa_1'], d_m['pcoa_2'],
                    marker = marker_map[mkr])
    else:
        plt.scatter(metadataWpcoa['pcoa_1'], metadataWpcoa['pcoa_2'])

    # Add the plot labels
    plt.title(title)
    plt.xlabel(
        'PCoA-1 (' + f'{float(pcoa_expvals[0])*100:.1f}'
        + '% variance explained)')
    plt.ylabel(
        'PCoA-2 (' + f'{float(pcoa_expvals[1])*100:.1f}'
        + '% variance explained)')
    
    # Fix the axis limits
    if xlim is not None and len(xlim) == 2:
        plt.xlim(xlim.iloc[0], xlim.iloc[1])
    if ylim is not None and len(ylim) == 2:
        plt.ylim(ylim.iloc[0], ylim.iloc[1])

    # Add the legends
    
    # Legend for marker
    if marker_map:
        marker_legend = [
            mlines.Line2D(
                [], [], color = 'black', marker = marker_map[mkr],
                linestyle = 'None', markersize = 8, label = mkr)
            for mkr in marker_map
        ]
        legend_mkr = plt.legend(
            handles = marker_legend, title = md_marker,
            loc='upper left', frameon = False, bbox_to_anchor = (1.01, 1))
        
    # Legend for color
    if color_map:
        color_legend = [
            mpatches.Patch(
                color = color_map[clr], label = clr)
            for clr in color_map
        ]

        legend_clr = plt.legend(
            handles = color_legend, title = md_color,
            loc='upper left', frameon = False, bbox_to_anchor = (1.01, 0.7))
        
    # Legend for size
    if size_map:
        size_legend = [
            plt.scatter([], [], s=size_map[sz], color='k', label=sz)
            for sz in size_map
        ]
        
        legend_size = plt.legend(
            handles = size_legend, title = md_size,
            loc='upper left', frameon = False, bbox_to_anchor = (1.01, 0.3))
    
        # Add the legends to the axes
        ax = plt.gca()
        ax.add_artist(legend_mkr)
        ax.add_artist(legend_clr)
    else:
        ax = plt.gca()
        ax.add_artist(legend_mkr)
        
    
    
    plt.subplots_adjust(right=0.8)

    plt.show()



#%%

####################
# MAIN FUNCTION
####################

# Import metadata
fname, directory, metadata = fileGet(
    "Select the metadata file", file_type='tsv', directory = os.getcwd())

# Import pcoa ordinates and prep metadata tables
pcoa_md_16S, expvals_16S, pcoa_minmax_16S = prepMetadataWpcoa(
    metadata, title = '16S', directory = directory)
pcoa_md_18S, expvals_18S, pcoa_minmax_18S = prepMetadataWpcoa(
    metadata, title = '18S', directory = directory)
pcoa_md_PP, expvals_PP, pcoa_minmax_PP = prepMetadataWpcoa(
    metadata, title = '18S Primary Producers', directory = directory)
dsetmap = {
    '16S' : [pcoa_md_16S, expvals_16S, pcoa_minmax_16S],
    '18S' : [pcoa_md_18S, expvals_18S, pcoa_minmax_18S],
    'Primary Producers' : [pcoa_md_PP, expvals_PP, pcoa_minmax_PP]
    }


####################
# Plot 1:   Samples from each fraction split by
#           month (color) and horizon (symbol)
#           to look for differences by month

# Seperate plots for each gene (16S, 18S)
for d in dsetmap:
    
    # Retrieve data for the gene
    md = dsetmap[d][0]
    expvals = dsetmap[d][1]
    xlim = dsetmap[d][2].loc['pcoa_1',['min','max']]
    ylim = dsetmap[d][2].loc['pcoa_2',['min','max']]
    
    # Only Chukchi Sea samples
    cs = md[md['loc'].isin(['CS', 'JY10', 'JY11'])]
    
    # Plot all fractions together
    plot2Dpcoa(
        cs, expvals, title = d + ' All fractions Weighted Unifrac',
        xlim = xlim, ylim = ylim,
        marker_map = marker_map_horizon, md_marker = 'horizon',
        color_map = color_map_month, md_color = 'month',
        #size_map = size_map_template, md_size = 'template'
        # for cDNA-only plots, comment out the line above
        )
    
    # Seperate plot for each fraction
    for m in materialmap:
        # Extract the samples
        df = cs[cs['material']==m]
    
        plot2Dpcoa(
            df, expvals,
            title = d + ' ' + materialmap[m] + ' Weighted Unifrac',
            xlim = xlim, ylim = ylim,
            marker_map = marker_map_horizon, md_marker = 'horizon',
            color_map = color_map_month, md_color = 'month',
            #size_map = size_map_template, md_size = 'template'
            # for cDNA-only plots, comment out the line above
            )


####################
# Plot 2:   Plot all months together, each fraction seperately with
#           month (color) and horizon (symbol)
#           to look for month differences

# Seperate plots for each gene (16S, 18S)
for d in dsetmap:
    
    md = dsetmap[d][0]
    expvals = dsetmap[d][1]
    xlim = dsetmap[d][2].loc['pcoa_1',['min','max']]
    ylim = dsetmap[d][2].loc['pcoa_2',['min','max']]
    
    # Plot all horizons together
    # Only Chukchi Sea samples
    cs = md[md['loc'].isin(['CS', 'JY10', 'JY11'])]
    plot2Dpcoa(
        cs, expvals, title = d + ' All months Weighted Unifrac',
        xlim = xlim, ylim = ylim,
        marker_map = marker_map_horizon, md_marker = 'horizon',
        color_map = color_map_fraction, md_color = 'material',
        #size_map = size_map_template, md_size = 'template'
        # for cDNA-only plots, comment out the line above
        )
    
    # Seperate plot for each month
    for f in monthmap:
        # Extract the samples
        df = cs[cs['loc_id']==f]
    
        plot2Dpcoa(
            df, expvals,
            title = d + ' ' + monthmap[f] + ' Weighted Unifrac',
            xlim = xlim, ylim = ylim,
            marker_map = marker_map_horizon, md_marker = 'horizon',
            color_map = color_map_fraction, md_color = 'material',
            #size_map = size_map_template, md_size = 'template'
            # for cDNA-only plots, comment out the line above
            )




####################
# EXAMPLE PLOTS

# Matplotlib 2D plot
# fig = plt.figure()
# plt.scatter(pcoa_coords[1], pcoa_coords[2])

# Matplotlib 3D plot
# Matplotlib: https://matplotlib.org/stable/gallery/mplot3d/scatter3d.html#sphx-glr-gallery-mplot3d-scatter3d-py
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(pcoa_coords[1], pcoa_coords[2], pcoa_coords[3])

# Using Plotly
# Plotly: https://plotly.com/python/3d-scatter-plots/
# Problem: symbols in 3D projection are really limited, making Plotly not
#           really more useful than EMPeror.
# import plotly.express as px
# fig = px.scatter_3d(
#    pcoa_md, x='pcoa_1', y='pcoa_2', z='pcoa_3',
#    color = 'month', symbol = 'horizon')
# fig.write_html(froot + '.html')

