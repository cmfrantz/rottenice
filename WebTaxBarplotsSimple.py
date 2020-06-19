#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""Rotten Ice Project: WebTaxBarplotsSimple

Created on Mon May 11 17:49:35 2020
@author: cariefrantz
@project: RottenIce

BUILDS INTERACTIVE (BOKEH) TAXONOMIC BARPLOTS: SIMPLE VERSION
This script builds an HTML webpage with interactive taxonomic bar plots.
It builds side-by-side plots displaying absolute and relative abundances.
Each taxonomic level is a tab in the webpage.

This version of the WebTaxBarplots script was used for visualizing
phytoplankton taxonomy data from microscopy counts.
By modifying variables at the top of the code, it can easily be adapted for
any OTU table.

There is another version (WebTaxBarplotsStacked) that was used for sequencing
data that produces stacked plots in order to compare DNA vs. cDNA results.

This script was created as part of the Rotten Ice Project


Arguments:  None

Requirements:      
    OTU table (csv) for algae taxonomy calls
        where rows = samples, columns = taxononomic call
        header row and index column are specified in the variables below

Example in command line:
    python WebTaxBarplotsSimple.py

Dependencies Install:
    sudo apt-get install python3-pip python3-dev
    pip install tkinter
    pip install progress
    pip install pandas
    pip install numpy
    pip install matplotlib
    pip install math

You will also need to have RottenIceModules.py downloaded to the same \
    directory as this script as it contains modules that this script calls.

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
from tkinter import *
from tkinter import filedialog
from progress.bar import Bar

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from math import pi

from bokeh.io import output_file, show, export_svgs
from bokeh.layouts import row
from bokeh.plotting import figure
from bokeh.models import Legend, LegendItem, Panel, Tabs

import RottenIceModules


####################
# VARIABLES
####################

# Input file variables
sample_row = 0      # Row in ASV table files containing unique sample names
OTU_col = 0         # Column in metadata file containing unique ASV names

# Set the max taxonomic depth to plot (e.g., deepest meaningful level)
max_level = 14

# Colors
cmap = 'viridis'    # Change this to change the default colormap
cmaptext = ''''Choose a colormap, e.g.,
viridis, plasma, prism, nipy_spectral, jet, gist_earth, etc.
See https://matplotlib.org/tutorials/colors/colormaps.html for a full list. 
> '''

# Plot formatting
pltht = 380         # height of the actual plot
pltw = 500          # width of the actual plot
minpadding = 5      # padding between plots and borders
pltbottompad = 100  # padding for x axis label
maxkeywidth = 925   # maximum width needed for the key
keyvalht = 24       # height taken up by each key entry

# Plot properties based on plot type
displayprops = {
    'L' : {
        'title'     : 'Absolute',
        'ylabel'    : 'Phytoplankton count (cells/mL)',
        'legend'    : False,
        'lpadding'  : minpadding,
        'rpadding'  : minpadding
        },
    'R' : {
        'title'     : 'Relative',
        'ylabel'    : 'Relative abundance',
        'legend'    : True,
        'lpadding'  : minpadding * 4,
        'rpadding'  : maxkeywidth
        }
    }


####################
# FUNCTIONS
####################



# Generate the data to use in the plots
def genPltdata(data, level, samples):
    '''Processes raw data to produce two datasets:
        - raw counts (absolute data)
        - relative abundances (relative data)'''
      
    # Reset taxonomic ID based on the selected level
    cols = levelCols(level)
    #taxlist = data['taxonomy'].values
    for val in list(data.index):
        #taxlist[i] = '>'.join(list(data.iloc[i][cols]))
        data.at[val, 'taxonomy'] = '>'.join(list(data.loc[val][cols]))
    
    # Create normalized data
    normdata = data.copy()
    for sample in samples:
        normdata[sample] = normdata[sample].values/np.sum(normdata[sample])
    
    # Save data dicts
    pltdata = genPltDict(data, samples)
    normpltdata = genPltDict(normdata, samples)
        
    return pltdata, normpltdata


# Create dict for plot
def genPltDict(data, samples):
    '''Builds dictionaries to hold the data and settings for each plot'''
    pltdata = {'samples' : samples}
    
    taxlist = genTaxlist(data['taxonomy'].values)
    for value in taxlist:
        rows = data.loc[data['taxonomy']==value]
        if rows.shape[0] == 1:
            pltdata[value] = list(rows[samples].to_numpy()[0])
        else:
            mat = rows[samples].to_numpy()
            pltdata[value] = list(np.sum(mat, axis=0))
     
    return pltdata


def buildPlot(pltdata, samples, level, colors, pltkey, legend=False):
    '''Builds the plots'''
    # Determine plot height
    minht = max(pltht + minpadding + pltbottompad, keyvalht * (len(pltdata)-1))
    
    # Build figure
    p = figure(x_range = samples, 
               plot_height = minht, 
               min_border_top = minpadding, 
               min_border_bottom = minht - pltht - minpadding,
               plot_width = (pltw + displayprops[pltkey]['lpadding']
                             + displayprops[pltkey]['rpadding']),
               min_border_left = displayprops[pltkey]['lpadding'],
               min_border_right = displayprops[pltkey]['rpadding'],
               title = ('Taxonomy level ' + str(level) + ' '
                        + displayprops[pltkey]['title']), 
               x_axis_label = 'Sample',
               y_axis_label = displayprops[pltkey]['ylabel'],
               toolbar_location = 'above',
               tools='hover,box_zoom,wheel_zoom,pan,reset,save', 
               tooltips='$name: @$name')
    
    # Add the bars corresponding to each taxon
    taxlist = list(pltdata)[1:]
    vbars = p.vbar_stack(taxlist, x='samples', width = 0.9, color = colors, 
                         source=pltdata)
    
    # Format plot
    p.y_range.start = 0
    p.x_range.range_padding = 0.1
    p.xgrid.grid_line_color = None
    p.xaxis.major_label_orientation = pi/4
    p.axis.minor_tick_line_color = None
    
    # Add legend
    if displayprops[pltkey]['legend'] == True:
        legend = genLegend(taxlist, vbars)
        p.add_layout(legend, 'right')
        
    # Save figure
    p.output_backend = 'svg'
    export_svgs(p, filename = (filename[0:-4] + '_barplot_'
                               + displayprops[pltkey]['title'] + '_L'
                               + str(level) + '.svg'))
    return p


def genLegend(taxlist, vbars):
    '''Draws a legend'''
    legenditems = []
    # Create a legend item for each taxon
    for i in np.arange(len(taxlist)):
        legenditems.append(LegendItem(label = taxlist[i],
                                      renderers = [vbars[i]]))
    # Build the legend from the legend items    
    legend = Legend(items = legenditems, location = 'top_right', 
                    glyph_height = 15, glyph_width = 15,
                    label_text_font_size = '9px',)
    return legend
#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    
    # Get and format file
    filename, filepath, data = RottenIceModules.fileGet(
        'Select OTU table file')
    data, samples = RottenIceModules.formatOTUtableData(
        data, max_level = max_level)
          
    # Get the colors (disable to use default colors)
    # cmap = input(cmaptext)
    
    # Build tabs of plots
    tabs = []
    for level in np.arange(1,max_level+1):
        print('Generating Level ' + str(level) + ' plot...')
        [pltdata, normpltdata] = genPltdata(data, level, samples)
        
        # Plot the data
        output_file(filename[0:-4] + '_barplots_L' + str(level) + '_' + cmap
                    +'.html', title = filename[0:-4] + ' barplots')
        colors = genColormap(len(pltdata)-1, cmap)
        lplt = buildPlot(pltdata, samples, level, colors, 'L')
        rplt = buildPlot(normpltdata, samples, level, colors, 'R')
        grid = row(lplt, rplt)
        
        # Enable this for troubleshooting to display individual plots,
        # but tabs won't stitch together into HTML file if it is enabled.
        # show(grid)
        
        # Create the tab
        panel = Panel(child = row(lplt, rplt), title = 'Level ' + str(level))
        tabs.append(panel)
    
    # Save tabs and HTML file
    print('Saving HTML file...')
    tabs = Tabs(tabs=tabs)
    show(tabs)