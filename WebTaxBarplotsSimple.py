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

This script was created as part of the Rotten Ice Project.


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
import numpy as np
from bokeh.io import output_file, show
from bokeh.layouts import row

import RottenIceModules
import RottenIceVars


####################
# VARIABLES
####################

# Set the max taxonomic depth to plot (e.g., deepest meaningful level)
max_level = 14

# Samples to plot, and how they should be sorted in the plot
samples = ['M-CS-HT', 'M-CS-HM', 'M-CS-HB', 'M-CS-SW',
           'JN-CS-HT', 'JN-CS-HM', 'JN-CS-HB', 'JN-CS-SW',
           'JY10-HT', 'JY10-HM', 'JY10-HB', 'JY10-Drain', 'JY10-PW',
           'JY11-HT', 'JY11-HM', 'JY11-HB',
           'JY11-IT', 'JY11-IM', 'JY11-IB', 'JY11-Drain', 'JY11-SW']

# Plot / HTML file title info
title = 'Phytoplankton taxonomy'
subtitle_text = ('Rotten ice project phytoplankton taxonomic ID barplots. '
                 + 'Taxonomic IDs by Sylvie Lessard (2016). '
                 + 'Taxonomy breakdown done by B. Tattersall (2017) using the '
                 + '<a href = "https://www.algaebase.org/AlgaeBase">'
                 + 'AlgaeBase database</a>. '
                 + 'Data processed by C. Frantz, May 2020 from cleaned tables '
                 + '(plausible swapped sample corrected) using the script '
                 + '<a href="https://github.com/cmfrantz/rottenice">'
                 + 'WebTaxBarplotsSimple.py</a>.')
    
filename_prefix = RottenIceVars.file_sets['algae_barplots']['pfx']

# Plot formatting
cmap = RottenIceVars.cmap # colormap to use; this uses the project default
plot_width = 500        # width of the actual plot
plot_height = 380       # height of the actual plot
min_padding = 5         # padding between plots and borders
plot_bottom_pad = 100   # padding for x axis label
max_key_width = 925     # maximum width needed for the key
key_val_height = 24     # height taken up by each key entry

# Bokeh plot formatters
# Values to the right of the ':' shown without quotes are based on
# plot formatting variables defined above. Those in quotes can be modified.
shared_kwargs = {
    'plot_height'       : plot_height,
    'min_border_top'    : min_padding,
    'x_axis_label'      : 'Sample',
    'toolbar_location'  : 'above',
    'tools'             : 'hover,box_zoom,wheel_zoom,pan,reset,save',
    'tooltips'          : '@samples $name: @$name'
    }

# Plot properties based on plot type
# Values to the right of the ':' shown without quotes are based on
# plot formatting variables defined above. Those in quotes can be modified.
plot_properties = {
    'L' : {
        'title'     : 'Absolute',
        'dataset'   : 'abs',
        'legend'    : False,
        'y_axis_label' : 'Phytoplankton counts (cells/ml)',
        'lpadding'  : min_padding,
        'rpadding'  : min_padding
        },
    'R' : {
        'title'     : 'Relative',
        'dataset'   : 'rel',
        'legend'    : True,
        'y_axis_label' : 'Relative abundance',
        'lpadding'  : min_padding * 4,
        'rpadding'  : max_key_width
        }
    }


#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    
    # Get file
    filename, filepath, data = RottenIceModules.fileGet(
        'Select OTU table file', tabletype = 'OTU_Table')
    # Format file
    data, sample_list = RottenIceModules.formatOTUtableData(
        data, max_level = max_level)
    
    # Find max absolute value for scaling axis
    max_val = np.max(data[samples].sum(axis = 0))
    max_ax_val = RottenIceModules.findMaxAxVal(max_val, 5)
          
    # Get the formatting info
    # Colors
    # cmap = input(cmaptext)    # Disable to use default colors
    tax_columns = RottenIceModules.levelCols(max_level) + ['taxonomy']
    colormap = RottenIceModules.colorMap2Tax(
        data[tax_columns].copy(), max_level, cmap)
    # Kwargs
    for plot in plot_properties:
        props = plot_properties[plot]
        kwargs = {
            'plot_width' : plot_width + props['lpadding'] + props['rpadding'],
            'min_border_left' : props['lpadding'],
            'min_border_right' : props['rpadding'],
            'y_axis_label' : props['y_axis_label'],
            'title'     : props['title']
            }
        kwargs.update(shared_kwargs)
        plot_properties[plot]['kwargs'] = kwargs
    
    # Create shared HTML navigation header
    filenames = RottenIceModules.genFilenamesByLevel(
        filename_prefix, max_level, cmap = cmap)
    nav_html = RottenIceModules.genPlotNavHTML(filenames)
    
    # Build each level's plot
    for level in np.arange(1,max_level+1):
        print('Generating Level ' + str(level) + ' plot...')
        # Prepare the data
        pltdata, normpltdata, colors = RottenIceModules.buildPlotDicts(
            data, level, samples, colormap)
        
        # Generate filename
        output_file(filenames[level-1] +'.html',
                    title = title + ' L' + str(level))
        
        # Create plots
        plots = []
        # Determine plot height
        min_height = max(plot_height + min_padding + plot_bottom_pad,
                    key_val_height * (len(pltdata)-1))
        # Build each plot
        for plot in plot_properties:
            # Gather plot formatting info
            props = plot_properties[plot]
            props['kwargs'].update({
                'plot_height'   : min_height,
                'min_border_bottom' : min_height - plot_height - min_padding
                })
            if props['dataset'] == 'abs':
                dset = pltdata
                ymax = max_ax_val
            elif props['dataset'] == 'rel':
                dset = normpltdata
                ymax = 1
            # Build plot
            plt = RottenIceModules.buildBarPlot(
                dset, colors, ymax, legend = props['legend'],
                **props['kwargs'])
            plots.append(plt)
                  
        # Save page
        grid = row(plots)
        print('Saving HTML file...')
        show(grid)
        
        # Add header to page
        RottenIceModules.addHTMLhead(
            filenames[level-1]+'.html',
            page_title = title,
            page_nav_html = nav_html,
            subtitle_text = ('<b>Taxonomic level ' + str(level) + '</b><br />'
                             + subtitle_text))
