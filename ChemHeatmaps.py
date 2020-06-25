#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""Rotten Ice Project: ChemHeatmaps

Created on Mon May  4 10:38:36 2020
@author: cariefrantz
@project: RottenIce

GENERATE HEATMAPS FROM CHEMICAL METADATA
This script generates a set of heatmaps for chemical metadata

This script was created as part of the Rotten Ice Project

Arguments:  None

Requirements:      
    Metadata table (csv)
        where rows = samples, columns = metadata characteristics
        header row and index column are specified in the variables below

Example in command line:
    python ChemHeatmaps.py


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
import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import RottenIceModules
import RottenIceVars

####################
# VARIABLES
####################

# Data to include
site = 'CS'
sample_types = ['HT','HM','HB',     # Sample types to plot (and their order)
                'IT','IM','IB',
                'BT','BM','BB',
                'PW','P1','P2','SW',
                'Blank','Blank2']
sample_groupsep = [3,6,9,13]        # Groups the sample types into subgroups

months = ['M','JN','JY11','JY10']   # Months to plot
                                    # Note: the order of the JY samples
                                    # is intentional. JY11 "clean" floe was
                                    # most similar to the ice character from
                                    # prior months so is presented first for
                                    # ease of comparison
                                    
plotvars = {                        # Variables to plot.
                                    # Each set plots in its own row.
        'physical parameters':  ['temperature', 'salinity_insitu', 'salinity',
                                 'bulk_density'],
        'sediment':             ['SPM', 'SedLoad', 'sterivex_vol_filtered'],
        'chemical parameters':  ['pH', 'nitrogen', 'CN'],
        'carbon fractions':     ['DOC', 'POC', 'pEPS', 'SPM'],
        'bacterial counts':     ['bact_cell_ct', 'CTC', 'bact_active_cell_ct'],
        'algae counts':         ['diatom_ct', 'phyto_ct_select',
                                 'phyto_ct_all', 'phyto_ct_other'],
        'photosynthesis':       ['Chl', 'Phaeo', 'FoFa', 'PAM']
        }

# Plot formatting variables
maxplotcols = 4     # Maximum number of columns to show on the  plot
scalewidth = 3.5    # Factor scales plot height based on number of plot cols
scaleheight = 7     # Factor scales plot width based on number of plot rows
cmap = RottenIceVars.cmap # Colormap to use; see:
                    # https://matplotlib.org/tutorials/colors/colormaps.html

# Output file info
file_info = RottenIceVars.file_sets['chem_heatmaps']
subtitle_text = ('Created from compiled project metadata using the script '
                 + '<a href="https://github.com/cmfrantz/rottenice">'
                 + 'ChemHeatmaps.py</a>. Analysis done by C. Frantz, '
                 + 'June 2020.')


####################
# FUNCTIONS
####################

def readMetadata():
    '''This function reads the sample metadata in from csv file'''
    filename, directory, metadata = RottenIceModules.fileGet(
        'Select metadata file', tabletype = 'metadata')
    # replace invalid values
    metadata = metadata.replace('na', np.nan)
    return metadata, filename, directory


def buildArray(data, var):
    '''Generates the array of measurements selected for the plot'''
    # Create a blank matrix and fill it with data for the selected variable
    valueMatrix = pd.DataFrame(columns = months, index = sample_types)
    if var in data.columns:
        # Loop through months
        for month in months:
            # Generate list of sample names
            samples = [RottenIceModules.genSampleName(month, fraction) 
                       for fraction in sample_types]
            # Grab and paste in the data for that month
            vals = [float(val) for val in data.loc[samples,var].values]
            valueMatrix.loc[:,month] = vals
    else: print(var + ' is not a valid variable')
    return valueMatrix
            
               
def getLineCoord(y):
    '''Generates the seperation lines between fraction types'''
    xvals = [-0.5, len(months)-0.5]
    yvals = [y-0.5, y-0.5]    
    return xvals, yvals


def buildSubplot(ax, var, matrix_val, matrix_std, cmap):
    '''Constructs each variable's heatmap and adds it to a subplot space'''
    # Plot heatmap
    ax.imshow(matrix_val, aspect='auto', cmap=cmap)
      
    # Annotate heatmap by looping through each dimension of the matrix plotted
    # Loop through columns (months)
    for i in range(len(months)):
        # Loop through rows (sample fractions)
        for j in range(len(sample_types)):
            # Retrieve values
            val = matrix_val.iloc[j][i]
            std = matrix_std.iloc[j][i]
            # Add the annotation text
            if      pd.isnull(val):    valtext = 'nan'
            elif    pd.isnull(std):    valtext = format(val,'g')
            else:   valtext = '{: 0.2g}\n\u00B1{:0.2g}'.format(val,std)
            ax.text(i, j, valtext,
                    ha="center", va="center", 
                    color="lightgray", fontsize=8)
    
    # Add lines to seperate sample groups
    for line in sample_groupsep:
        ax.plot(getLineCoord(line)[0], getLineCoord(line)[1],
                color = 'k', linewidth = 1)
    
    # Add subplot title
    ax.set_title(var)

 #%% 
####
# MAIN FUNCTION
if __name__ == '__main__':
    
    # Get data from the user
    data, infilename, directory = readMetadata()
        
    # Plot the data
    # Set up the subplots
    print('Setting up plot field...')
    cmap = plt.get_cmap(cmap)     # Load colormap
    fig, axs = plt.subplots(len(plotvars), maxplotcols, 
                            sharey=True, sharex=True, 
                            figsize = [scalewidth*len(months),
                                       scaleheight*len(plotvars)])

    # Generate each subplot
    # Loop through each variable set (plot rows)
    for i, varset in enumerate(plotvars):
        # Loop through each variable (plot columns)
        for j, var in enumerate(plotvars[varset]):
            print('Plotting ' + var)
            # Build the matrices of averages & standard deviations
            matrix_val = buildArray(data, var)
            matrix_std = buildArray(data, var + '_std')
            # Build subplot
            buildSubplot(axs[i][j], var, matrix_val, matrix_std, cmap)
    
    # Set up x and y ticks for each month and sample fraction ...
    print('Saving image...')
    axs[0,0].set_xlim(-0.5,len(months)-0.5)
    axs[0,0].set_xticks(np.arange(len(months)))
    axs[0,0].set_yticks(np.arange(len(sample_types)))
    # ... and label them with the respective month and sample fraction
    axs[0,0].set_xticklabels(months)
    axs[0,0].set_yticklabels(sample_types)
    # Rotate the tick labels and set their alignment
    plt.setp(axs[0,0].get_xticklabels(), ha="right",
             rotation=90, rotation_mode="anchor")
    
    # Display the image
    fig.tight_layout()
    plt.show()

    # Save the image
    filename = file_info['pfx']
    fig.savefig(directory + '\\' + filename + ".svg", transparent=True)
    fig.savefig(directory + '\\' + filename + ".pdf", transparent=True)
    
    # Generate HTML
    RottenIceModules.genHTMLfile(
        directory + '\\' + file_info['land_page'],
        file_info['title'], subtitle_text,
        filename + '.svg',
        alt_text = 'Heatmap table of metadata')
        