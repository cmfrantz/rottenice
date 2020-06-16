# -*- coding: utf-8 -*-
"""
Created on Mon May  4 10:38:36 2020

@author: cariefrantz


Script generates a set of heatmaps from chemical metadata
It was created as part of the Rotten Ice Project
It uses the Master Bio-Chem Sample Sheet Excel file 'Metadata (Raw)' tab as input

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
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from tkinter import filedialog
from tkinter import *

# This tells matplotlib to save text in a format that Adobe Illustrator can read
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


#####
# VARIABLES

excelsheet = 'Metadata (Raw)'   # Tab in the excel sheet where the data is compiled

site = 'CS'
sample_types = ['HT','HM','HB','IT','IM','IB','BT','BM','BB','PW','P1','P2','SW','Blank','Blank2']  # Sample types to plot (and their order)
sample_groupsep = [3,6,9,13]  # Groups the sample types into subgroups on plot

#months = ['M','JN','JY10','JY11']
months = ['M','JN','JY11','JY10']  # Months to plot

plotvars = {
        'physical parameters':  ['temperature', 'salinity_insitu', 'salinity', 'bulk_density'],
        'sediment':             ['SPM', 'SedLoad', 'sterivex_vol_filtered'],
        'chemical parameters':  ['pH', 'nitrogen', 'CN'],
        'carbon fractions':     ['DOC', 'POC', 'pEPS', 'SPM'],
        'bacterial counts':     ['bact_cell_ct', 'CTC', 'bact_active_cell_ct'],
        'algae counts':         ['diatom_ct', 'phyto_ct_select', 'phyto_ct_all', 'phyto_ct_other'],
        'photosynthesis':       ['Chl', 'Phaeo', 'FoFa', 'PAM']
        }

maxplotcols = 4
scalewidth = 3.5      # Factor scales the plot height based on number of plot columns
scaleheight = 7     # Factor scales the plot width based on number of plot rows
clrmap = 'viridis'    # Colormap to use, see https://matplotlib.org/tutorials/colors/colormaps.html

filename = 'RottenIce_ChemHeatmaps_All'

####
# FUNCTIONS

# Import chemistry file
def getFile():
    while True:
        try:
            root = Tk()                     # Opens user input window
            root.filename = filedialog.askopenfilename(initialdir = "/", title = 'Select file', filetypes = [('Excel', '*.xls, *.xlsx')])     # Ask user for file
            filename = root.filename        # Retrieves the selected filename
            print (filename)                # Prints the filename
            root.destroy()                  # Closes the user input window
            print('Loading file...')    
            data = pd.read_excel(filename, sheet_name = excelsheet, header = 2, index_col = 0)        # Loads in data from filename
            break
        except ValueError:
            print("Something was wrong with the file " + filename + ".  Try again...")

    return data, filename


# Build the measurement array
def buildArray(data, var):
    valueMatrix = pd.DataFrame(columns = months, index = sample_types)
    if var in data.columns:
        # Loop through months
        for month in months:
            # Generate list of sample names
            samples = [genSampleName(month, site, x) for x in sample_types]
            valueMatrix[month] = data.loc[samples][var].values
                        
        # replace any invalid values
        valueMatrix = valueMatrix.replace('na', np.nan)
    else: print(var + ' is not a valid variable')
    
    return valueMatrix
                

# Generates the samename
def genSampleName(month, site, sample):
    if 'Blank' in sample:
        if 'JY' in month:   name = 'JY-' + sample
        else:               name = month + '-' + sample
    else:
        if 'JY' in month:  
            if sample in ['BT', 'BM', 'BB']:    sample = 'B'
            elif sample in ['P1', 'P2']:        sample = 'Drain'
            name = month + '-' + sample
        else:               name = month + '-' + site + '-' + sample
    return name
            


# Generates the seperation lines between fraction types                 
def getLineCoord(y):
    xvals = [-0.5, len(months)-0.5]
    yvals = [y-0.5, y-0.5]    
    return xvals, yvals



# Build the plots
def buildSubplot(ax, var, matrix_val, matrix_std):
    ax.imshow(matrix_val, aspect='auto', cmap=cmap)
      
    # Loop over data dimensions and create text annotations.
    for i in range(len(months)):
        for j in range(len(sample_types)):
            val = matrix_val.iloc[j][i]
            std = matrix_std.iloc[j][i]
            if      pd.isnull(val):    valtext = 'nan'
            elif    pd.isnull(std):    valtext = format(val,'g')
            else:   valtext = '{: 0.2g}\n\u00B1{:0.2g}'.format(val,std)
            ax.text(i, j, valtext,
                    ha="center", va="center", 
                    color="lightgray", fontsize=8)
    
    # Add lines to seperate sample groups
    for line in sample_groupsep:
        ax.plot(getLineCoord(line)[0], getLineCoord(line)[1], color = 'k', linewidth = 1)
        
    ax.set_title(var)

    
####
# MAIN FUNCTION
if __name__ == '__main__':
    
    # Get data from the user
    [data, infilename] = getFile()
    
        
    # Plot the data
    # Set up the subplots
    cmap = plt.get_cmap(clrmap)
    fig, axs = plt.subplots(len(plotvars), maxplotcols, 
                            sharey=True, sharex=True, 
                            figsize = [scalewidth*len(months), scaleheight*len(plotvars)])

    # Generate each subplot
    for i in np.arange(len(plotvars)):
        varset = list(plotvars)[i]
        for j in np.arange(len(plotvars[varset])):
            var = list(plotvars[varset])[j]
            # Build the matrices
            matrix_val = buildArray(data, var)
            matrix_std = buildArray(data, var+'_std')
            # Build subplot
            buildSubplot(axs[i][j], var, matrix_val, matrix_std)
    
    # Set up ticks...
    axs[0,0].set_xlim(-0.5,len(months)-0.5)
    axs[0,0].set_xticks(np.arange(len(months)))
    axs[0,0].set_yticks(np.arange(len(sample_types)))
    # ... and label them with the respective list entries
    axs[0,0].set_xticklabels(months)
    axs[0,0].set_yticklabels(sample_types)
    # Rotate the tick labels and set their alignment.
    plt.setp(axs[0,0].get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
    
    
    fig.tight_layout()
    plt.show()

    # Save the image
    fig.savefig(os.getcwd() + '\\' + filename + ".svg", transparent=True)
    fig.savefig(os.getcwd() + '\\' + filename + ".pdf", transparent=True)
        