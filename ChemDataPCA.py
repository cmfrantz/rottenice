#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""Rotten Ice Project: ChemDataPCA

Created on Wed Jun 26 15:50:09 2019
@author: cariefrantz
@project: RottenIce

CREATE 2D PCA PLOTS FROM IMPORTED METADATA
This script creates 2D PCA plots from imported metadata

This script was created as part of the Rotten Ice Project

Arguments:  None

Requirements:      
    Metadata table (csv)
        where rows = samples, columns = metadata characteristics
        header row and index column are specified in the variables below

Example in command line:
    python ChemDataPCA.py

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


Copyright (C) 2019  Carie M. Frantz

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

from matplotlib import pyplot as plt
import matplotlib
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import RottenIceModules


####################
# VARIABLES
####################

# List of parameters that can be used for PCA
featlist = ["temperature", "salinity", "bulk_density",
            "DOC", "POC", "pEPS", "nitrogen", "CN",
            "Chl", "Phaeo", "FoFa", "PAM",
            "bact_cell_ct", "CTC", "bact_active_cell_ct",
            "diatom_ct", "phyto_ct_all", "phyto_ct_other"]

# List of months to include
months = ['M','JN','JY11']
# List of fractions to include
fractions = ['HT', 'HM', 'HB']

#%%

####################
# FUNCTIONS
####################

def prepData(metadata):
    '''Prepare data from user-selected features'''
    # Ask user for features to plot
    features = input('Enter the list of features to include in the PCA. '
                     + 'Seperate items with commas as in the example below: '
                     + '\n > temperature, salinity, DOC, nitrogen \n'
                     + 'Options include: ' + ', '.join(featlist)
                     + '\n > ')
    features = features.split(', ')    
    return metadata[features], features


def calcPCA(data, features):
    '''Perform principal components analysis'''
    # Normalize data for PCA (PCA is affected by scale)
    vals = StandardScaler().fit_transform(data[features])
    
    # Calculate 2D PCA
    pca = PCA(n_components = 2)     # 2-component PCA
    principalComponents = pca.fit_transform(vals)
    explained_var = pca.explained_variance_ratio_
    data.loc[:,'PC1'] = principalComponents[:,0]
    data.loc[:,'PC2'] = principalComponents[:,1]
        
    return data, explained_var
        

def plotResults(data, features, datafmt, explained_var, dirPath):
    # Create the plot & label axes
    fig = plt.figure(figsize = (6,4))
    ax = fig.add_subplot(1,1,1)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
    
    # Add title
    ax.set_title('2 component PCA:\n' + ', '.join(features))
    
    # Add axis titles
    ax.set_xlabel('PC1 (' + str(round(explained_var[0]*100, 1)) + '%)',
                  fontsize = 15)
    ax.set_ylabel('PC2 (' + str(round(explained_var[1]*100, 1)) + '%)',
                  fontsize = 15)
    
    # Plot each point  
    for sample in data.index:
        ax.plot(data.loc[sample, 'PC1'],                 # x value
                data.loc[sample, 'PC2'],                 # y value
                linestyle = 'None',
                **datafmt[sample])                       # marker size
    
    # Add legend
    ax.legend(data.index, loc = 'upper left', bbox_to_anchor = (1, 1))
    
    # Show and save figure
    plt.show()
    fig.savefig(dirPath + "\\metadata_PCA_" + '-'.join(features) + ".pdf", transparent=True) # Save figure

#%%
####################
# MAIN FUNCTION
####################
if __name__ == '__main__':
    # Generate list of samples and their plot formatting info
    samples, markers = RottenIceModules.genSamples_w_Markers(months, fractions)
    datafmt = dict(zip(samples, markers))
    
    # Retrieve the metadata file
    filename, dirpath, metadata = RottenIceModules.fileGet(
        'Select metadata file', tabletype = 'metadata')
    metadata = metadata.replace('na', np.nan).loc[samples]
    
    # Loop for as many plots as user wants to build
    while True:
        # Prep the data for PCA based on user input for which variables to use
        data, features = prepData(metadata)
        # Perform PCA
        data, explained_var = calcPCA(data.copy(), features)
        # Plot results
        plotResults(data, features, datafmt, explained_var, dirpath)
        # Ask user if they want to build another plot
        if input('Build another plot? Y/N  > ') == 'N':
            break
