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
    pip install os
    pip install tkinter
    pip install pandas
    pip install numpy
    pip install matplotlib
    pip install sklearn


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
import os
from tkinter import filedialog
from tkinter import *
import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import matplotlib
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


####################
# VARIABLES
####################

# List of parameters that could be used for PCA
featlist = ["temperature", "salinity", "bulk_density",
            "DOC", "POC", "pEPS", "nitrogen", "CN",
            "Chl", "Phaeo", "FoFa", "PAM",
            "bact_cell_ct", "CTC", "bact_active_cell_ct",
            "diatom_ct", "phyto_ct_all", "phyto_ct_other"]

# Input file variables
metadata_row = 2                    # Row containing unique metadata names
sample_col = 0                      # Column containing unique sample names

# Plot formatting variables
site = 'CS' # Restrict analysis to CS samples
# Assign color by sample month
# (this also serves as the list of months to plot)
months = {
    'M'     : 'blue',
    'JN'    : 'green',
    'JY11'  : 'red'
    }
# Assign marker by sample horizon
# (this also serves as the list of horizons to plot)
horizons = {
    'HT'    : '^',
    'HM'    : 'o',
    'HB'    : 'v'
    }


####################
# FUNCTIONS
####################
def genSampleName(month, horizon):
    '''Generates the sample name from month and horizon'''
    if 'JY' in month:
        sample = month + '-' + horizon
    else:
        sample = month + '-' + site + '-' + horizon
    return sample


def readMetadata():
    '''This function reads the sample metadata in from csv file'''
    # User input dialog to locate the metadata file
    root = Tk()                     # Opens user input window
    root.filename = filedialog.askopenfilename(     # Ask user for file
                            initialdir = "/",
                            title = 'Select metadata file',
                            filetypes = [('CSV', '*.csv')])
    filename = root.filename        # Retrieves the selected filename
    root.destroy()                  # Closes the user input window
    dirPath = os.path.dirname(filename)     # Directory
    print('Loading ' + filename + '...')
    
    # Read metadata tab
    metadata = pd.read_csv(filename, header = metadata_row,
                           index_col = sample_col)
    
    # replace invalid values
    metadata = metadata.replace('na', np.nan)

    return metadata, filename, dirPath


def prepData(metadata):
    '''Prepare data from user-selected features'''
    # Ask user for features to plot
    features = input('''
*** Enter the list of features to include in the PCA ***
Seperate items with commas as in the example below:
> temperature, salinity, DOC, nitrogen
                     
Options include: ''' + ', '.join(featlist) + '''
> ''')

    features = features.split(', ')
    metadata = metadata[features]
    
    # Generate list of samples and their plot formatting info
    samples = []
    datafmt = pd.DataFrame(index = metadata.index, columns=['color', 'marker'])
    for month in months:
        for horizon in horizons:
            sample = genSampleName(month, horizon)
            datafmt.loc[sample,'color'] = months[month]
            datafmt.loc[sample,'marker'] = horizons[horizon]
            samples.append(sample)
    data = metadata.merge(datafmt, left_index = True, right_index = True)
    data = data.loc[samples]
   
    return data, features


def calcPCA(data, features):
    '''Perform principal components analysis'''
    # Normalize data for PCA (PCA is affected by scale)
    vals = StandardScaler().fit_transform(data[features])
    
    # Calculate 2D PCA
    pca = PCA(n_components = 2)     # 2-component PCA
    principalComponents = pca.fit_transform(vals)
    explained_var = pca.explained_variance_ratio_
    data['PC1'] = principalComponents[:,0]
    data['PC2'] = principalComponents[:,1]
        
    return data, explained_var
        

def plotResults(data, features, explained_var, dirPath):
    # Create the plot & label axes
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    
    # Add title
    ax.set_title('2 component PCA: ' + ', '.join(features))
    
    # Add axis titles
    ax.set_xlabel('PC1 (' + str(round(explained_var[0]*100, 1)) + '%)',
                  fontsize = 15)
    ax.set_ylabel('PC2 (' + str(round(explained_var[1]*100, 1)) + '%)',
                  fontsize = 15)
    
    # Plot each point  
    for sample in data.index:
        ax.scatter(data.loc[sample, 'PC1'],                 # x value
                   data.loc[sample, 'PC2'],                 # y value
                   c = data.loc[sample, 'color'],           # marker color
                   marker = data.loc[sample, 'marker'],     # marker shape
                   s=50)                                    # marker size
    
    # Add legend
    ax.legend(data.index)
    
    # Show and save figure
    plt.show()
    fig.savefig(dirPath + "\\metadata_PCA_" + '-'.join(features) + ".pdf", transparent=True) # Save figure

#%%
####################
# MAIN FUNCTION
####################
if __name__ == '__main__':
    # Retrieve the metadata file
    metadata, filename, dirPath = readMetadata()
    
    # Loop for as many plots as user wants to build
    while True:
        # Prep the data for PCA based on user input for which variables to use
        data, features = prepData(metadata)
        # Perform PCA
        data, explained_var = calcPCA(data, features)
        # Plot results
        plotResults(data, features, explained_var, dirPath)
        # Ask user if they want to build another plot
        if input('Build another plot? Y/N  > ') == 'N':
            break
    
    