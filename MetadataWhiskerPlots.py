#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'
"""Rotten Ice Project: MetadataWhiskerPlots

Created on Mon Jun 22 09:36:30 2020
@author: cariefrantz
@project: RottenIce

CREATES WHISKER BOX PLOTS COMPARING METADATA VALUES
This script creates whisker box plots from specified metadata values and
    overlays all datapoints in order to enable metadata comparisons.
This script was created as part of the Rotten Ice Project.

Arguments:    None

Requirements:     
    Sample Metadata (measurements) table (csv)
        Where rows = samples, columns = metadata characteristics
        Header row and index column are specified in the variables below
        This is NOT the same metadata table used in QIIME2, which provides
        metadata for each sample-gene-template sequenced. It is the master
        metadata table of environmental measurements for each sample.

Example in command line:
    python MetadataWhiskerPlots.py

Dependencies Install:
    sudo apt-get install python3-pip python3-dev
    pip install numpy
    pip install matplotlib

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

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import RottenIceModules
import RottenIceVars


####################
# VARIABLES
####################

monthcode = {
    'M': 'May',
    'JN': 'June',
    'JY10': 'July 10\nDirty Floe',
    'JY11': 'July 11\nClean Floe'
}

# Variables defining what samples and metadata parameters to plot
# sets = [['M'],['JN'],['JY10'],['JY11'],['M','JN'],['JY10','JY11']]
sets = [['M'],['JN'],['JY10'],['JY11']]
fractions = ['HT','HM','HB']
varlist = [['temperature','salinity_comb','bulk_density'],
           ['SPM', 'nitrogen','CN'],
           ['DOC','POC','pEPS'],
           ['Chl','Phaeo'],
           ['FoFa','PAM'],
           ['bact_cell_ct', 'CTC', 'bact_active_cell_ct'],
           ['phyto_ct_all', 'diatom_ct', 'phyto_ct_other']]
'''
varlist = [['temperature','salinity_comb','bulk_density'],
           ['SPM', 'nitrogen','CN'],
           ['DOC','POC','pEPS'],
           ['Chl','Phaeo'],
           ['FoFa','PAM'],
           ['bact_cell_ct', 'CTC', 'bact_active_cell_ct'],
           ['phyto_ct_all', 'diatom_ct', 'phyto_ct_other'],
           ['gels_total','gels_sm','gels_md'],
           ['gels_lg','gels_xl']]
'''

# Variables defining plot parameters
# Plot size
pltht = 3
pltw = 5
# Variables with numbers that require scientific notation
vars2format = ['bact_cell_ct', 'bact_active_cell_ct']

subtitle_text = ('Metadata from whole-core melt fractions. '
                 + 'Created from compiled project metadata using the script '
                 + '<a href="https://github.com/cmfrantz/rottenice">'
                 + 'MetadataWhiskerPlots.py</a>. Analysis done by C. Frantz, '
                 + 'June 2025.')



#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    # Generate sample list and marker parameters for each month set
    samples = []
    markerprops = []
    for monthlist in sets:
        samplelist, markers = RottenIceModules.genSamples_w_Markers(
            monthlist, fractions)
        samples.append(samplelist)
        markerprops.append(markers)
    
    # Import metadata file
    filename, directory, metadata = RottenIceModules.fileGet(
        'Select sample metadata (environmental measurements) file',
        tabletype = 'metadata-sample')
    
    # Trim & prep metadata file
    # Flatten lists of lists
    allvars = np.unique([var for varset in varlist for var in varset])
    allsamples = np.unique([sample for sset in samples for sample in sset])
    metadata = metadata.loc[allsamples, allvars]
    # replace invalid values
    metadata = metadata.replace('na', np.nan)
    metadata = metadata.astype(float)
    
    # Set up figure
    rows = len(varlist)
    cols = max([len(list) for list in varlist])
    fig, axs = plt.subplots(rows, cols, figsize = (cols*pltw, rows*pltht))
    
    # Loop through metadata parameters to create each figure
    for row, varset in enumerate(varlist):
        for col in range(cols):
            # Build a blank plot if no variable for this column
            if col >= len(varset):
                axs[row,col].axis('off')
            else:
                # Grab data for each set of months
                var = varset[col]
                comb_data = []
                for n, monthset in enumerate(sets):
                    # Grab combined data for the set of months
                    # and add it as a list item in comb_data
                    comb_data.append(metadata.loc[samples[n], var])
                
                # Make whisker plot from non-nan data, hiding outliers
                boxplotvals = []
                for dataset in comb_data:
                    data = dataset[~np.isnan(dataset)]
                    boxplotvals.append(data)
                #axs[row,col].set_title(var)
                axs[row,col].boxplot(boxplotvals, showfliers = False)
                
                # Layer on individual points
                for n, dataset in enumerate(comb_data):
                    for p, datapoint in enumerate(dataset):
                        axs[row,col].plot([n+1], datapoint,
                                          **markerprops[n][p])
                        
                # Add axis labels
                '''
                xticklabels = [
                    '-'.join(months) + '\n(n=' + str(len(boxplotvals[i])) + ')'
                    for i,months in enumerate(sets)]
                '''
                xticklabels = [
                    '\n'.join([monthcode.get(m, m) for m in months])
                    for i, months in enumerate(sets)
                ]
                axs[row,col].set_xticklabels(xticklabels)
                axs[row,col].set_ylabel(RottenIceVars.metadataFullTitle[var])
                if var in vars2format:
                    axs[row,col].yaxis.set_major_formatter(
                        mtick.FormatStrFormatter('%.2e'))
    
    # Show and save the plot
    file_info = RottenIceVars.file_sets['metadata_boxplots']
    filename = file_info['pfx']
    page_title = file_info['title']
    html_filename = file_info['land_page']
    fig.tight_layout()
    fig.show()
    fig.savefig(directory + '\\' + filename + '.pdf', transparent = True)
    fig.savefig(directory + '\\' + filename + '.svg', transparent = True)
    
    # Generate HTML
    RottenIceModules.genHTMLfile(
        directory + '\\' + html_filename,
        page_title, subtitle_text,
        filename + '.svg',
        alt_text = 'Whisker plots for each metadata variable')
