# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 10:24:04 2025

@author: cariefrantz
@project: RottenIce

BUILDS PCA BIPLOTS OF WHOLE-HORIZON MELT METADATA
This script calculates principal components and builds biplots for metadata
from the whole-horizon melt samples, with metadata contributing significantly
to PC space plotted as vectors.

This script was created as part of the Rotten Ice Project.


Arguments:  None

Requirements:      
    Metadata table (tsv) for the rotten ice samples.
        This is the same metadata table format as used in QIIME2
        where rows = sample names
        columns = metadata parameters

Example in command line:
    python MetadataPCA.py

Dependencies Install:
    pip install numpy

You will also need to have the following files
    in the same directory as this script.
They contain modules and variables that this script calls.
    RottenIceModules.py
If you get an error indicating that one of these modules is not found,
    change the working directory to the directory containing these files.

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
import numpy as np
import RottenIceModules

####################
# VARIABLES
####################

# Colors for coding the different months / floes
color_col = 'loc_id'    # Defines the column in the metadata table that contains the information mapped to colors below
color_map = {      # Colors based on the viridis palette
	'M-CS' :   '#441a54',  # Purple
	'JN-CS':   '#3c538c',  # Blue
	'JY10':    '#abd03a',  # Light green
	'JY11':    '#fce61f'   # Yellow
    }

# Markers for coding the different horizons
marker_col = 'fraction' # Defines the column in the metadata table that contains the information mapped to markers below
marker_map = {
    'HT'    : '^',
    'HM'    : 's',
    'HB'    : 'v'
    }

# Sets of variables to plot, each is a different PCA analysis
varsets = {
    'physical parameters'   : ['temperature','salinity','bulk_ice_density'],
    'habitat'               : ['temperature','salinity','bulk_ice_density',
                               'SPM','SedLoad','DOC','POC','nitrogen','pEPS'],
    'nutrients'             : ['DOC','POC','pEPS','nitrogen','CN'],
    'carbon'                : ['DOC','POC','pEPS','gels_total'],
    'cells'                 : ['cell_ct','CTC','diatom_ct'],
    'photosynthesis'        : ['chl','phaeo','FoFa','PAM',
                               'diatom_ct','phyto_ct_all','phyto_ct_select',
                               'phyto_ct_other'],
    'all_numeric'           : ['date_num','lat','lon','depth_in_ice',
                               'salinity','SPM','DOC','chl','phaeo','FoFa',
                               'PAM','SedLoad','cell_ct','CTC','pEPS','POC',
                               'nitrogen','CN','temperature',
                               'salinity_bulk_phys','bulk_ice_density',
                               'diatom_ct','phyto_ct_all','phyto_ct_select',
                               'phyto_ct_other','gels_total'],
    'all_numeric_noNaN'     :  ['date_num','lat','lon','depth_in_ice',
                               'salinity','DOC','chl','phaeo','FoFa',
                               'PAM','cell_ct','CTC','pEPS','POC',
                               'nitrogen','CN',
                               'diatom_ct','phyto_ct_all','phyto_ct_select',
                               'phyto_ct_other']
    }


#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    
    # Import Metadata table
    filename, directory, metadata, = RottenIceModules.fileGet(
        'Select metadata table', tabletype = 'metadata')
    metadata = metadata.dropna(how = 'all')
    metadata = metadata.replace('na', np.nan)
    
    # Trim the metadata table to only the samples of interest,
    # as defined in the color and marker maps
    # Since the metadata table duplicates values for the two genes (16S/18S),
    #  templates (DNA/cDNA), and replicates, trim to just 16S DNA values
    meta_trimmed = metadata[
        (metadata[color_col].isin(list(color_map))) &
        (metadata[marker_col].isin(list(marker_map))) &
        (metadata['gene'] == '16S') &
        (metadata['template'] == 'DNA') &
        (metadata['replicate'] == 1)
        ]
    samples = list(meta_trimmed.index)
    
    for varset in varsets:
        print('Plotting ' + varset + 'PCA biplot')
        variables = varsets[varset]
        
        # Run PCA & Biplot script
        scores, loadings, biplot_fig = RottenIceModules.biplot(
            meta_trimmed, samples, variables, 'Metadata PCA Biplot: ' + varset,
            color_col, color_map, marker_col, marker_map,
            n_arrows = 10, marker_size = 100)
        
        # Save the biplot figure
        biplot_fig.savefig(
            directory + '\\' + 'Metadata_Biplot ' + varset + '.png', transparent = True)
        biplot_fig.savefig(
            directory + '\\' + 'Metadata_Biplot ' + varset + '.pdf', format = 'pdf')
    
