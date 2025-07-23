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
    Sample Metadata (measurements) table (csv)
        Where rows = samples, columns = metadata characteristics
        Header row and index column are specified in the variables below
        This is NOT the same metadata table used in QIIME2, which provides
        metadata for each sample-gene-template sequenced. It is the master
        metadata table of environmental measurements for each sample.

Example in command line:
    python MetadataPCA.py

Dependencies Install:
    pip install numpy

You will also need to have the following files
    in the same directory as this script.
They contain modules and variables that this script calls.
    RottenIceModules.py
    RottenIceVars.py
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
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import RottenIceModules
import RottenIceVars

####################
# VARIABLES
####################

# Colors for coding the different months / floes
color_col = 'loc_id'    # Defines the column in the metadata table that contains the information mapped to colors below
color_map = {      # Pulls color hex values based on the standard color set
	'M-CS'     : RottenIceVars.plotColorsByMonth['M'],
	'JN-CS'    : RottenIceVars.plotColorsByMonth['JN'],
	'JY10'     : RottenIceVars.plotColorsByMonth['JY10'],
	'JY11'     : RottenIceVars.plotColorsByMonth['JY11']
    }

# Markers for coding the different horizons
marker_col = 'fraction' # Defines the column in the metadata table that contains the information mapped to markers below
marker_map = RottenIceVars.plotMarkersByFraction

# Sets of variables to plot, each is a different PCA analysis
varsets = {
    'physical'      : {
        'title'     : 'Physical parameters',
        'varlist'   : ['temperature','salinity','bulk_ice_density'],
        'arrows'    : True
        },
    'habitat'       : {
        'title'     : 'Habitat parameters',
        'varlist'   : ['temperature','salinity','bulk_ice_density',
                       'SPM','SedLoad','DOC','POC','nitrogen','pEPS'],
        'arrows'    : True
        },
    'nutrients'     : {
        'title'     : 'Nutrients',
        'varlist'   : ['DOC','POC','pEPS','nitrogen','CN'],
        'arrows'    : True
        },
    'carbon'        : {
        'title'     : 'Carbon fractions',
        'varlist'   : ['DOC','POC','pEPS','gels_total'],
        'arrows'    : True
        },
    'cells'         : {
        'title'     : 'Cell counts and activity',
        'varlist'   : ['cell_ct','CTC','phyto_ct_all','diatom_ct'],
        'arrows'    : True
        },
    'photosynthesis' : {
        'title'     : 'Measures of photosynthesis',
        'varlist'   : ['chl','phaeo','FoFa','PAM',
                       'diatom_ct','phyto_ct_all','phyto_ct_select',
                       'phyto_ct_other'],
        'arrows'    : True
        },
    'all-numeric'   : {
        'title'     : 'All numeric metadata',
        'varlist'   : ['date_num','lat','lon','depth_in_ice',
                       'salinity','SPM','DOC','chl','phaeo','FoFa',
                       'PAM','SedLoad','cell_ct','CTC','pEPS','POC',
                       'nitrogen','CN','temperature',
                       'salinity_bulk_phys','bulk_ice_density',
                       'diatom_ct','phyto_ct_all','phyto_ct_select',
                       'phyto_ct_other','gels_total'],
        'arrows'    : 10
        },
    'all-numeric-noNaN' : {
        'title'     : ('All numeric metadata for which all analyzed samples '
                       + 'have measured values (no NaN values)'),
        'varlist'   : ['date_num','lat','lon','depth_in_ice',
                       'salinity','DOC','chl','phaeo','FoFa',
                       'PAM','cell_ct','CTC','pEPS','POC',
                       'nitrogen','CN',
                       'diatom_ct','phyto_ct_all','phyto_ct_select',
                       'phyto_ct_other'],
        'arrows'    : 10
        },
    'interesting' : {
        'title'     : ('Metadata parameters of particular interest'),
        'varlist'   : ['date_num','salinity',
                       'DOC','chl','phaeo','FoFa','PAM','cell_ct','CTC',
                       'pEPS','POC','nitrogen','CN',
                       'diatom_ct','phyto_ct_other'],
        'arrows'    : True
        }
    }

# Info about file naming
file_info = RottenIceVars.file_sets['metadata_biplots']
out_filename = file_info['pfx']
title = file_info['title']

# HTML code for the for the HTML files
# Subheading text (Part II after the level number)
subtitle_text = '''
<p>Principal Correlation Analysis biplot showing PCA of metadata parameters
(markers) along with how different metadata parameters influence PCA space
(vectors). PCA calculations performed using the 
<a href="https://scikit-learn.org/">scikit-learn package</a> (1.7) for python. 
Created using the script <a href="https://github.com/cmfrantz/rottenice">
MetadataPCA.py</a>. Analysis done by C. Frantz, July 2025.</p>
'''

#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    
    # Import Metadata table
    filename, directory, metadata, = RottenIceModules.fileGet(
        'Select metadata file (from qiime)',
        tabletype = 'metadata-qiime')
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
    
    # Compute and save metadata correlations
    numeric_df = meta_trimmed.select_dtypes(include='number')
    corr_matrix = numeric_df.corr(method='spearman')
    corr_matrix.to_csv(
        directory + '\\' + out_filename + '_correlation-matrix.csv')
    plt.figure(figsize=(12,10))
    sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap = RottenIceVars.cmap,
                center=0)
    plt.title('Metadata Spearman Correlation Matrix')
    plt.tight_layout()
    plt.show()
    plt.savefig(directory + '\\' + out_filename + '_correlation-heatmap.svg')
    
    # Set up file navigation HTML
    filenames = [(varset, out_filename + '_' + varset) for varset in varsets]
    nav_html = RottenIceVars.nav_html_start
    for varset in filenames:
        nav_html += f' <a href="{varset[1]}.html">{varset[0]}</a> /'
    nav_html = nav_html[:-2]
    
    # Plot Set 1
    # Build the biplots for each set of variables defined
    for varset in varsets:
        print('Plotting ' + varset + 'PCA biplot')
        variables = varsets[varset]['varlist']
        
        # Run PCA & Biplot script
        scores, loadings, biplot_fig, arrow_key = RottenIceModules.biplot(
            meta_trimmed, samples, variables,
            'Metadata PCA Biplot: ' + varsets[varset]['title'],
            color_col, color_map, marker_col, marker_map,
            arrows = varsets[varset]['arrows'], marker_size = 100)
        
        # Save the biplot figure
        biplot_fig.savefig(
            directory + '\\' + out_filename + '_' + varset + '.png',
            transparent = True, bbox_inches='tight')
        biplot_fig.savefig(
            directory + '\\' + out_filename + '_' + varset + '.pdf',
            format = 'pdf')
   
        # Generate HTML per dataset
        html_path = f"{directory}\\{out_filename}_{varset}.html"
        html_title = file_info['title']
        html_body = subtitle_text
        vartext = ''
        for var in variables:
            vartext = vartext + (
                '<li>' + var + ': '
                + RottenIceModules.latex_to_html(
                    RottenIceVars.metadataFullTitle[var])
                + '</li>')
        html_body = (html_body
                     + '<p><b>Variables included in this analysis:</b><ul>'
                     +  vartext + '</ul></p>')

        RottenIceModules.genHTMLfile(
            html_path,
            html_title,
            html_body,
            image_filepaths=[out_filename + '_' + varset + '.png'],
            alt_text=['PCA Biplot: ' + varset],
            page_nav_html=nav_html
        )

        print(f'Saved HTML and images for {varset} to {directory}')