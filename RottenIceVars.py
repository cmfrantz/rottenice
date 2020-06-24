# -*- coding: utf-8 -*-
"""RottenIceVars
Created on Fri Jun 19 11:01:36 2020
@author: cariefrantz@author: cariefrantz
@project: RottenIce

VARIABLES SHARED BETWEEN DIFFERENT SCRIPTS IN THE ROTTEN ICE PROJECT
This document contains a collection of variables that are used in multiple \
    scripts in the Rotten Ice project.

Arguments:  None

Example in external script (must be either in the PATH or same directory):
    import RottenIceVars

Dependencies Install: None


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

# Colors
cmap = 'viridis'    # Change this to change the default colormap
cmaptext = ''''Choose a colormap, e.g.,
viridis, plasma, prism, nipy_spectral, jet, gist_earth, etc.
See https://matplotlib.org/tutorials/colors/colormaps.html for a full list. 
> '''

# Plot formatters
plotColorsByMonth = {      # Colors based on the viridis palette
	'M' :      '#441a54',  # Purple
	'JN':      '#3c538c',  # Blue
	'JY10':    '#abd03a',  # Light green
	'JY11':    '#fce61f'   # Yellow
    }

plotMarkersByFraction = {
	'HT' : '^',
	'HM' : 'o',
	'HB' : 'v'
    }

plotMarkerLineProperties = {
    'markeredgewidth'   : 0.5,
    'markeredgecolor'   : 'k'
    }

plotMarkerBorderByMonth = {
    'M':    True,
    'JN':   True,
    'JY10': True,
    'JY11': True
    }


# Data file format info
# Metadata file
metadata_head_row = 2       # Row in table containing unique variable names
metadata_index_col = 0      # Row in table containing unique sample names
metadata_filetype = 'csv'
# OTU (ASV) tables
OTU_table_head_row = 0      # Row in table containing unique sample names
OTU_table_index_col = 0     # Columnin table containing unique OTU IDs
OTU_table_filetype = 'csv'


# Metadata variable full titles
metadataFullTitle = {
    'day_num'       : 'Calendar day',
    'coord_lat'     : 'Latitude',
    'coord_lon'     : 'Longitude',
    'pH'            : 'pH',
    'temperature'   : 'Temperature (°C)',
    'bulk_density'  : 'Bulk ice density (' + r'$g \cdot cm^{3}$' + ')',
    'salinity_insitu' : 'Salinity (from 5 cm pucks; ppt)',
    'salinity'      : 'Salinity (ppt)',
    'SPM'           : 'SPM (' + r'$mg \cdot ml^{-1}$' + ')',
    'DOC'           : 'DOC (' + r'$mg C \cdot l^{-1}$' + ')',
    'Chl'           : 'Chlorophyll (' + r'$mg \cdot m^{3}$' + ')',
    'Phaeo'         : 'Phaeopigments (' + r'$mg \cdot m^{3}$' + ')',
    'Chl_Phaeo'     : '[Chl]/[Chl+Phaeo]',
    'FoFa'          : r'$F_o/F_a$',
    'PAM'           : 'PAM read',
    'SedLoad'       : 'Visible sediment load (1 low, 2 med, 3 high)',
    'bact_cell_ct'  : 'Bacteria (' + r'$cells \cdot ml^{-1}$' + ')',
    'CTC'           : 'CTC %',
    'bact_active_cell_ct' :  'active cells (' + r'$cells \cdot ml^{-1}$' + ')',
    'diatom_ct'     : 'Diatoms (' + r'$cells \cdot ml^{-1}$' + ')',
    'phyto_ct_select' : ('Phytoplankton - excluding ciliates & others)'
                         + ' \n(' + r'$cells \cdot ml^{-1}$' + ')'),
    'phyto_ct_all'  : 'Phytoplankton (' + r'$cells \cdot ml^{-1}$' + ')',
    'phyto_ct_other' : ('Phytoplankton - excluding diatoms)'
                        + '\n(' + r'$cells \cdot ml^{-1}$' + ')'),
    'pEPS'          : 'pEPS (' + r'$\mu g glucose \cdot ml^{-1}$' + ')',
    'POC'           : 'POC (' + r'$\mu g \cdot ml^{-1}$' + ')',
    'nitrogen'      : 'PN (' + r'$\mu g \cdot ml^{-1}$' + ')',
    'CN'            : 'C/N',
    'sterivex_vol_filtered' : 'Sterivex volume filtered (ml)',
    'Chl_pEPS'      : 'Ratio of Chl to pEPS',
    'CTC_pEPS'      : 'Ratio of %CTC to pEPS'
    }


# HTML Header
nav_html_start = '<b>Plot navigation: </b>'
html_head = '''
<h1>Rotten Ice Bio Data Analysis Navigation</h1>
<b>Metadata whisker plots: </b>
<a href = "http://faculty.weber.edu/cariefrantz/RottenIce
/metadata_whiskerplots.html">Link</a><br />
<b>Algae ID bar plots: </b>
<a href = "http://faculty.weber.edu/cariefrantz/RottenIce/
AlgaeID_barplots_L15_viridis.html">Link</a><br />
<b>ASV bar plots: </b>
<a href = "http://faculty.weber.edu/cariefrantz/RottenIce/
ASV-barplots/16S_barplots_L1_viridis.html">Link</a><br />
<b>Beta diversity clustered heatmaps: </b>
<a href = "http://faculty.weber.edu/cariefrantz/RottenIce/
B-div-heatmaps/heatmaps_16S_all.html">Link</a><br />
<b>Metadata vs. taxonomic group Spearman correlation analysis: </b>
<a href = "http://faculty.weber.edu/cariefrantz/RottenIce/
spearman/spearman_16S-DNA.html">Link</a><br />
<b>Python code: </b>
<a href="https://github.com/cmfrantz/rottenice">Link</a><br />
'''