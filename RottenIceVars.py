# -*- coding: utf-8 -*-
"""RottenIceVars
Created on Fri Jun 19 11:01:36 2020
@author: cariefrantz@author: cariefrantz
@project: RottenIce

VARIABLES SHARED BETWEEN DIFFERENT SCRIPTS IN THE ROTTEN ICE PROJECT
This document contains a collection of variables that are used in multiple
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

#################
# SAMPLE METADATA

months = {
    'May'       : 'M',
    'June'      : 'JN',
    'July 10'   : 'JY10',
    'July 11'   : 'JY11'
    }

# Available samples collected each month for the CS dataset
fraction_sets_CS = {
    'May'       : ['HT', 'HM', 'HB',
                   'IT', 'IM', 'IB',
                   'BT', 'BM', 'BB',
                   'P1', 'P2', 'SW'],
    'June'      : ['HT', 'HM', 'HB',
                   'IT', 'IM', 'IB',
                   'BT', 'BM', 'BB',
                   'P1', 'P2', 'SW'],
    'July 10'   : ['HT', 'HM', 'HB',
                   'IT', 'IM', 'IB',
                   'B', 'Drain', 'PW'],
    'July 11'   : ['HT', 'HM', 'HB',
                   'IT', 'IM', 'IB',
                   'Drain', 'SW']
    }

# Indicators of 'outgroup' samples to exclude from some analyses
other_samples = ['Blank', 'Negative', 'Unknown', 'EL']


#################
# BIO ANALYSIS VARIABLES

templates = ['DNA', 'cDNA']

# Taxonomic level definitions in order
tax_levels = ['domain','phylum','class','order','family','genus','species']

# Genes sequenced in the project
genes = {
    '16S'   : {
        'max_level'     : 7,
        'tax_reassign_list'     : {
            ('Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
             + 'Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
             + 'Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
             + 'Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
             + 'Ambiguous_taxa; Ambiguous_taxa; D_14__')            :
                ('Unassigned'),
            ('D_0__Bacteria')                                       :
                ('Bacteria; Other')
            }
        },
    '18S'   : {
        'max_level'     : 11,
        'tax_reassign_list'     : {}
        }
    }

# Alpha diversity metrics calculated
alpha_diversity_metrics = {
    'pielou_e'  : {
        'title'     : "Pielou's evenness",
        'yrange'    : (0,1)
        },
    'shannon'   : {
        'title'     : "Shannon's diversity index",
        'yrange'    : (0,8)
        },
    'faith_pd'  : {
        'title'     : "Faith's phylogenetic diversity",
        'yrange'    : (0,50)
        }
    }
    

#################
# COLORS & MARKERS
cmap = 'viridis'    # Change this to change the default colormap
cmaptext = ''''Choose a colormap, e.g.,
viridis, plasma, prism, nipy_spectral, jet, gist_earth, etc.
See https://matplotlib.org/tutorials/colors/colormaps.html for a full list. 
> '''

# Categorical colors to identify the months
# Based on the viridis color palette
plotColorsByMonth = {      
	'M'    :      '#441a54',  # Purple
	'JN'   :      '#3c538c',  # Blue
	'JY10' :    '#abd03a',  # Light green
	'JY11' :    '#fce61f'   # Yellow
    }

# Categorical color palette to identify the fractions
'''
# Based on the 'iceberg1' (cool) palette:
# https://www.color-hex.com/color-palette/91134
plotColorsByFraction = {
    # Ice-only melts in tints of #4F518B
    'IT'    : [149,150,185],
    'IM'    : [114,115,162],
    'IB'    : [79,81,139],
    # Whole-core melts in tints of #51729c
    'HT'    : [150,170,195],
    'HM'    : [115,142,175],
    'HB'    : [81,114,156],
    # Brines in tints of #4b8ba7
    'BT'    : [147,185,202],
    'BM'    : [110,162,184],
    'BB'    : [75,139,167],
    'B'     : [56,104,125],
    # Sackhole percolates and drains in tints of #7fb8b1
    'P1'    : [178,212,208],
    'P2'    : [152,198,192],
    'Drain' : [127,184,177],
    # Other fluids
    'PW'    : [207,237,212],
    'SW'    : [176,225,184],
    'Blank' : [255,255,255]
    }

for f in fraction_cmap:
    fraction_cmap[f] = [c/256 for c in fraction_cmap[f]]
'''
# Based on the 'Seasons of the Dead - Faded Spectrum' (muted rainbow) palette:
# https://www.color-hex.com/color-palette/14353
plotColorsByFraction = {
    # Ice-only melts in tints of #5aa095
    'IT'    : '#7ab3aa',
    'IM'    : '#5aa095', # blue
    'IB'    : '#488077',
    # Whole-core melts in tints of #ddd481
    'HT'    : '#e3dc9a',
    'HM'    : '#ddd481', # yellow
    'HB'    : '#b0a967',
    # Brines in tints of #b25959
    'BT'    : '#c17a7a',
    'BM'    : '#b25959', # red
    'BB'    : '#8e4747',
    'B'     : '#6a3535',
    # Sackhole percolates and drains in tints of #dd9d72
    'P1'    : '#e3b08e',
    'P2'    : '#dd9d72', # orange
    'Drain' : '#b07d5b',
    # Other fluids in tints of #83ca82
    'PW'    : '#b4dfb4', # light green
    'SW'    : '#5b8d5b', # dark green
    'BW'    : '#83ca82', # green
    'Blank' : '#ffffff'
    }

# Base the material categorical colormap on the fraction colormap used
plotColorsByMaterial = {
    'I'     : plotColorsByFraction['IM'],
    'H'     : plotColorsByFraction['HM'],
    'B'     : plotColorsByFraction['BM'],
    'P'     : plotColorsByFraction['P1'],
    'D'     : plotColorsByFraction['Drain'],
    'PW'    : plotColorsByFraction['PW'],
    'SW'    : plotColorsByFraction['SW'],
    'BW'    : plotColorsByFraction['BW'],
    'Blank' : plotColorsByFraction['Blank']
    }

plotMarkersByHorizon = {
    'T' : '^', # Top = triangle up
    'M' : 's', # Mid = square
    'B' : 'v', # Bottom = triangle down
    'A' : 'h', # Whole-core = hexagon
    'W' : 'o'  # Waters = circle
    }

plotMarkersByFraction = {
	'HT' : plotMarkersByHorizon['T'],
	'HM' : plotMarkersByHorizon['M'],
	'HB' : plotMarkersByHorizon['B']
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
data_table_fmts = {
    'metadata-qiime'  : {   # This is the metadata file used in QIIME2
                            # For each sample - gene - template sequenced
        'head_row'  : 0,
        'index_col' : 0,
        'filetype'  : 'tsv'
        },
    'metadata-sample' : {   # This is the master sample measurement table csv
        'head_row'  : 2,
        'index_col' : 0,
        'filetype'  : 'csv'
        },
    'ASV-table' : {
        'head_row'  : 0,
        'index_col' : 0,
        'filetype'  : 'tsv'
        },
    'OTU-table' : {     # This is the format exported from QIIME2's barplots
        'head_row'  : 0,
        'index_col' : 0,
        'filetype'  : 'csv'
        },
    'ASV-table-fromBIOM' : {
        # This is the format using the series of commands to add taxonomy to
        # QIIME feature-table qza files using BIOM intermediates.
        # Unlike the barplot exports, it retains the feature IDs in addition
        # to the taxonomy.
        'head_row'  : 1,
        'index_col' : 0,
        'filetype'  : 'tsv'
        },
    'alpha-div' : {
        'head_row'  : 0,
        'index_col' : 0,
        'filetype'  : 'csv'
        }
    }


# Statistical preferences
p_cutoff = 0.05

# Metadata variable full titles
metadataFullTitle = {
    'day_num'       : 'Calendar day',
    'date_num'      : 'Calendar day',
    'coord_lat'     : 'Latitude',
    'coord_lon'     : 'Longitude',
    'lat'           : 'Latitude',
    'lon'           : 'Longitude',
    'depth_in_ice'  : 'Depth below the ice surface (cm)',
    'pH'            : 'pH',
    'temperature'   : 'Temperature (' + r'$\degree$' + 'C)',
    'bulk_ice_density'  : ('Bulk ice density (g ' + r'$\cdot$'
                       + ' cm' + r'$^{3}$' + ')'),
    'salinity_direct' : 'Salinity (from 5 cm puck direct melts; ppt)',
    'salinity_bulk_phys' : 'Salinity (from 5 cm puck direct melts; ppt)',
    'salinity_lab'  : ('Salinity (from lab melts; whole-horizon melts are '
                       + 'mixed with ASW; ppt)'),
    'salinity_comb' : 'Salinity (ppt)',
    'salinity'      : 'Salinity (ppt)',
    'SPM'           : 'SPM (mg ' + r'$\cdot$' + ' ml' + r'$^{-1}$' + ')',
    'DOC'           : 'DOC (mg C ' + r'$\cdot$' + ' l' + r'$^{-1}$' + ')',
    'chl'           : 'Chlorophyll (mg ' + r'$\cdot$' + ' m' + r'$^{3}$' + ')',
    'phaeo'         : ('Phaeopigments (mg ' + r'$\cdot$'
                       + ' m' + r'$^{3}$' + ')'),
    'Chl_Phaeo'     : '[Chl]/[Chl+Phaeo]',
    'FoFa'          : r'$F_o/F_a$',
    'PAM'           : 'PAM read' + r'$F_v/F_m$',
    'SedLoad'       : 'Visible sediment load (1 low, 2 med, 3 high)',
    'cell_ct'       : ('Bacteria (cells ' + r'$\cdot$' + ' ml'
                       + r'$^{-1}$' + ')'),
    'CTC'           : 'CTC %',
    'bact_active_cell_ct' :  ('Active cells (cells ' + r'$\cdot$' + ' ml'
                              + r'$^{-1}$' + ')'),
    'cells_act'     : ('Active cells (cells ' + r'$\cdot$' + ' ml'
                              + r'$^{-1}$' + ')'),
    'diatom_ct'     : ('Diatoms (cells ' + r'$\cdot$' + ' ml'
                       + r'$^{-1}$' + ')'),
    'phyto_ct_select' : ('Phytoplankton - excluding ciliates & others'
                         + '\n(cells ' + r'$\cdot$'
                         + ' ml' + r'$^{-1}$' + ')'),
    'phyto_ct_all'  : ('Phytoplankton (cells ' + r'$\cdot$'
                       + ' ml' + r'$^{-1}$' + ')'),
    'phyto_ct_other' : ('Phytoplankton - excluding diatoms'
                        + '\n(cells ' + r'$\cdot$' + ' ml' + r'$^{-1}$' + ')'),
    'pEPS'          : ('pEPS (' + r'$\mu$' + 'g glucose ' + r'$\cdot$'
                       + ' ml' + r'$^{-1}$' + ')'),
    'POC'           : ('POC (' + r'$\mu$' + 'g ' + r'$\cdot$'
                       + ' ml' + r'$^{-1}$' + ')'),
    'nitrogen'      : ('PN (' + r'$\mu$' + 'g ' + r'$\cdot$' + ' ml'
                       + r'$^{-1}$' + ')'),
    'CN'            : 'C/N',
    'sterivex_vol_filtered' : 'Sterivex volume filtered (ml)',
    'Chl_pEPS'      : 'Ratio of Chl to pEPS',
    'CTC_pEPS'      : 'Ratio of %CTC to pEPS',
    'gels_total'    : 'Gel count per mL',
    'gels_sm'       : 'Gel % < 0.5 ' + r'$\mu$' + 'm',
    'gels_md'       : 'Gel % < 1 ' + r'$\mu$' + 'm',
    'gels_lg'       : 'Gel % > 1 ' + r'$\mu$' + 'm',
    'gels_xl'       : 'Gel % > 10 ' + r'$\mu$' + 'm'    
    }

# Types of files and info about them
file_sets = {
    'chem_heatmaps'         : {
        'title'     : 'Metadata heatmaps',
        'pfx'       : 'chem_heatmaps',
        'land_page' : 'chem_heatmaps.html'
        },
    'metadata_boxplots'     : {
        'title'     : 'Metadata box plots',
        'pfx'       : 'metadata_whiskerplots',
        'land_page' : 'metadata_whiskerplots.html'
        },
    'metadata_biplots'      : {
        'title'     : 'Metadata PCA biplots',
        'pfx'       : 'metadata_biplot',
        'land_page' : 'metadata_biplots_physical.html'
        },
    'algae_barplots'        : {
        'title'     : 'Phytoplankton ID bar plots',
        'pfx'       : 'PhytoID_barplots',
        'land_page' : 'Algae/PhytoID_barplots_L1_viridis.html'
        },
    'ASV_barplots'          : {
        'title'     : 'Sequencing ASV bar plots',
        'pfx'       : 'ASV_barplots',
        'land_page' : 'ASV-barplots/ASV_barplots_16S_L1_viridis.html'
        },
    'alpha_metrics'         : {
        'title'     : 'Alpha diversity metrics',
        'pfx'       : 'alpha_boxplots',
        'land_page' : 'alpha-div/alpha_boxplots_16S_DNA.html'
        },
    'beta_cluster'          : {
        'title'     : 'Beta diversity sample cluster heatmaps',
        'pfx'       : 'beta_clustermap',
        'land_page' : 'B-div-heatmaps/beta_clustermap_16S_all.html'
        },
    'spearman_taxonomy'     : {
        'title'     : ('Metadata vs. taxonomic group '
                       + 'Spearman correlation analysis'),
        'pfx'       : 'spearman',
        'land_page' : 'spearman/spearman_16S-DNA.html'
        }
    }


    
# HTML Header for generating webpages
# Page header text
web_nav_header = '<h1>Rotten Ice Bio Data Analysis Navigation</h1>'
# Root directory where html files are stored
web_dir = 'http://faculty.weber.edu/cariefrantz/RottenIce/'
# Other stuff to be added to the navigation list
other_nav = ('<b>Python code: </b>'
             + '<a href="https://github.com/cmfrantz/rottenice">Link</a>')
# Text denoting local navigation
nav_html_start = '<b>Plot navigation: </b>'