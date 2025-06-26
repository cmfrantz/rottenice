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

# Datasets
templates = ['DNA', 'cDNA']
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
other_samples = ['Blank', 'EL']

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
data_table_fmts = {
    'metadata'  : {
        'head_row'  : 2,
        'index_col' : 0,
        'filetype'  : 'csv'
        },
    'OTU-table' : {
        'head_row'  : 0,
        'index_col' : 0,
        'filetype'  : 'csv'
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
    'coord_lat'     : 'Latitude',
    'coord_lon'     : 'Longitude',
    'pH'            : 'pH',
    'temperature'   : 'Temperature (' + r'$\degree$' + 'C)',
    'bulk_density'  : ('Bulk ice density (g ' + r'$\cdot$'
                       + ' cm' + r'$^{3}$' + ')'),
    'salinity_direct' : 'Salinity (from 5 cm puck direct melts; ppt)',
    'salinity_lab'  : ('Salinity (from lab melts; whole-horizon melts are '
                       + 'mixed with ASW; ppt)'),
    'salinity_comb' : 'Salinity (ppt)',
    'salinity'      : 'Salinity (ppt)',
    'SPM'           : 'SPM (mg ' + r'$\cdot$' + ' ml' + r'$^{-1}$' + ')',
    'DOC'           : 'DOC (mg C ' + r'$\cdot$' + ' l' + r'$^{-1}$' + ')',
    'Chl'           : 'Chlorophyll (mg ' + r'$\cdot$' + ' m' + r'$^{3}$' + ')',
    'Phaeo'         : ('Phaeopigments (mg ' + r'$\cdot$'
                       + ' m' + r'$^{3}$' + ')'),
    'Chl_Phaeo'     : '[Chl]/[Chl+Phaeo]',
    'FoFa'          : r'$F_o/F_a$',
    'PAM'           : 'PAM read',
    'SedLoad'       : 'Visible sediment load (1 low, 2 med, 3 high)',
    'bact_cell_ct'  : ('Bacteria (cells ' + r'$\cdot$' + ' ml'
                       + r'$^{-1}$' + ')'),
    'CTC'           : 'CTC %',
    'bact_active_cell_ct' :  ('Active cells (cells ' + r'$\cdot$' + ' ml'
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
        'land_page' : 'chem_heatmaps.html'},
    'metadata_boxplots'     : {
        'title'     : 'Metadata box plots',
        'pfx'       : 'metadata_whiskerplots',
        'land_page' : 'metadata_whiskerplots.html'},
    'algae_barplots'        : {
        'title'     : 'Phytoplankton ID bar plots',
        'pfx'       : 'PhytoID_barplots',
        'land_page' : 'Algae/PhytoID_barplots_L1_viridis.html'},
    'ASV_barplots'          : {
        'title'     : 'Sequencing ASV bar plots',
        'pfx'       : 'ASV_barplots',
        'land_page' : 'ASV-barplots/ASV_barplots_16S_L1_viridis.html'
        },
    'alpha_metrics'         : {
        'title'     : 'Alpha diversity metrics',
        'pfx'       : 'alpha_boxplots',
        'land_page' : 'alpha-div/alpha_boxplots_16S_DNA.html'},
    'beta_cluster'          : {
        'title'     : 'Beta diversity sample cluster heatmaps',
        'pfx'       : 'beta_clustermap',
        'land_page' : 'B-div-heatmaps/beta_clustermap_16S_all.html'},
    'spearman_taxonomy'     : {
        'title'     : ('Metadata vs. taxonomic group '
                       + 'Spearman correlation analysis'),
        'pfx'       : 'spearman',
        'land_page' : 'spearman/spearman_16S-DNA.html'
        },
    'spearman_diversity'    : {
        'title'     : ('Metadata vs. diversity '
                       + 'Spearman correlation analysis'),
        'pfx'       : 'spearman_div',
        'land_page' : 'spearman/spearman_div.html'
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