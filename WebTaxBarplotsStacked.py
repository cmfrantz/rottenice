#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'

"""Rotten Ice Project: WebTaxBarplotsStacked
Created on Tue Jun 23 15:14:02 2020

@author: cariefrantz
@project: RottenIce

BUILDS INTERACTIVE (BOKEH) TAXONOMIC BARPLOTS: STACKED VERSION
This script builds an HTML webpage with interactive taxonomic bar plots.
It builds stacked plots displaying absolute and relative sequence abundances
for both DNA and cDNA, enabling direct comparisons of DNA vs. cDNA results.

Each taxonomic level is saved as a seperate webpage.

This version of the WebTaxBarplots script was used for sequencing data.
By modifying variables at the top of the code, it can easily be adapted for
any paired sets of OTU tables.

There is another version (WebTaxBarplotsSimple) that was used for
visualizing phytoplankton taxonomy data from microscopy counts.
It produces side-by-side plots and is best for a low number of samples.

This script was created as part of the Rotten Ice Project.


Arguments:  None

Requirements:      
    ASV tables (csv) for every sequence dataset analyzed
        where rows = ASVs, columns = samples
        header row and index column are specified in the variables below

Example in command line:
    python WebTaxBarplotsStacked.py

Dependencies Install:
    sudo apt-get install python3-pip python3-dev
    pip install tkinter
    pip install progress
    pip install numpy
    pip install pandas
    pip install matplotlib
    pip install math
    pip install bokeh

You will also need to have the following files \
    in the same directory as this script. \
They contain modules and variables that this script calls.
    RottenIceModules.py
    RottenIceVars.py
If you get an error indicating that one of these modules is not found, \
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
import os
from progress.bar import Bar

import numpy as np

from bokeh.io import output_file, show
from bokeh.layouts import row, column
from bokeh.models import FactorRange

import RottenIceModules
import RottenIceVars


####################
# VARIABLES
####################

# Set the max taxonomic depth to plot (e.g., deepest meaningful level)
max_level = 14

# Plot / HTML file title info
title = 'Sample taxonomy from amplicon sequencing'
subtitle_text = ('Rotten ice project Illumina sequencing data '
                 + 'processed by C. Frantz, May 2020 using pipeline-finished, '
                 + 'clean ASV tables using the script '
                 + '<a href="https://github.com/cmfrantz/rottenice">'
                 + 'WebTaxBarplotsStacked.py</a>.')
    
filename_prefix = RottenIceVars.file_sets['ASV_barplots']['pfx']

# List of major sample groups and 'others' for sample grouping by prefix
# Sample names with the prefixes on the left
#   will be given the category label on the right
groups = {
    'M'     : 'May',
    'JN'    : 'June',
    'JY10'  : 'July 10',
    'JY11'  : 'July 11',
    'Other' : 'Other'
    }
# Delimiter splitting the sample group prefix from the sample name
grpdelim = '-'   
# Sample names including any of these strings will be grouped in the 'Other'
# category and will be excluded from the focus taxa analysis
other_list = ['space', 'Blank', 'EL']

# List of taxonomic values to reassign,
# with their taxonomic level and new classification (after ':')
tax_reassign_list = {
    ('Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
     + 'Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
     + 'Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; '
     + 'Ambiguous_taxa; Ambiguous_taxa; D_14__')
            : 'Unassigned',
    ('Bacteria')
            : 'Bacteria; Other'}

# List of genes sequenced with maximum meaningful taxonomic level
# Each of these will be on its own page
genes       = {
    '16S'   : {
        'max_level'     : 7
        },
    '18S'   : {
        'max_level'     : 11
        }
    }

# Number of taxa plotted
max_taxa = 100   # Max number of taxa to show in a plot
f_cutoff = 0.1   # Taxa with relative abundance greater than this value
                 # in any sample in the dataset will be displayed
     
# Color
cmap = RottenIceVars.cmap

# Plot dimensions
plot_height = 220     # Height of the plot box
plot_width = 1500     # Width given to the actual plot
min_padding = 5  # Padding around plot
min_y_step = 5    # Number of y axis ticks (~5)
plot_bottom_pad = 100 # Padding given to x axis label
key_max_width = 925 # Max width needed for key
key_val_ht = 24   # Height taken up by each key entry
key_top_pad = 40  # Height taken up by key title

# Bokeh plot formatters
# Values to the right of the ':' shown without quotes are based on
# plot formatting variables defined above. Those in quotes can be modified.
shared_kwargs = {
    'x_axis_label'      : 'Sample',
    'toolbar_location'  : 'right',
    'tooltips'          : '@samples $name: @$name'
    }

# Plot details
# The types of plots to generate for each dataset
plttype = {
    'abs' : {
        'title'     : 'absolute',
        'data'      : 'pltdata',
        'xaxis'     : True,
        'ymax'      : 'auto',
        'widgets'   : 'hover,box_zoom,reset',
        'legend'    : False
        },
    'rel' : {
        'title'     : 'relative',
        'data'      : 'normpltdata',
        'xaxis'     : False,
        'ymax'      : 1,
        'widgets'   : 'hover,box_zoom,wheel_zoom,pan,reset',
        'legend'    : False
        }
    }

# All of the plots to draw in order and what goes in them
pltlist = {
    1 : {
        'template'  : 'DNA',
        'type'      : 'abs',
        'xlabelloc' : 'above',
        'ylabel'    : 'DNA ASV count \u27F6',
        'invert'    : False
        },
    2 : {
        'template'  : 'DNA',
        'type'      : 'rel',
        'xlabelloc' : 'below',
        'ylabel'    : 'DNA Relative abundance \u27F6',
        'invert'    : False
        },
    3 : {
        'template'  : 'cDNA',
        'type'      : 'rel',
        'xlabelloc' : 'below',
        'ylabel'    : '\u27F5 cDNA Relative abundance',
        'invert'    : True
        },
    4 : {
        'template'  : 'cDNA',
        'type'      : 'abs',
        'xlabelloc' : 'below',
        'ylabel'    : '\u27F5 cDNA ASV count',
        'invert'    : True
        }
    }


# PROGRAM VARIABLES
# (Modifying these may cause the program to break... proceed with caution)
# Data types
templates   = ['DNA', 'cDNA']


####################
# FUNCTIONS
####################

def loadFiles():
    '''Loads user-selected files for each plot dataset'''
    # Load the ASV data files
    directory = os.getcwd()
    datasets = {}
    
    # Loop through each gene to grab files for each gene
    for gene in genes:
        # Loop through each template type to grab files for each template
        for template in templates:
            dsetID = gene + '_' + template
            datasets[dsetID]={
                'gene'      : gene,
                'template'  : template,
                'title'     : gene + ' ' + template,
                'max_level' : genes[gene]['max_level']
                }
            filename, directory, data = RottenIceModules.fileGet(
                'Select ' + datasets[dsetID]['title'] + ' OTU Table CSV',
                directory = directory, tabletype = 'OTU_Table')
            datasets[dsetID]['data'] = data
            
    return datasets


def formatData(datasets):
    '''Formats datasets for use'''
    print('\n***************************\n'
          'Loading and formatting data\n'
          '***************************')   
    for ds in datasets:
        print(ds)
        data, samples = RottenIceModules.formatOTUtableData(
            datasets[ds]['data'], max_level = datasets[ds]['max_level'],
            tax_reassign_list = tax_reassign_list)
        grouped_samples = groupSamples(samples)
        mapper = dict(zip(samples, grouped_samples))
        datasets[ds]['data'] = data.rename(columns = mapper)
        datasets[ds]['samples'] = grouped_samples
    return datasets    


def groupSamples(samples):
    '''Group samples as outlined in the initial variables'''
    sampleset=[]
    for sample in samples:
        if any(other in sample for other in other_list):
            sampleset.append(('Other', sample))
        else:
            ssplit = sample.split('-')
            sampleset.append((groups[ssplit[0]], '-'.join(ssplit[-2:])))  
    return sampleset


def genFilenames():
    '''Generates list of output filename prefixes and the html header text'''
    filename_sets = {}
    for gene in genes:
         # Create filenames and html header navigation info for each gene    
        filepfx = filename_prefix + '_' + gene
        filenames = RottenIceModules.genFilenamesByLevel(
            filepfx, genes[gene]['max_level'], cmap = cmap)
        genes[gene]['filenames'] = filenames
        filename_sets[gene]=filenames
    nav_html = RottenIceModules.genPlotNavHTML_Sets(filename_sets)
    return filename_sets, nav_html
    

def updateKwargs():
    '''Updates bokeh plot kwargs based on variables defined above'''
    for pltnum in pltlist:
        pltlist[pltnum]['kwargs'] = {
            'plot_height'       : plot_height + min_padding * 2,
            'plot_width'        : plot_width + min_padding *2,
            'min_border_top'    : min_padding,
            'min_border_bottom' : min_padding,
            'min_border_right'  : min_padding,
            'min_border_left'   : min_padding,
            'x_axis_location'   : pltlist[pltnum]['xlabelloc'],
            'y_axis_label'      : pltlist[pltnum]['ylabel'],
            'tools'             : plttype[pltlist[pltnum]['type']]['widgets']
            }
        pltlist[pltnum]['kwargs'].update(shared_kwargs)


def samplesExclOther(samples):
    '''Remove 'others' from sample list containing multiple header levels'''
    samples_excl_other = []
    for sample in samples:
        if not any(other in sample[1] for other in other_list):
            samples_excl_other.append(sample)
    return samples_excl_other


def groupByAbundantTaxa(datasets, dslist):
    '''Group data according to most abundant 'focus' taxa'''
    dsets = []
    samples = []
    for ds in dslist:
        dsets.append(datasets[ds]['data'])
        samples.append(samplesExclOther(datasets[ds]['samples']))
    
    # Find abundant taxa
    focus_taxa = RottenIceModules.findAbundantTaxa(
        dsets, samples, fcutoff = f_cutoff, max_taxa = max_taxa)
    
    # Loop through each dataset and group taxa not in focus taxa
    for ds in dslist:
        # Rename all members of taxonomic groups that are in the same
        # L-1 taxonomic group as a focus taxon 'Other'
        print('Grouping ' + ds + ' ASVs at each taxonomic level...')
        grouped_data, others = RottenIceModules.groupTaxa(
            datasets[ds]['data'], focus_taxa)
        datasets[ds]['grouped_data'] = grouped_data
        datasets[ds]['excluded_others'] = others
        
        # Get max values
        samples = datasets[ds]['samples']
        sums = datasets[ds]['data'][samples].sum(axis = 0)
        max_count = np.max(sums)
        datasets[ds]['y_max'] = RottenIceModules.findMaxAxVal(max_count, 5)
        
    return datasets


def buildPlots(gene, datasets):
    '''Builds all of the plots at a given level'''
    # Build each plot in the list
    plots = []
    bar = Bar('', max = len(pltlist))
    for pltnum in pltlist: 
        
        # Retrieve the corresponding dataset
        dataset = datasets[gene + '_' + pltlist[pltnum]['template']]
        if pltlist[pltnum]['type'] == 'abs':
            data = dataset['pltdata']
            y_max = dataset['y_max']
        elif pltlist[pltnum]['type'] == 'rel':
            data = dataset['normpltdata']
            y_max = 1
            
        # Link x axes in the plot so all plots are on same x scale
        if not plots:
            sampleset = dataset['samples']
            xrange = FactorRange(*sampleset)
        else: xrange = plots[0].x_range
        
        # Build the plot
        p = RottenIceModules.buildBarPlot(
            data, colors, y_max, save_image = True,
            filename = (fname + '_'
                        + pltlist[pltnum]['template'] + '_'
                        + pltlist[pltnum]['type']),
            invert = pltlist[pltnum]['invert'],
            show_x_ax = plttype[pltlist[pltnum]['type']]['xaxis'],
            x_range = xrange,
            legend = False,
            **pltlist[pltnum]['kwargs'])
        plots.append(p)
        bar.next()
    bar.finish()
    
    return plots

#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    
    # Set up variables from variables defined in VARIABLE header
    # Generate output file filenames and HTML headers
    filename_sets, nav_html = genFilenames()
    updateKwargs()
    
   # Get and format ASV files from user
    datasets = loadFiles()
    datasets = formatData(datasets)
        
    # Build plot sets for each gene
    for gene in list(genes):
        print('\n***************************\n'
              'Processing ' + gene + ' data\n'
              '***************************')
        # Get list of datasets for the gene
        dslist = [dset for dset in datasets if datasets[dset]['gene']==gene]
        
        # Group data according to abundant taxa
        datasets = groupByAbundantTaxa(datasets, dslist)
        
        # Get the color map for the plot
        max_level = genes[gene]['max_level']
        tax_columns = RottenIceModules.levelCols(max_level) + ['taxonomy']
        colormap = RottenIceModules.colorMap2Tax(
            datasets[dslist[0]]['grouped_data'][tax_columns].copy(),
            max_level, cmap)
        
        # Build the plots
        for level in range(1,max_level+1):
            print('Generating Level ' + str(level) + ' plots...')
            
            # Prepare the data
            for ds in dslist:
                [
                    datasets[ds]['pltdata'],
                    datasets[ds]['normpltdata'],
                    colors
                    ] = RottenIceModules.buildPlotDicts(
                        datasets[ds]['grouped_data'].copy(),
                        level, datasets[ds]['samples'], colormap)
            
            # Set up HTML file
            fname = genes[gene]['filenames'][level-1]
            output_file(fname + '.html',
                        title = (datasets[ds]['title']
                                 + ' ASV barplots L' + str(level)))

            # Build plots and html files
            plots = buildPlots(gene, datasets)
            taxlist = list(datasets[ds]['pltdata'])
            # Build the legend
            taxlist.remove('samples')
            legend = RottenIceModules.genLegendOutside(
                taxlist, colors, key_val_ht, key_top_pad)
        
            # Build the HTML file
            grid = row([column(plots), legend])
            # Save and show HTML file
            show(grid)
        
            # Add header to page
            RottenIceModules.addHTMLhead(
                fname +'.html',
                page_title = title,
                page_nav_html = nav_html,
                subtitle_text = ('<b>Taxonomic level ' + str(level)
                                 + '</b><br />' + subtitle_text))