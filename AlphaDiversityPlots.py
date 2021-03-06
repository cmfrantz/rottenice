#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'
"""Rotten Ice Project: DiversityWhiskerPlots

Created on Thu Jul  9 17:17:31 2020
@author: cariefrantz
@project: RottenIce

CREATES WHISKER PLOTS COMPARING DIVERSITY METRICS
This script creates whisker box plots for alpha diversity metrics calculated
for different samples: Shannon's Index, Faith's Phylogenetic Diversity, and
Pielou's Evenness.
This script was created as part of the Rotten Ice Project.

Arguments:  None

Requirements:
    16S Alpha diversity metrics table (csv)
    18S Alpha diversity metrics table (csv)

Example in command line:
    python AlphaDiversityPlots.py
    
Dependencies Install:
    sudo apt-get install python3-pip python3-dev
    pip install numpy
    pip install matplotlib
    
You will also need to have the following files in the same directory as this
script. They contain modules and variables that this script calls.
    RottenIceModules.py
    RottenIceVars.py
If you get an error indicating that one of these modules is not found, change
the working directory to the directory containing these files.

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

import matplotlib
import matplotlib.pyplot as plt
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import RottenIceModules
import RottenIceVars


####################
# VARIABLES
####################

# Variables defining what samples and metadata parameters to plot
genes = RottenIceVars.genes
templates = RottenIceVars.templates
metrics = RottenIceVars.alpha_diversity_metrics

# Variables defining plot parameters
# Plot size
pltht = 5
pltw = 10

page_title = RottenIceVars.file_sets['alpha_metrics']['title']
file_pfx = RottenIceVars.file_sets['alpha_metrics']['pfx']
alt_text_pfx = 'Alpha diversity metric boxplots: ' 

subtitle_text = ('Diversity metrics for pipeline-finished, rarefied sequencing'
                 'results of Rotten Ice Illumina sequence data. '
                 + 'Created from compiled alpha diversity metrics using the '
                 + 'script <a href="https://github.com/cmfrantz/rottenice">'
                 + 'AlphaDiversityPlots.py</a>. Analysis done by C. Frantz, '
                 + 'July 2020.')



####################
# FUNCTIONS
####################

def buildWhiskerPlot(ax, data, samples, metric, title=''):
    '''Builds whisker plot for the given metric'''
    print('Plotting ' + metric + '...')
    # Prep the data
    sample_data = []
    for sample in samples:
        sample_data.append(data.loc[data['sample']==sample,
                                    metric])
    # Plot the data
    ax.boxplot(sample_data)
    ax.set_title(title)
    ax.set_ylim(ymin = metrics[metric]['yrange'][0],
                ymax = metrics[metric]['yrange'][1])
    ax.set_xticklabels(samples,
                              rotation = 'vertical', ha = 'center')
    
#%%

####################
# MAIN FUNCTION
####################

if __name__ == '__main__':
    
    # Get list of samples
    sample_list = RottenIceModules.genSampleListCS()
    
    # Set up page navigation HTML
    page_nav_html = '<p>' + RottenIceVars.nav_html_start
    for gene in genes:
        page_nav_html = page_nav_html + '<em>' + gene + ':  </em>'
        for template in templates:
            page_nav_html = (page_nav_html 
                             + '<a href="' + file_pfx + '_'
                             + gene + '_' + template + '.html">' + template
                             + '</a> ')
    page_nav_html = page_nav_html + '</p>'
        
    # Get data from the user
    directory = os.getcwd()
    divdata = {}
    for gene in genes:
        filename, directory, divdata[gene] = RottenIceModules.fileGet(
        'Select ' + gene + ' alpha diversity metrics table',
        tabletype = 'alpha-div', directory = directory)
    
    # Set up master (combined) figure
    rows = len(genes) * len(metrics)
    cols = len(templates)
    plt.figure(1)
    m_fig, m_axs = plt.subplots(rows, cols,
                               figsize = (cols * pltw, rows * pltht),
                               sharey = 'row')
        
    # Make a seperate page for each gene and template
    for i, gene in enumerate(genes):
        for col, template in enumerate(templates):
            print('Building ' + gene + ' ' + template + ' page...')
            filename = file_pfx + '_' + gene + '_' + template
            
            # Grab the data
            data = divdata[gene].copy()
            temp_list = [sample for sample in data.index
                         if '-'+template in sample]
            data = data.loc[temp_list]
            # Get the metadata
            samplelist = [sample.split('.')[0] for sample in temp_list]
            data['sample'] = ['-'.join(sample.split('-')[:-1])
                                   for sample in samplelist]
            
            # Make whisker plots for each metric
            figlist = []
            alt_text = []
            for m, metric in enumerate(metrics):
                fig, ax = plt.subplots(figsize = (pltw,pltht))
                buildWhiskerPlot(
                    ax, data, sample_list, metric,
                    title = (gene + ' ' + template + ': ' 
                             + metrics[metric]['title']))

                # Save the figure
                fig.tight_layout()
                figname = filename + '_' + metric
                fig.savefig(directory + '\\' + figname + '.pdf',
                            transparent = True)
                fig.savefig(directory + '\\' + figname + '.svg',
                            transparent = True)
                figlist.append(figname + '.svg')
                alt_text.append(alt_text_pfx + metrics[metric]['title'])
                
                # Copy the plot to the combined figure
                plt.figure(1)
                m_axs[(i*len(metrics) + m)][col] = buildWhiskerPlot(
                    m_axs[(i*len(metrics) + m)][col],
                    data, sample_list, metric,
                    title = (gene + ' ' + template + ': '
                             + metrics[metric]['title']))
                
            # Make HTML page for each gene-template pair
            RottenIceModules.genHTMLfile(
                directory + '\\' + filename + '.html', page_title,
                ('<p><b>' + gene + ' ' + template + '</b><br />'
                 + subtitle_text + '</p>'),
                figlist, alt_text = alt_text, page_nav_html = page_nav_html)
            
    # Save the master (combined) figure
    print('Saving master figure...')
    plt.figure(1)
    m_fig.tight_layout()
    m_fig.savefig(directory + '\\' + file_pfx + '_combined.pdf',
                  transparent = True)
    m_fig.savefig(directory + '\\' + file_pfx + '_combined.svg',
                  transparent = True)