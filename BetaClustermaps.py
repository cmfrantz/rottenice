#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Carie Frantz'
__email__ = 'cariefrantz@weber.edu'
"""Rotten Ice Project: BetaClustermaps
Created on Wed Jun 26 15:50:09 2019
@author: cariefrantz
@project: RottenIce

BUILDS SAMPLE SIMILARITY HEATMAPS WITH SAMPLE CLUSTERING
This script produces sample clustering heatmaps from sample distance matrices
based on sequencing data. It generates several sample groupings as 
webpages.

This script was created as part of the Rotten Ice Project.


Arguments:  None

Requirements:      
    ASV tables (csv) for every sequence dataset analyzed
        where rows = ASVs, columns = samples
        header row and index column are specified in the variables below

Example in command line:
    python BetaClustermaps.py

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
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

import matplotlib
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import RottenIceModules
import RottenIceVars


####################
# VARIABLES
####################

# Outputdirectory
# Define this if there is a directory where files should go.
# If None, the default is the working directory
outdir = None

# Plot / HTML file title info
title = RottenIceVars.file_sets['beta_cluster']['title']
subtitle_text = ('Heatmaps generated from weighted UniFrac distance matrices '
                 + 'calculated on quality-controlled Rotten Ice sample '
                 + 'amplicon sequences. Samples were grouped using '
                 + 'nearest-point clustering on correlation distances using '
                 + 'the <a href="https://seaborn.pydata.org/generated/'
                 + 'seaborn.clustermap.html">seaborn.clustermap</a> package '
                 + '(v. 0.12.2) for Python. Data processed using the script '
                 + '<a href="https://github.com/cmfrantz/rottenice">'
                 + 'BetaClustermaps.py</a>. C. Frantz, July 2025</p>')
file_pfx = RottenIceVars.file_sets['beta_cluster']['pfx']
alt_text_pfx = 'Sample clustermap: '

# Correct sample info when sample names are misleading
correct_sample_info = {
    'JN-SW' : 'BW'
    }

# Plot variables
cmap = RottenIceVars.cmap
genes = {
    '16S'   : {},
    '18S'   : {},
    '18S-PP' : {}
    }
templates = ['DNA', 'cDNA']
sets = {
    'all'       : {'title' : 'all'},
    'template'  : {'title' : 'by template'},
    'month'     : {'title' : 'by month'},
    'fraction'  : {'title' : 'by fraction'}
    }

# Categorical color palette to identify the months
# Based on the standard colors used in RottenIceVars
month_cmap = {
    'M'     : RottenIceVars.plotColorsByMonth['M'],
    'JN'    : RottenIceVars.plotColorsByMonth['JN'],
    'JY10'  : RottenIceVars.plotColorsByMonth['JY10'],
    'JY11'  : RottenIceVars.plotColorsByMonth['JY11'],
    'Blank' : '#FFFFFF'     # White
    }
 
# Categorical color palette to identify the fractions
# Based on the standard colors used in RottenIceVars
fraction_cmap = RottenIceVars.plotColorsByFraction

months = {
    'M'     : 'May',
    'JN'    : 'June',
    'JY10'  : 'July 10 "dirty" floe',
    'JY11'  : 'July 11 "clean" floe'
    }

fractions = {
    'I' : {
        'title'     : 'Ice-only melts',
        'subset'    : ['IT','IM','IB']
        },
    'H' : {
        'title'     : 'Whole-horizon melts',
        'subset'    : ['HT','HM','HB']
        },
    'F' : {
        'title'     : 'Fluids',
        'subset'    : ['BT', 'BM', 'BB', 'B','P1','P2','PW','BW','SW','Drain']
        }
    }
    
#%%

####################
# FUNCTIONS
####################

def parseSamples(distMatrix):
    '''Parses sample names into subcategories'''
    columns = list(distMatrix.columns)
    headerlist = []
    for sname in columns:
        # Split out sample name, gene, and template
        ssplit = sname.split('.')
        sample = ssplit[0]
        template = ssplit[-1]
        [gene, template] = template.split('-',2)
        
        # Split sample name into month and fraction
        ssplit = sample.split('-')
        if 'Blank' in ssplit:
            month = 'Blank'
            fraction = 'Blank'
        else:
            month = ssplit[0]
            fraction = ssplit[-2]    
        replicate = ssplit[-1]
        
        if month + '-' + fraction in list(correct_sample_info):
            fraction = correct_sample_info[month + '-' + fraction]
        
        headerlist.append((template, month, fraction, replicate))
        
    distMatrix.columns = pd.MultiIndex.from_tuples(
        headerlist, names = ['Template','Month','Fraction','Replicate'])
    distMatrix.index = distMatrix.columns
    
    return distMatrix


def buildPlot(data, color_category, figsize, filename):
    '''Produces clustered heatmap plot'''
    
    '''
    Old code: shows either month or fraction along top axis
    
    if color_category == 'month':
        # Map months to colors
        months = data.columns.get_level_values('Month')
        color_set = pd.Series(months, index = data.columns).map(month_cmap)
    elif color_category == 'fraction':
    # Map fractions to colors
        fractions = data.columns.get_level_values('Fraction')
        color_set = pd.Series(fractions, index = data.columns).map(
            fraction_cmap)
    '''
    
    # Map months to colors along the top of the plot
    months = data.columns.get_level_values('Month')
    col_colors = pd.Series(months, index = data.columns).map(month_cmap)
    
    # Map fractions to colors along the side of the plot
    fractions = data.index.get_level_values('Fraction')
    row_colors = pd.Series(fractions, index = data.index).map(fraction_cmap)
    
    # Convert distance matrix to condensed form that seaborn expects
    condensed = squareform(data.values)
    # Generate linkage for rows and columns
    row_linkage = linkage(condensed, method = 'single')
    col_linkage = linkage(condensed.T, method = 'single')
    
    # Heatmap
    p = sns.clustermap(
        data, figsize = (figsize,figsize),
        row_linkage = row_linkage,
        col_linkage = col_linkage,
        cmap = cmap, col_colors = col_colors, row_colors = row_colors,
        metric=None # Don't need to re-calculate distances since it's a distance matrix
        )
    
    p.ax_heatmap.set_xlabel('Sample')
    p.ax_heatmap.set_ylabel('Sample')
    p.cax.set_title('Distance')

    
    # Save heatmap
    p.savefig(filename + '.svg', transparent = True)
    
    return p



def genHTMLPage(headstr, filepath):   
    '''Builds HTML pages for all of the images'''    
    for gene in genes:
        for s in sets:
            htmltxt = (headstr
                       + '<h1>' + gene + ' clustered heatmaps '
                       + '(' + sets[s]['title'] + ')</h1>'
                       + '<p><img src = "MonthFractionKey_horiz.png"'
                       + ' height = 150></p>'
                       )
            
            if s == 'all':
                htmltxt = (htmltxt
                           + '<p><img src = "'
                           + file_pfx + '_' + gene + '_all.svg", alt = "'
                           + gene + ' all"></p>')
           
            elif s == 'template':
                for template in templates:
                    htmltxt = (htmltxt
                               + '<h2>' + gene + ' ' + template
                               + '</h2><p><img src = "' + file_pfx + '_'
                               + gene + '_' + template + '.svg", alt = "'
                               + gene + ' ' + template + '"></p>')
            
            elif s == 'month':
                htmltxt = htmltxt + '<table>'
                for month in months:
                    htmltxt = htmltxt + '<tr>'
                    for template in templates:
                        htmltxt = (htmltxt
                                   + '<td><h2>' + gene + ' '
                                   + months[month] + ' ' + template
                                   + '</h2><img src = "' + file_pfx + '_'
                                   + gene + '_' + month + '_' + template
                                   + '.svg", alt = "' + gene + ' ' + month
                                   + ' ' + template + '"></td>')
                    htmltxt = htmltxt + '</tr>'
                htmltxt = htmltxt + '</table>'
            
            elif s == 'fraction':
                htmltxt = htmltxt + '<table>'
                for fraction in fractions:
                    htmltxt = htmltxt + '<tr>'
                    for template in templates:
                        htmltxt = (htmltxt
                                   + '<td><h2>' + gene + ' ' 
                                   + fractions[fraction]['title'] + ' '
                                   + template + '</h2><img src = "'
                                   + file_pfx + '_' + gene + '_' + fraction
                                   + '_' + template + '.svg", alt = "'
                                   + gene + ' ' + fractions[fraction]['title']
                                   + ' ' + template + ' ' + '"></td>')
                    htmltxt = htmltxt + '</tr>'
                htmltxt = htmltxt + '</table>'
            
            hfile = open(filepath + '_' + gene + '_' + s + '.html', 'w')
            hfile.write(htmltxt)
            hfile.close()
    
#%%
# Main function


# Load the distance matrices
directory = os.getcwd()
for gene in genes:
    filename, directory, distance_matrix = RottenIceModules.fileGet(
        'Select ' + gene + ' distance matrix',  file_type = 'tsv',
        header_row = 0, index_col = 0, directory = directory)
    genes[gene]['distMatrix'] = parseSamples(distance_matrix)

# Determine the filepath
if outdir != None: # If directory is specified in the initial variables
    filepath_pfx = outdir + '\\' + file_pfx
else: # Otherwise use the current working directory
    filepath_pfx = os.getcwd() + '\\' + file_pfx

# Run the analysis and build the plots
for gene in genes:
    # Get gene's info
    filename = filepath_pfx + '_' + gene
    distMatrix = genes[gene]['distMatrix']
    
    # Plot set 1: all data
    buildPlot(distMatrix, 'month', 50, filename + '_all')
    
    # Plot set 2: cDNA and DNA in seperate figures
    # Split datasets into cDNA and DNA
    datasets = {
        'DNA'   : {
            'all'   : distMatrix.loc['DNA']['DNA']
            },
        'cDNA'  : {
            'all'   : distMatrix.loc['cDNA']['cDNA']
            }
        }
    for ds in datasets:
        # Build and save plot of all data for the template
        print('Building ' + gene + ' ' + ds + ' plots...')
        buildPlot(datasets[ds]['all'], 'month', 20, filename + '_' + ds)
        
        # Split out different months
        for month in months:
            print(month)
            datasets[ds][month] = datasets[ds]['all'].loc[
                datasets[ds]['all'].index.get_level_values('Month') == month,
                datasets[ds]['all'].columns.get_level_values('Month') == month
            ]
            buildPlot(datasets[ds][month], 'fraction', 10,
                      filename + '_' + month + '_' + ds)
        
        # Split out different fractions
        for ftype in fractions:
            print(ftype)
            f = fractions[ftype]['subset']
            # Make sure the fractions are represented
            mask_rows = datasets[ds]['all'].index.get_level_values('Fraction').isin(f)
            mask_cols = datasets[ds]['all'].columns.get_level_values('Fraction').isin(f)
            # Split out the dataset to plot
            datasets[ds][ftype] = datasets[ds]['all'].loc[mask_rows, mask_cols]
            buildPlot(datasets[ds][ftype], 'month', 10,
                      filename + '_' + ftype + '_' + ds)
        

# Generate HTML pages
# Build the file navigation HTML
nav_html = RottenIceVars.nav_html_start
for gene in genes:
    set_html = '  <em>' + gene + ': </em>'
    for s in sets:
        set_html = (set_html + ' <a href = "' + file_pfx + '_' + gene + '_' 
                    + s + '.html">' + sets[s]['title'] + '</a> /')
    nav_html = nav_html + set_html[:-1] + '  //  '
nav_html = nav_html[:-6] 
# Add the header information
headstr = RottenIceModules.genHTMLhead(
    title, page_nav_html = nav_html, subtitle_text = subtitle_text)
genHTMLPage(headstr, filepath_pfx)
        

    
    
    

