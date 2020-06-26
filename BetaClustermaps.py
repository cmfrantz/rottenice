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
import pandas as pd
import seaborn as sns

import matplotlib
# Define the font type to make exported plots editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import RottenIceModules
import RottenIceVars


####################
# VARIABLES
####################

# Plot / HTML file title info
title = RottenIceVars.file_sets['beta_cluster']['title']
subtitle_text = ('Heatmaps generated from weighted UniFrac distance matrices '
                 + 'calculated on quality-controlled Rotten Ice sample '
                 + 'amplicon sequences. Samples were grouped using '
                 + 'nearest-point clustering on correlation distances using '
                 + 'the <a href="https://seaborn.pydata.org/generated/'
                 + 'seaborn.clustermap.html">seaborn.clustermap</a> package '
                 + '(v. 0.10.1) for Python. Data processed using the script '
                 + '<a href="https://github.com/cmfrantz/rottenice">'
                 + 'BetaClustermaps.py</a>. C. Frantz, May 2020</p>')
file_pfx = RottenIceVars.file_sets['beta_cluster']['pfx']
alt_text_pfx = 'Sample clustermap: '

# Samples to include
months = ['M', 'JN', 'JY10', 'JY11']

fractions = {
    'I' : {
        'title'     : 'Ice-only melts',
        'subset'    : ['IT','IM','IB']
        },
    'H' : {
        'title'     : 'Whole-horizon melts',
        'subset'    : ['HT','HM','HB']
        },
    'B' : {
        'title'     : 'Brines',
        'subset'    : ['BT', 'BM', 'BB', 'B']
        },
    'O': {
        'title'     : 'Other fluids',
        'subset'    : ['P1','P2','PW','SW','Drain']
        }
    }

# Plot variables
cmap = RottenIceVars.cmap
genes = {
    '16S'   : {},
    '18S'   : {}
    }
templates = ['DNA', 'cDNA']
sets = {
    'all'       : {'title' : 'all'},
    'template'  : {'title' : 'by template'},
    'month'     : {'title' : 'by month'},
    'fraction'  : {'title' : 'by fraction'}
    }

# Categorical color palette to identify the months
# Based on this palette: https://www.color-hex.com/color-palette/22257
month_cmap = {
    'M'     : RottenIceVars.plotColorsByMonth['M'],
    'JN'    : RottenIceVars.plotColorsByMonth['JN'],
    'JY10'  : RottenIceVars.plotColorsByMonth['JY10'],
    'JY11'  : RottenIceVars.plotColorsByMonth['JY11'],
    'Blank' : '#FFFFFF'     # White
    }
    
# Categorical color palette to identify the fractions
# Based on this palette: https://www.color-hex.com/color-palette/91134
fraction_cmap = {
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
    'B'     : [75,139,167],
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
    'B' : {
        'title'     : 'Brines',
        'subset'    : ['BT', 'BM', 'BB', 'B']
        },
    'O': {
        'title'     : 'Other fluids',
        'subset'    : ['P1','P2','PW','SW','Drain']
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
        name = ssplit[-1]
        headerlist.append((template, month, fraction, name))
        
    distMatrix.columns = pd.MultiIndex.from_tuples(
        headerlist, names = ['Template','Month','Fraction','Sample'])
    distMatrix.index = distMatrix.columns
    
    return distMatrix


def buildPlot(data, color_category, figsize, filename):
    '''Produces clustered heatmap plot'''
        
    if color_category == 'month':
        # Map months to colors
        months = data.columns.get_level_values('Month')
        color_set = pd.Series(months, index = data.columns).map(month_cmap)
    elif color_category == 'fraction':
    # Map fractions to colors
        fractions = data.columns.get_level_values('Fraction')
        color_set = pd.Series(fractions, index = data.columns).map(
            fraction_cmap)
    
    # Heatmap
    p = sns.clustermap(data, metric="correlation", method="single",
                       cmap=cmap, figsize=(figsize,figsize),
                       col_colors = color_set)
    
    # Save heatmap
    p.savefig(filename + '.svg', transparent = True)
    
    return p


def genHTMLPage(headstr):   
    '''Builds HTML pages for all of the images'''    
    for gene in genes:
        for s in sets:
            htmltxt = (headstr
                       + '<h1>' + gene + ' clustered heatmaps '
                       + '(' + sets[s]['title'] + ')</h1>')
            
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
            
            hfile = open(file_pfx + '_' + gene + '_' + s + '.html', 'w')
            hfile.write(htmltxt)
            hfile.close()
    
#%%
# Main function

# Build the file navigation HTML
nav_html = RottenIceVars.nav_html_start
for gene in genes:
    set_html = '  <em>' + gene + ': </em>'
    for s in sets:
        set_html = (set_html + ' <a href = "' + file_pfx + '_' + gene + '_' 
                    + s + '.html">' + sets[s]['title'] + '</a> /')
    nav_html = nav_html + set_html[:-1] + '  //  '
nav_html = nav_html[:-6] 

# Load the distance matrices
for gene in genes:
    filename, directory, distance_matrix = RottenIceModules.fileGet(
        'Select ' + gene + ' distance matrix',  file_type = 'csv',
        header_row = 0, index_col = 0)
    genes[gene]['distMatrix'] = parseSamples(distance_matrix)

filepath_pfx = directory + '\\' + file_pfx

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
            datasets[ds][month] = datasets[ds]['all'].loc[month][month]
            buildPlot(datasets[ds][month], 'fraction', 10,
                      filename + '_' + month + '_' + ds)
        
        # Split out different fractions
        for ftype in fractions:
            print(ftype)
            f = fractions[ftype]['subset']
            fdf = datasets[ds]['all'].loc[:, pd.IndexSlice[:,f]]
            datasets[ds][ftype] = fdf.loc[pd.IndexSlice[:,f],:]
            buildPlot(datasets[ds][ftype], 'month', 10,
                      filename + '_' + ftype + '_' + ds)
        

# Generate HTML pages
headstr = RottenIceModules.genHTMLhead(
    title, page_nav_html = nav_html, subtitle_text = subtitle_text)
genHTMLPage(headstr)
        

    
    
    

