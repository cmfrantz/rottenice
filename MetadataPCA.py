# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 10:24:04 2025

@author: cariefrantz
"""

####################
# IMPORTS
####################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import RottenIceModules

####################
# VARIABLES
####################

# Define the samples to analyze
samples = [
    'M-CS-HT','M-CS-HM','M-CS-HB',
    'JN-CS-HT','JN-CS-HM','JN-CS-HB',
    'JY10-HT','JY10-HM','JY10-HB',
    'JY11-HT','JY11-HM','JY11-HB'
    ]

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
# FUNCTIONS
####################

def biplot(scores, loadings, meta_df, color_col, color_map,
           marker_col, marker_map, variable_names, title,
           n_arrows=None, marker_size = 10):
    '''
    Creates a PCA biplot of metadata

    Parameters
    ----------
    scores : Array of float with dimension [n samples, m pcoa axes]
        PCA-transformed sample data (from pca.transform()).
    loadings : Array of float with dimension [n metadata parameters, m pcoa axes]
        PCA component loadings (from pca.components).T).
    meta_df : pandas.DataFrame
        Full metadata DataFrame (samples as rows, metadata as columns) including
        the columns used to map color and marker types (can be non-numeric).
        Note that this does NOT need to be the same dataframe used for PCA
        (which can only use numeric data) so long as the sample and metadata
        (column) names used in PCA are consistent.
    color_col : str
        Metadata column that determines the marker colors (used in color_map).
    color_map : dict of str
        Dictionary mapping values in color_col to colors.
    marker_col : str
        Metadata column that determines the marker shape (used in marker_map).
    marker_map : dict of str
        Dictionary mapping values in marker_col to marker shapes.
    variable_names : list of str
        List of the variable names (used for labeling the variable vector
        arrows) in the same order as in the PCA input. Should have the same
        dimension as the length of loadings.
    title : str
        Title for labeling the plot
    n_arrows : int, optional
        The maximum number of variable vector arrows to display.
        If defined, the plot will display only the N most important variables.
        The default is None.
    arrow_scale : float, optional
        Factor used to scale the variable vector arrows. The default is 1.0.
        Increasing this makes the arrows longer/larger.
    marker_size : int, optional
        Size of the sample markers in the PCA plot. The default is 10.

    Returns
    -------
    fig : matplotlib figure
        handle for the figure

    '''

    # Set up plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Extract groupings
    color_vals = meta_df[color_col]
    marker_vals = meta_df[marker_col]

    # Plot samples
    for c_val in color_map:
        for m_val in marker_map:
            idx = (color_vals == c_val) & (marker_vals == m_val)
            if idx.any():
                ax.scatter(scores[idx, 0], scores[idx, 1],
                           color=color_map[c_val],
                           marker=marker_map[m_val], s = marker_size,
                           label=f"{c_val} / {m_val}",
                           edgecolor='black', alpha=0.8)
    
    # Compute the vector magnitude (importance)
    vector_lengths = np.linalg.norm(loadings, axis=1)
    
    # Select the top N variables if n_arrows is capped
    if n_arrows is not None and n_arrows < len(loadings):
        top_indices = np.argsort(vector_lengths)[-n_arrows:] # N largest
    else:
        top_indices = np.arange(loadings.shape[0]) # show all variables
        
    # Scale the arrows based on the plot size
    # Determine the range of PCA scores (sample scatter points)
    x_range = scores[:, 0].max() - scores[:, 0].min()
    y_range = scores[:, 1].max() - scores[:, 1].min()
    max_plot_radius = 0.3 * max(x_range, y_range)  # Half the width/height
    # Determine max length of loading vector
    arrow_lengths = np.linalg.norm(loadings, axis=1)
    max_arrow_length = arrow_lengths.max()
    # Compute scale factor to make longest arrow ~half the plot
    arrow_scale = max_plot_radius / max_arrow_length
    
    # Plot arrows for variable loadings    
    for i in top_indices:
        ax.arrow(0, 0,
                 loadings[i, 0] * arrow_scale,
                 loadings[i, 1] * arrow_scale,
                 color='gray', alpha=0.8,
                 head_width=0.03, head_length=0.05)
        ax.text(loadings[i, 0]  * arrow_scale * 1.15,
                loadings[i, 1] * arrow_scale * 1.15,
                variable_names[i], color='gray',
                ha='center', va='center', fontsize=9)

    # Label axes with explained variance
    pca_var = np.var(scores, axis=0) / np.sum(np.var(scores, axis=0))
    ax.set_xlabel(f"PC1 ({pca_var[0]*100:.1f}%)", fontsize=10, color = 'k')
    ax.set_ylabel(f"PC2 ({pca_var[1]*100:.1f}%)", fontsize=10, color = 'k')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, frameon=False)
    ax.set_title(title, fontsize=12)
    ax.set_facecolor('white')
    ax.grid(False)
    
    # Ensure all axis box lines are visible
    for side in ['top', 'right', 'bottom', 'left']:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(1.2)
        ax.spines[side].set_color('black')
        
    ax.tick_params(direction='out', length=6, width=1, colors='k')
    
    plt.tight_layout()
    plt.show()
    
    return fig




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
    
    # Run analysis and build plots for each set of variables
    for varset in varsets:
        
        # Build matrix with only the indicated metadata values
        meta_df = meta_trimmed[varsets[varset]]
        
        # Identify and remove non-numeric columns
        non_numeric_cols = meta_df.select_dtypes(exclude=['number']).columns
        for col in non_numeric_cols:
            print(f"Column '{col}' is non-numeric and will be removed.")
        meta_df_clean = meta_df.drop(columns=non_numeric_cols)
                
        # Delete any samples/rows with missing data
        meta_df_clean = meta_df_clean.dropna()
        clean_index = meta_df_clean.index
        
        # Subset and scale the data
        X = meta_df_clean.values
        X_scaled = StandardScaler().fit_transform(X)
        
        # Perform 2-dimensional PCA analysis
        pca = PCA(n_components=2)
        scores = pca.fit_transform(X_scaled)
        loadings = pca.components_.T  # shape (n_features, 2)
        expl_var = pca.explained_variance_ratio_
        
        # Align full metadata table to samples used for PCA
        meta_for_plot = meta_trimmed.loc[clean_index]
        
        # Define plot information
        var_names = list(meta_df_clean.columns)
        plot_title = 'Metadata PCA Biplot: ' + varset
        
        # Build the biplot
        biplot_fig = biplot(
            scores, loadings, meta_for_plot,
            color_col, color_map, marker_col, marker_map,
            var_names, plot_title,
            n_arrows = 10, marker_size = 100
            )
        
        # Save the biplot
        biplot_fig.savefig(
            directory + '\\' + 'Metadata_Biplot ' + varset + '.png', transparent = True)
        biplot_fig.savefig(
            directory + '\\' + 'Metadata_Biplot ' + varset + '.pdf', format = 'pdf')
        
        
