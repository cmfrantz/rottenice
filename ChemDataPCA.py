<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:50:09 2019

@author: cariefrantz

This script creates 2D PCA plots from imported metadata
It was created as part of the Rotten Ice Project

Copyright (C) 2019  Carie M. Frantz

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

import pandas as pd
import numpy as np
import os
from tkinter import filedialog
from tkinter import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


# %%

# Read the bio-chem data in from master Excel file
root = Tk()
root.filename = filedialog.askopenfilename(initialdir = "/", title = "Select bio-chem raw data master file (*.xls)")
filename = root.filename
print (filename)                        # Prints the filename
root.destroy()                          # Closes the user input window
dirPath = os.path.dirname(os.path.abspath(filename)) # Get the file directory from the filename
dataSet = pd.read_excel(filename,sheet_name="meta_for_python")  # Reads in all of the data from the "meta_for_python" sheet
dataSet.columns = dataSet.iloc[1]   # Sets the column names to the short form
dataSet = dataSet.reindex(dataSet.index.drop([0,1])) # Drops the extra columns / header rows
# dataSet = dataSet.set_index('SampleID') # Use sample names as indices
# %%
# Define and normalize the data
# (PCA is affected by scale)
# Update this if spreadsheet is updated with all numerical, scalable values
#List of parameters that could be used: "temperature", "Salnity", "bulk_density", "DOC", "Chl", "Phaeo", "FoFa", "PAM", "bact_cell_ct", "CTC", "bact_active_cell_ct", "pEPS", "POC", "Nitrogen", "CN"
features = ["DOC","POC","pEPS","Nitrogen"]
plotTags = ['month','SampleID','plot_color','plot_shape']
x = dataSet.loc[:,features].values  # Grabs the numberical data
y = pd.DataFrame(dataSet.loc[:,plotTags].values) # Defines the target (the dependent variable)
y.columns = plotTags    # Add column headers
x = StandardScaler().fit_transform(x)   # Normalize the data

# %%
# PCA
pca = PCA(n_components=2)   # 2-D PCA
principalComponents = pca.fit_transform(x)  # Perform PCA on the data
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1','principal component 2']) # Add headers to the result
finalDf = pd.concat([principalDf, y], axis=1)   # Add sample info & plot settings to the result
finalDf = finalDf.set_index('SampleID') # Use sample names as indices
explained_var = pca.explained_variance_ratio_
# %%
# Plot

# Create the plot & label axes
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1 (' + str(round(explained_var[0],2)*100) + '%)', fontsize = 15)
ax.set_ylabel('Principal Component 2 (' + str(round(explained_var[1],2)*100) + '%)', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)

# Plot the data
samples = y['SampleID'].values  # List of samples to plot

for sample in samples:  # Plot each sample individually
    ax.scatter(finalDf.loc[sample, 'principal component 1']     # x value
    , finalDf.loc[sample, 'principal component 2']  # y value
    , c=finalDf.loc[sample,'plot_color']    # marker color
    , s=50  # marker size
    , marker=finalDf.loc[sample,'plot_shape'])  # marker shape

ax.legend(samples)  # add legend
# ax.grid() # add grid (optional)
fig.savefig(dirPath + "\\metadata_PCA_" + '-'.join(features) + ".pdf", transparent=True) # Save figure

=======
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:50:09 2019

@author: cariefrantz

This script creates 2D PCA plots from imported metadata
It was created as part of the Rotten Ice Project

Copyright (C) 2019  Carie M. Frantz

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

import pandas as pd
import numpy as np
import os
from tkinter import filedialog
from tkinter import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


# %%

# Read the bio-chem data in from master Excel file
root = Tk()
root.filename = filedialog.askopenfilename(initialdir = "/", title = "Select bio-chem raw data master file (*.xls)")
filename = root.filename
print (filename)                        # Prints the filename
root.destroy()                          # Closes the user input window
dirPath = os.path.dirname(os.path.abspath(filename)) # Get the file directory from the filename
dataSet = pd.read_excel(filename,sheet_name="meta_for_python")  # Reads in all of the data from the "meta_for_python" sheet
dataSet.columns = dataSet.iloc[1]   # Sets the column names to the short form
dataSet = dataSet.reindex(dataSet.index.drop([0,1])) # Drops the extra columns / header rows
# dataSet = dataSet.set_index('SampleID') # Use sample names as indices
# %%
# Define and normalize the data
# (PCA is affected by scale)
# Update this if spreadsheet is updated with all numerical, scalable values
#List of parameters that could be used: "temperature", "Salnity", "bulk_density", "DOC", "Chl", "Phaeo", "FoFa", "PAM", "bact_cell_ct", "CTC", "bact_active_cell_ct", "pEPS", "POC", "Nitrogen", "CN"
features = ["DOC","POC","pEPS","Nitrogen"]
plotTags = ['month','SampleID','plot_color','plot_shape']
x = dataSet.loc[:,features].values  # Grabs the numberical data
y = pd.DataFrame(dataSet.loc[:,plotTags].values) # Defines the target (the dependent variable)
y.columns = plotTags    # Add column headers
x = StandardScaler().fit_transform(x)   # Normalize the data

# %%
# PCA
pca = PCA(n_components=2)   # 2-D PCA
principalComponents = pca.fit_transform(x)  # Perform PCA on the data
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1','principal component 2']) # Add headers to the result
finalDf = pd.concat([principalDf, y], axis=1)   # Add sample info & plot settings to the result
finalDf = finalDf.set_index('SampleID') # Use sample names as indices
explained_var = pca.explained_variance_ratio_
# %%
# Plot

# Create the plot & label axes
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1 (' + str(round(explained_var[0],2)*100) + '%)', fontsize = 15)
ax.set_ylabel('Principal Component 2 (' + str(round(explained_var[1],2)*100) + '%)', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)

# Plot the data
samples = y['SampleID'].values  # List of samples to plot

for sample in samples:  # Plot each sample individually
    ax.scatter(finalDf.loc[sample, 'principal component 1']     # x value
    , finalDf.loc[sample, 'principal component 2']  # y value
    , c=finalDf.loc[sample,'plot_color']    # marker color
    , s=50  # marker size
    , marker=finalDf.loc[sample,'plot_shape'])  # marker shape

ax.legend(samples)  # add legend
# ax.grid() # add grid (optional)
fig.savefig(dirPath + "\\metadata_PCA_" + '-'.join(features) + ".pdf", transparent=True) # Save figure

>>>>>>> 342ad789dfd262dcea074b050d0df35a278136d9
