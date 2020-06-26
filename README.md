# rottenice
Data processing python scripts used for the Rotten Ice project
by Carie Frantz cariefrantz@weber.edu

## List of scripts

<table>
<tr><th>Script</th><th>Description</th><th>Data files used</th></tr>
<tr><td>BetaClustermaps.py</td><td>Used to generate sample clustermaps from sample distance matrices</td><td>16S_dist_matrix.csv, 18S_dist_matrix.csv</td></tr>
<tr><td>ChemDataPCA.py</td><td>Used to generate PCA plots of physical and chemical habitat parameters</td><td>metadata.csv</td></tr>
<tr><td>ChemHeatmaps.py</td><td>Used to generate heatmaps of sample metadata (physical, chemical, and biological measurements)</td><td>metadata.csv</td></tr>
<tr><td>SpearmanGrid.py</td><td>Used to generate a heatmap grid of Spearman correlation coefficients for metadata vs. taxonomic group abundance</td><td>metadata.csv, 16S_cDNA_ASV-table.csv, 18S_cDNA_ASV-table.csv</td></tr>
<tr><td>WebTaxBarplotsSimple.py</td><td>Used to generate interactive side-by-side bokeh taxonomic bar plots (web-based) from algae taxonomy data. Results at different taxonomic levels are built in as tabs.</td><td>AlgaeIDs.csv</td></tr>
<tr><td>WebTaxBarplotsStacked.py</td><td>Used to generate interactive stacked bokeh taxonomic bar plots (web-based) from Illumina sequencing data. Plots for absolute and relative abundance for both DNA and cDNA are stacked for ease of comparison.</td><td>16S_DNA_ASV-table.csv, 16S_cDNA_ASV-table.csv, 18S_DNA_ASV-table.csv, 18S_cDNA_ASV-table.csv</td></tr>
</table>

## Setting up
The code for this project requires the following list of packages in order to run.
<ul>
<li>tkinter</li>
<li>progress</li>
<li>numpy</li>
<li>pandas</li>
<li>math</li>
<li>scipy</li>
<li>sklearn</li>
<li>seaborn</li>
<li>matplotlib</li>
<li>bokeh</li>
</ul>

To install using conda, execute the command:

	conda install tkinter
	conda install progress
	
...and so on

To install using pip, execute the command:

	pip install tkinter
	pip install progress
	
...and so on

## To do

<ul>
<li>WebTaxBarplotsSimple: Add grouping module</li>
<li>ChemDataPCA: Determine which metadata parameters are available (no nan) based on sample set</li>
<li>BetaClustermap: Fix key colors</li>
</ul>
