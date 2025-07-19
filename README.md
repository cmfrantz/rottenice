# rottenice
Data processing python scripts used for the Rotten Ice project
by Carie Frantz cariefrantz@weber.edu

## List of scripts
<table>
<tr><th>Script</th><th>Description</th><th>Data files used</th></tr>
<tr><td>AlphaDiversityPlots.py</td><td>Used to generate alpha diversity metric boxplots from rarefied sequencing sample metrics</td><td>alpha_diversity_metrics_16S.csv, alpha_diversity_metrics_18S.csv</td></tr>
<tr><td>BetaClustermaps.py</td><td>Used to generate sample clustermaps from sample distance matrices</td><td>16S_dist_matrix.csv, 18S_dist_matrix.csv</td></tr>
<tr><td>ChemDataPCA.py</td><td>Used to generate PCA plots of physical and chemical habitat parameters</td><td>metadata.csv</td></tr>
<tr><td>ChemHeatmaps.py</td><td>Used to generate heatmaps of sample metadata (physical, chemical, and biological measurements)</td><td>metadata.csv</td></tr>
<tr><td>SpearmanGrid.py</td><td>Used to generate a heatmap grid of Spearman correlation coefficients for metadata vs. taxonomic group abundance</td><td>metadata.csv, 16S_cDNA_ASV-table.csv, 18S_cDNA_ASV-table.csv</td></tr>
<tr><td>WebTaxBarplotsSimple.py</td><td>Used to generate interactive side-by-side bokeh taxonomic bar plots (web-based) from algae taxonomy data. Results at different taxonomic levels are built in as tabs.</td><td>AlgaeIDs.csv</td></tr>
<tr><td>WebTaxBarplotsStacked.py</td><td>Used to generate interactive stacked bokeh taxonomic bar plots (web-based) from Illumina sequencing data. Plots for absolute and relative abundance for both DNA and cDNA are stacked for ease of comparison.</td><td>16S_DNA_ASV-table.csv, 16S_cDNA_ASV-table.csv, 18S_DNA_ASV-table.csv, 18S_cDNA_ASV-table.csv</td></tr>
<tr><td>CommunityPCoA.py</td><td>Generates PCoA plots split out and coded by metadata from pcoa ordinate data.</td><td>ordination.txt for different datasets of interest</td></tr>
<tr><td>displot_TMB.py</td><td>Generates box-and-whisker plots comparing community distances in different ice horizons during different months to assess vertical homogenization.</td><td>distance-matrix.tsv</td></tr>
<tr><td>MetadataPCA.py</td><td>Similar to ChemDataPCA, but generates PCA biplots comparing samples based on their metadata. Displays vectors showing how different metadata parameters influence the principal components.</td><td>metadata.tsv</td></tr>
<tr><td>formatASVtable.py</td><td>Compiles ASV and taxonomy tables exported from QIIME2 as an excel file with tabs for counts and normalized data that is sorted by sample and taxonomy.</td><td>asv-table.tsv (taxonomy-added asv table from BIOM export) and otu-table-L7.csv (from QIIME2 View barplot export)</td></tr>
<tr><td>TaxBarplotPretty.py</td><td>Creates stacked barplots for publication from annotated relative abundance tables.</td><td>ASV-table-formatted-annotated.csv for datasets of interest</td></tr>
</table>

## Shared scripts
The following scripts are called by several (most) of the main scripts listed above. Download these to the same directory as the scripts above.
<table>
<tr><th>Script</th><th>Description</th></tr>
<tr><td>RottenIceModules.py</td><td>Collection of scripts shared by multiple other scripts in this project</td></tr>
<tr><td>RottenIceVars.py</td><td>Collection of variables shared by multiple other scripts in this project</td></tr>
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

## Running
Once python and the packages listed above are installed, to run a script from command line execute the command:

	python ChemDataPCA.py
	python WebTaxBarplotsSimple.py
	
...and so on
