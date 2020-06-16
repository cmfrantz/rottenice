# rottenice
Data processing python scripts used for the Rotten Ice project
by Carie Frantz cariefrantz@weber.edu

List of scripts:

<table>
<tr><th>Script</th><th>Description</th><th>Data files used</th></tr>
<tr><td>ChemDataPCA.py</td><td>Used to generate PCA plots of physical and chemical habitat parameters</td><td>Master Bio-Chem Sample Sheet (Metadata).xls</td></tr>
<tr><td>ChemHeatmaps.py</td><td>Used to generate heatmaps of sample metadata (physical, chemical, and biological measurements)</td><td>Master Bio-Chem Sample Sheet (Metadata).xls</td></tr>
<tr><td>SpearmanGrid.py</td><td>Used to generate a heatmap grid of Spearman correlation coefficients for metadata vs. taxonomic group abundance</td><td>metadata.csv, 16S_cDNA_ASV-table.csv, 18S_cDNA_ASV-table.csv</td></tr>
</table>

To do:

<ul>
<li>Add scripts for...
	<ul>
	<li>Clustermap</li>
	<li>BarPlots</li>
	<li>AlgaeBarPlots</li>
	</ul>
</li>
<li>Update any scripts...
	<ul>
	<li>that use metadata xls to use final metadata csv</li>
	<li>that use ASV tables other than the final csv tables</li>
	</ul>
</li>
</ul>
