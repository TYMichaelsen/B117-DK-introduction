# B117-DK-introduction
All code used to run the phylogenetic analysis, Poisson regression models, and generate the visualizations used in Michaelsen et al. 2021

**NOTE** The code here is not self-contained, as the raw genomic and epidemiological data used in this work is person-sensitive and we are prohibited to make it publicly available according to Danish legislation. SARS-CoV-2 consensus genome sequences associated with this work have been uploaded to the GISAID database in accordance with Danish law (no. 285 vers. 2021-02-27). As consequence, dates are binned by week and maximum spatial resolution is at regional level. 

**phylogenetic-analysis.sh** contains all code used to perform the phylogenetic analysis and ancestral state reconstruction. 

**analysis.Rmd** contains the R-code used to analyse all epidemiological data and the output from the ancestral state reconstruction. Contains all code used to create visualizations for the manuscript.
