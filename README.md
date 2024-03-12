# MardenPark

## Code and attached data for plot creation and analysis of allele frequency shifts
C.L.Metheringham and R.A.Nichols

### create_figure1.R
Uses the data in Supplementary Data 1 to create the stacked bar plot in Figure 1, showing the proportion of trees in each health catagory in 2019 and 2021.

### bglr.R
Test if the sites identified as having large effects in the field trial (10) contributed to the visually assessed dieback damage at Marden Park, a new genomic prediction was carried out using the BGLR (V1.1)(17) package’s implementation of the BayesB algorithm (50,000 iterations with a burn in of 2000). 

### pca.R
Create Extended Data Figure 2

## GEBV_pheno_corr.R
Create Extended Data Figure 3

###  gebv_extremes_plots.Rmd
Create Extended Data Figure 4

###  adult_juv_gebv_shift.R
Are juvenile GEBVs significantly different from adults
Accounting for ancestry of juvenile trees by unlinked sites

### allele_shift_correlation.R

### sequoia_parentage.R 
Used to obtain estimates of most likely parentage using sequoia 

###  parent_offspring_allele_shift.R
Uses effect_sizes.csv, gebv_site_frequencies.csv and related_trees.csv to test if allele frequencies differ from that expected from likely frequency

###  parent_offspring_correlation.R
Tests and plots correlation of adult GEBV and score with offspring health

## Data
* effect_sizes.csv
* gebv_site_frequencies.csv
* related_trees.csv

* See Suppliment 2 for simulated analysis
