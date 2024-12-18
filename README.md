# Rapid polygenic adaptation in a wild population of ash trees under a novel fungal epidemic

## Code and attached data for plot creation and analysis of allele frequency shifts
C.L.Metheringham and R.A.Nichols

## Main text figure creation

### Figure1.R
Uses the data in Supplementary Data 1 to create the stacked bar plot in Figure 1, showing the proportion of trees in each health catagory in 2019 and 2021.

### Figure2.Rmd
Tests and plots correlation of adult GEBV and score with mean offspring health, weighted by number of offspring per adult tree

### Figure3.Rmd

Uses the GEBV of adult and juvenile trees as listed in Supplementary Data 1 and part of the code from GEBV-regression.Rmd (see below) to plot Figure 3

### allele_shift_correlation.R
Correlation of effect size with alle frequency shift.

Based on methods and code by Richard A. Nichols.

Creates figure 4: To visualize the trend more clearly we grouped the loci into 200 bins (quantiles of gAf(1-f) ). The means of each bin are plotted in with an area of each point being inversely proportional to the variance in f, to convey the relative precision of the mean festimate. The line is the fitted linear regression (carried out on the individual f values). 

# Additional analysis and creation of supplementary figures

### gebv_mp_bglr.R and bayesb_results.R
Test if the sites identified as having large effects in the field trial (10) contributed to the visually assessed dieback damage at Marden Park, a new genomic prediction was carried out using the BGLR (V1.1)(17) package’s implementation of the BayesB algorithm (50,000 iterations with a burn in of 2000). 

### pca_plots.Rmd
Create Extended Data Figure 2 from PLINK output

### GEBV_pheno_corr.R
Create Extended Data Figure 3 showing the correlation og GEBV with phenotype for both the adult and juvenile age cohorts

###  gebv_extremes_plots.Rmd
Create Extended Data Figure 4

### GEBV-regression.Rmd
Richard A. Nichols.

Analysis of allele frequency change between adults and juveniles

To assess whether the change in GEBV can be attributed to selection rather than genetic drift we adopt two approaches that exploit the genotypes at unlinked sites (sites not used in the GEBV calculations).

Firstly we use these loci to calculate a matrix of relatedness among the pairs of plants. This matrix is then used to predict the GEBV of juveniles, from that of related adults.  The regression between observed and predicted GEBV allows for genetic drift; i.e. the differential success the adults. The action of selection would be apparent in in intercept: a positive intercept would be expected with selection for higher GEBV. This effect would occur if the siblings of the surviving juveniles, which had lower GEBV had succumbed to the fungus.

### sequoia_pedigree_subset.Rmd 
Notebook used to estimate most likely parentage using sequoia.  

Parentage assignment was performed with the sequoia R package (V2.3.1) (19), using a randomly selected set of 1,000 SNPs having read depth > 20, minor allele frequency > 0.4, and an estimated error rate < 0.01. Parentage assignment was run using the hermaphrodite “B” mode. DBH was used to create a proxy for the birth year of adult trees (birth year = 100 - DBH) and juvenile trees were assigned a birth year of 100. This allowed for larger, and therefore presumably older, adult trees to be assigned as parents of younger adult trees, while disallowing juveniles from being considered as parent trees. The maximum age of parents was set as 99, allowing all adult trees to be considered as parents of the juveniles. The proxy years were not intended to be an accurate estimate of tree age. Confidence of parentage assignment was estimated within sequoia by simulating genotype data based on the estimated pedigree, recalculating parentage assignment based upon simulated data and comparing the recalculated pedigree to the original pedigree. The proportion of correct assignments over multiple runs of the process gives an estimate of confidence in our parentage assignments. 

###  parent_offspring_allele_shift.R
Uses effect_sizes.csv, gebv_site_frequencies.csv and related_trees.csv to test if allele frequencies differ from that expected from likely frequency

### GreenupRate.R
Richard A. Nichols.
Calculate the rate of greenup and plot Supplementary Figure 5

### Simulations.Rmd

## Data
* effect_sizes.csv
* gebv_site_frequencies.csv
* related_trees.csv
* unlinked_sites.csv

* See Supplement 2 for details on simulated analysis and Zenodo 10.5281/zenodo.10808942 for additional files
