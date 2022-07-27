#Are alleles associated with ADB found at different than expected 
#frequencies in offsping trees compared to likely parents

#Import csv files with data on the GEBV sites and unlinked sites
gt <- read.csv(file = "gebv_site_frequencies.csv")
#Read in the related trees
related <- read.csv(file = "related_trees.csv")

#Pull out the frequencies of related trees
#Format allele frequencies as tables
juv_df <- matrix(, ncol = nrow(gt), nrow = nrow(related))
dam_df <- matrix(, ncol = nrow(gt), nrow = nrow(related))
sire_df <- matrix(, ncol = nrow(gt), nrow = nrow(related))
for (i in 1:nrow(related)) {
  juv_df[i,] <- gt[, which(colnames(gt) == related$id[i])]
  dam_df[i,] <- gt[, which(colnames(gt) == related$dam[i])]
  if (is.na(related$sire[i])) {
    sire_df[i,] <- rowMeans(gt)
  }
  else{
    sire_df[i,] <- gt[, which(colnames(gt) == related$sire[i])]
  }
}

#Get the expected allele frequency for the juveniles
exp_freq <- matrix(, nrow = nrow(gt), ncol = nrow(related))
freq_diff <- matrix(, nrow = nrow(gt), ncol = nrow(related))
mean_exp_freq <- c()
for (i in 1:ncol(juv_df)) {
  exp_freq[i,] <- (dam_df[, i] + sire_df[, i]) / 2
  mean_exp_freq[i] <- mean(exp_freq[i,], na.rm = T)
  #Difference between expected and observed at this site
  freq_diff[i,] <-  juv_df[,1] - exp_freq[i,]
}

#Test if there is a significant difference between expected and observed frequencies
diff_test <- c()
for (i in 1:ncol(juv_df)) {
  diff_test[i] <- wilcox.test(exp_freq[[i]], juv_df[, i])$p.value
}

#adjust p values and reformat as table
#Get the mean difference between observed and expected at each site
diff_df <-
  data.frame(
    SNP = row.names(gt),
    p = diff_test,
    q = p.adjust(diff_test),
    mean_exp = mean_exp_freq,
    mean_obs = colMeans(juv_df, na.rm = T),
    mean_diff = rowMeans(freq_diff)
  )

#Import effect sizes for each site - Supplementary Data 2
es_file <- read.csv("MP_effects_MIA_and_MAA.csv")
es <- data.frame(SNP = es_file$SNP, ES = with(es_file, EES.MIA - EES.MAA))
#Test for correlation between mean_diff and effect size
ees_diff <- merge(diff_df, es, by = "SNP")
cor.test(ees_diff$ES, ees_diff$mean_diff)
#Sites with negative effects more often seen at lower than expected
# frequencies in juvenile trees compared to their assigned parents

#How many sites have significant difference between expected
#and observed allele frequencies
sum(p.adjust(diff_test) < 0.05)