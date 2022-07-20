#Correlation of effect size with alle frequency shift
#Based on methods and code by Richard A. Nichols

#Import effect sizes
es_file <- read.csv("MP_effects_MIA_and_MAA.csv")
es <- with(es_file, EES.MIA - EES.MAA)

#Import csv files with data on the GEBV sites and unlinked sites
gebv_sites <- read.csv(file = "gebv_sites.csv")
#Extract the genotype matrices
GTgebv <- gebv_sites[, 4:ncol(gebv_sites)]
GTAgebv <- gebv_sites[gebv_sites$Age=="Adult", 4:ncol(gebv_sites)]
GTJgebv <- gebv_sites[gebv_sites$Age=="Juv", 4:ncol(gebv_sites)]
#Get frequencies
fAll <- colMeans(GTgebv)/2
fAdult <- colMeans(GTAgebv)/2
fJuvenile <- colMeans(GTJgebv)/2
#Shift to use minor allele frequency

for (snp in 1:length(fAll)){
  if(fAll[snp] > 0.5){
    fAdult[snp] <- 1 - fAdult[snp]
    fJuvenile[snp] <- 1 - fJuvenile[snp]
    fAll[snp] <- 1 - fAll[snp]
    es[snp] <- -es[snp]
  }
}

#Change in frequency between adults and juveniles
pdiff <- fJuvenile - fAdult

# calculate mean allele frequency for each locus and hence p^2 pq and q^2
# and fitnesses for each genotype (assuming additive effect ~es)
pmean <- fAll
pq <- pmean*(1-pmean)
p2 <- pmean^2
q2 <- (1-pmean)^2
pqes <- pq*es
whet <- 1 + es
whom <- 1 + 2*es

# expected difference in allele frequency (should be proportional to pq)
dexp <- pmean-(2*p2 +2*pq*whet)/(2*p2 + 4*pq*whet + 2*q2*whom)

mod0 <- lm(pdiff ~ 1)
mod1 <- lm(pdiff ~ es)
mod2 <- lm(pdiff ~ dexp)
mod3 <- lm(pdiff ~ es + dexp)

anova(mod0, mod1)
anova(mod0, mod2)
anova(mod2, mod3)
summary(mod2)

# Plot relationship of pdiff with es*pq
# Reduce noise by getting average change in allele frequency 
# in 200 equal quantiles of effect size
bins <-200
# find that number of es*pq quantiles by finding bin boundaries
brks <- quantile(es*pq, seq(0,1, length.out = bins+1))
# set up variates to hold data about the bins
dmean <- espqmean <- ese <- rep(0,bins)
for (i in 1:bins) {
  # Identify loci in quantile i
  choice <- (es*pq)>brks[i] & (es*pq)<brks[i+1]
  # Calculate the mean of the allele frequency difference
  dmean[i] <- mean(pdiff[choice])
  # Calculate the mean effect size
  espqmean[i] <- mean((es*pq)[choice])
  # Calculate the ES on mean change in allele frequency
  ese[i] <- sd(pdiff[choice])/sqrt(sum(choice))
}

#Format binned data for plotting
binned_table <- data.frame(ESPQ = espqmean, DIFF = dmean, ESE = ese)

#Get slope and intercept of model2
int <- coefficients(mod2)[1]
slope <- coefficients(mod2)[2]

library(ggplot2)
fig5 <- ggplot(data = binned_table, aes(x = ESPQ, y = DIFF))+
  geom_point(aes(alpha = 0.2, size = 0.7/ese))+
  geom_abline(intercept = int, slope = slope, colour = "red")+
  theme_minimal()+
  theme(legend.position = "none")+
  ylab("Mean change in allele frequency")+
  xlab("Effect size x pq (200 quantiles)")
  
tiff("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_fig5.tiff", units="mm", width=180, height=100, res=300)
fig5
dev.off()
png("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_fig5.png", units="mm", width=180, height=100, res=300)
fig5
dev.off()
jpeg("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_fig5.jpeg", units="mm", width=180, height=100, res=300)
fig5
dev.off()