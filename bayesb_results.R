#Process BayesB result

numSites <- 1000
reps <- 100

randomJ <- snp_effect_tableJ[1:numSites,]
large_effectJ <- snp_effect_tableJ[(numSites+1):(numSites*2),]

table((randomJ^2)>0.0000001)/(numSites*reps)
table((large_effectJ^2)>0.0000001)/(numSites*reps)

hist(abs(unlist(randomJ)), breaks = 100, xlab = "Absolute effect size", main = "Juveniles - random")
hist(abs(unlist(large_effectJ)), breaks = 100, xlab = "Absolute effect size", main = "Juveniles - large effect")
t.test(abs(unlist(randomJ)), abs(unlist(large_effectJ)))

res.aov.J <- aov(abs(unlist(randomJ)) ~ abs(unlist(large_effectJ)))
summary(res.aov.J)

randomA <- snp_effect_tableA[1:numSites,]
large_effectA <- snp_effect_tableA[(numSites+1):(numSites*2),]

table((randomA^2)>0.00001)/(numSites*reps)
table((large_effectA^2)>0.00001)/(numSites*reps)

hist(abs(unlist(randomA)), breaks = reps, xlab = "Absolute effect size", main = "Adults - random")
hist(abs(unlist(large_effectA)), breaks = reps, xlab = "Absolute effect size", main = "Adults - large effect")
t.test(abs(unlist(randomA)), abs(unlist(large_effectA)))

#Anova

res.aov <- aov(abs(unlist(randomA)) ~ abs(unlist(large_effectA)))
summary(res.aov)

#correlation with prev effect sizes
ees.A <- EES[which(large_effect_snps %in% rownames(GT3A)),]
cor.test(colMeans(large_effectA), -ees.A$EES.MIA)

var.test(abs(unlist(randomA)), abs(unlist(large_effectA)))
var.test(abs(unlist(randomJ)), abs(unlist(large_effectJ)))

var.test(unlist(large_effectJ), (unlist(randomJ)))
var.test(unlist(large_effectA), (unlist(randomA)))

var.test(rowMeans(large_effectJ), (unlist(randomJ)))
var.test(rowMeans(large_effectA), (unlist(randomA)))

FVals_Juv <- data.frame(rep(0,reps))
pVals_Juv <- data.frame(rep(0,reps))

for (i in 1:reps) {
  pVals_Juv[i] <-
    var.test((unlist(large_effectJ[[i]])), (unlist(randomJ[[i]])), alternative = "greater")$p.value
  FVals_Juv[i] <-
    var.test((unlist(large_effectJ[[i]])), (unlist(randomJ[[i]])), alternative = "greater")$statistic
}
hist(as.numeric(unlist(pVals_Juv)), breaks = 50)
hist(as.numeric(unlist(FVals_Juv)), breaks = 50)

FVals_Adult <- data.frame(rep(0,reps))
pVals_Adult <- data.frame(rep(0,reps))

for (i in 1:reps) {
  pVals_Adult[i] <-
    var.test((unlist(large_effectA[[i]])), (unlist(randomA[[i]])), alternative = "greater")$p.value
  FVals_Adult[i] <-
    var.test((unlist(large_effectA[[i]])), (unlist(randomA[[i]])), alternative = "greater")$statistic
}
hist(as.numeric(unlist(pVals_Adult)), breaks = 50)
hist(as.numeric(unlist(FVals_Adult)), breaks = 50)

f_pJ <-  fisher.method(pVals_Juv, p.corr = "none")
f_pA <- fisher.method(pVals_Adult, p.corr = "none")
head(f_pA)
head(f_pJ)

#I think you should be using var.test(candidate_effectsize, control_effectsize)

#It would return the ratio of effect sizes, a confidence interval (most useful for biological interpretation, I think) and a p value

#Here is some dummy data, where the ratio of sds is 3 (variances ratio 9). You can see the true value is in the CI, and p is off the scale.

candidate_e <- rnorm(1000, sd=1.1)
control_e <- rnorm(1000, sd=1)
var.test(candidate_e, control_e)
# 
# F test to compare two variances
# 
# data:  candidate_e and control_e
# F = 8.6874, num df = 999, denom df = 999, p-value < 2.2e-16 alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   7.673547 9.835149
# sample estimates:
#   ratio of variances 
# 8.687374 


# pvals <- rep(0,100)
# 
# for (j in 1:100) pvals <-  var.test(gebv_effects[,j], control_effects[,j], alternative = 'greater')$p.value
# 
# library(metaseqR)
# fisher.method(pvalues, p.corr = "none") 

