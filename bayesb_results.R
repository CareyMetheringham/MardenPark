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
