setwd("~/MardonPark/Data")

#Load required modules
library(vcfR)
library(BGLR)
library(dplyr)
library(stringr)

convert_GT <- function(GT, convert_na = F) {
  #replace values with numeric
  GT[GT == "0/0"] <- as.numeric(0)
  GT[GT == "1/0"] <- as.numeric(1)
  GT[GT == "0/1"] <- as.numeric(1)
  GT[GT == "0/2"] <- as.numeric(1)
  GT[GT == "0/3"] <- as.numeric(1)
  GT[GT == "1/1"] <- as.numeric(2)
  GT[GT == "1/2"] <- as.numeric(2)
  GT[GT == "2/2"] <- as.numeric(2)
  GT[GT == "0|0"] <- as.numeric(0)
  GT[GT == "1|0"] <- as.numeric(1)
  GT[GT == "0|1"] <- as.numeric(1)
  GT[GT == "0|2"] <- as.numeric(1)
  GT[GT == "0|3"] <- as.numeric(1)
  GT[GT == "1|1"] <- as.numeric(2)
  GT[GT == "1|2"] <- as.numeric(2)
  GT[GT == "2|2"] <- as.numeric(2)
  GT[GT == "3|3"] <- as.numeric(2)
  GT[GT == "3/3"] <- as.numeric(2)
  GT[GT == "1/3"] <- as.numeric(2)
  GT[GT == "2/3"] <- as.numeric(2)
  if (convert_na == T) {
    GT[is.na(GT)] <- as.numeric(-9)
  }
  return(GT)
}

#Get vcf file
VCF <- read.vcfR("~/MardonPark/Data/selected_snps.maf01.vcf")
#Convert to gt matrix
GT <- convert_GT(extract.gt(VCF, element = "GT"))
rm(VCF)
#Remove rows with NA
class(GT) <- "numeric"
#Fix names
split_id <-
  sapply(row.names(GT), function(x)
    str_split(x, pattern = "_"))
rownames(GT) <- paste(sapply(split_id, "[[", 1), sapply(split_id, "[[", 2), sep = "_")
split_id_col <-
  sapply(colnames(GT), function(x)
    str_split(x, pattern = "_"))
colnames(GT) <- sapply(split_id_col, "[[", 1)

VCF2 <- read.vcfR("~/MardonPark/Data/maf4.pass3.miss25.snps.only.LD.vcf")
#Convert to gt matrix
GT2 <- convert_GT(extract.gt(VCF2, element = "GT"))
rm(VCF2)
#Fix names
split_id <-
  sapply(row.names(GT2), function(x)
    str_split(x, pattern = "_"))
rownames(GT2) <- paste(sapply(split_id, "[[", 1), sapply(split_id, "[[", 2), sep = "_")
split_id_col <-
  sapply(colnames(GT2), function(x)
    str_split(x, pattern = "_"))
colnames(GT2) <- sapply(split_id_col, "[[", 1)

#into rows, Adult (rows containing adults) Young (rows containing young).
pheno <- read.csv("~/MardonPark/Data/mp_mastersheet.10.03.2022_tidy_gebv.csv")
#Use only invividuals with genotype
pheno <- pheno[which(pheno$Label %in% colnames(GT2)),]
adults <- pheno[which(pheno$Type =="Adult" & !is.na(pheno$PercentScore_2019)),]
adults <- adults[order(adults$Label),]
#GTA <- GTF[,which(colnames(GTF) %in% adults$Label)]
juvs <- pheno[which(pheno$Type =="Juvenile" & !is.na(pheno$Score_2019)),]
juvs <- juvs[order(juvs$Label),]

#Define the number of sites to use
site_num <- 7985/2

# Get a list of effect sizes and find the sites with the largest absolute effect sizes
EES <- read.csv("~/MardonPark/Data/new_ees.table_10000", sep="")
large_effect <- EES[order((EES$EES.MIA^2), decreasing = T),][1:site_num,]
large_effect_snps <- paste(large_effect$CONTIG, large_effect$SNP, sep = "_")

#JUVENILES ONLY

#Pull those sites from the GT table
GTF <- GT[which(rownames(GT) %in% large_effect_snps), which(colnames(GT) %in% juvs$Label)]
GTF[is.na(GTF)] <- as.numeric(0)

#Pull an equal number of sites at random from the non GEBV sites
GTR <- GT2[sample(nrow(GT2), nrow(GTF)), which(colnames(GT2) %in% juvs$Label)]
GTR[is.na(GTR)] <- as.numeric(0)

#Combine the two sets of SNPs
GT3 <- rbind(GTR, GTF)
#Check size of matrix
dim(GT3)
#Order the GT matrix by sample
GT3 <- GT3[, order(colnames(GT3))]
class(GT3) <- "numeric"

#Create the ETA object
ETAJ <-list(MRK=list(X=t(GT3), model='BayesB'))

fmj <- BGLR(y = juvs$Score_2019, ETA = ETAJ, nIter=20000, burnIn=2000, response_type = "ordinal")

#1# Estimated Marker Effects & posterior SDs
bHat<- fmj$ETA[[1]]$b
SD.bHat<- fmj$ETA[[1]]$SD.b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')
abline(v=nrow(GTF),col=2)
#2# Predictions
# Total prediction
yHat<-fmj$yHat
plot(y = yHat, x =juvs$Score_2019,xlab='Observed',ylab='Predicted',col=2); abline(a=0,b=1,col=4,lwd=2)
# Just the genomic part
gHat<-t(GT3)%*%fmj$ETA[[1]]$b
plot(y = gHat, x =juvs$Score_2019, xlab='Phenotype',
     ylab='Predicted Genomic Value',col=2); abline(a=0,b=1,col=4,lwd=2)
#3# Godness of fit and related statistics
fmj$fit

#Are the GEBV sites more likely to have effect assigned to them?
table(fmj$ETA[[1]]$b[1:nrow(GTF)] >0.0000001)
table(fmj$ETA[[1]]$b[nrow(GTF):nrow(GT3)] >0.0000001)
wilcox.test(fmj$ETA[[1]]$b[1:nrow(GTF)]^2, fmj$ETA[[1]]$b[nrow(GTF):nrow(GT3)]^2)

#ADULTS ONLY

#Pull those sites from the GT table
GTFA <- GT[which(rownames(GT) %in% large_effect_snps), which(colnames(GT) %in% adults$Label)]
GTFA[is.na(GTFA)] <- as.numeric(0)

#Pull an equal number of sites at random from the non GEBV sites
GTRA <- GT2[sample(nrow(GT2), nrow(GTF)), which(colnames(GT2) %in% adults$Label)]
GTRA[is.na(GTRA)] <- as.numeric(0)

#Combine the two sets of SNPs
GT3A <- rbind(GTRA, GTFA)
#Check size of matrix
dim(GT3A)
#Order the GT matrix by sample
GT3A <- GT3A[, order(colnames(GT3A))]
class(GT3A) <- "numeric"

#Create the ETA object
ETAA <-list(MRK=list(X=t(GT3A), model='BayesB'))

fma <- BGLR(y = adults$PercentScore_2019, ETA = ETAA, nIter=20000, burnIn=2000)

#1# Estimated Marker Effects & posterior SDs
bHat<- fma$ETA[[1]]$b
SD.bHat<- fma$ETA[[1]]$SD.b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')
abline(v=nrow(GTF),col=2)
#2# Predictions
# Total prediction
yHata<-fma$yHat
plot(y = yHata, x =adults$PercentScore_2019,xlab='Observed',ylab='Predicted',col=2); abline(a=0,b=1,col=4,lwd=2)
# Just the genomic part
gHata<-t(GT3A)%*%fma$ETA[[1]]$b
plot(y = gHata, x =adults$PercentScore_2019, xlab='Phenotype',
     ylab='Predicted Genomic Value',col=2); abline(a=0,b=1,col=4,lwd=2)
#3# Godness of fit and related statistics
fma$fit
fma$varEa # compare to var(y)
# Residual variance
varEa<-scan('varE.dat')
plot(varEa,type='o',col=2,cex=.5,ylab=expression(var[e]));
abline(h=fma$varEa,col=4,lwd=2);
abline(v=fma$burnIn/fma$thin,col=4)

#Are the GEBV sites more likely to have effect assigned to them?
table(fma$ETA[[1]]$b[1:nrow(GTF)] >0.0000001)
table(fma$ETA[[1]]$b[nrow(GTF):nrow(GT3)] >0.0000001)
wilcox.test(fma$ETA[[1]]$b[1:nrow(GTF)]^2, fma$ETA[[1]]$b[nrow(GTF):nrow(GT3)]^2)