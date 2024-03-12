#Load required modules
library(Rcpp)
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

#First, we should filter the dataset to remove all the GEBV loci and,
#all loci likely to be in LD with them (i.e. mask out those adjacent areas of the genome)
#- because we want to make use of the GENOME WIDE relatedness in predicting GEBV,
#then comparing observed GEBV with predictions depending on the generation (young old). 

#Get vcf file with maf of at least 0.1
#VCF <- read.vcfR("/data/SBCS-BuggsLab/MardenPark_Ash/08_Filtered_VCF/LD_Filtered_VCF_Files/maf01.pass3.miss25.snps.only.LD.vcf")
VCF <- read.vcfR("C:/Users/carey/Documents/MardonPark/Data/selected_snps.maf01.vcf")
#VCF <- read.vcfR("C:/Users/carey/Documents/MardonPark/Data/maf4.pass3.miss25.snps.only.LD.vcf")
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

#write.csv(GT, file = "maf01.pass3.miss25.snps.only.LD.gt")

#Remove rows with NA
#GT <- GT[which(!(is.na(apply(GT, 1, sum)))),]
GT[is.na(GT)] <- as.numeric(0)

#Lets divide up the genetic data (a matrix of -1, 0, 1)
GT <- GT - 1

#into rows, Adult (rows containing adults) Young (rows containing young).
pheno <- read.csv("C:/Users/carey/Documents/MardonPark/Data/mp_mastersheet.10.03.2022_tidy_gebv.csv")
#Use only invividuals with genotype
pheno <- pheno[which(pheno$Label %in% colnames(GT)),]
adults <- pheno[which(pheno$Type =="Adult" & !is.na(pheno$PercentScore_2019)),]
GTA <- GT[,which(colnames(GT) %in% adults$Label)]
juvs <- pheno[which(pheno$Type =="Juvenile" & !is.na(pheno$Score_2019)),]
GTJ <- GT[,which(colnames(GT) %in% juvs$Label)]

#Order the df and matrix
GTA <- GTA[, order(colnames(GTA))]
adults <- adults[order(adults$Label),]
GTJ <- GTJ[, order(colnames(GTJ))]
juvs <- juvs[order(juvs$Label),]

ETA_J<-list( list(~factor(Sample),
                data=juvs,model='FIXED'),
           list(X=t(GTJ), model='BayesB'))

ETA_J<-list(list(X=t(GTJ), model='BayesB'))

fm_j <- BGLR(y = juvs$Score_2019, ETA = ETA_J, nIter=20000, burnIn=2000, response_type = "ordinal")

#1# Estimated Marker Effects & posterior SDs
bHat<- fm_j$ETA[[1]]$varB
SD.bHat<- fm_j$ETA[[1]]$SD.b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')
#2# Predictions
# Total prediction
yHat<-fm_j$yHat
plot(y = yHat, x =juvs$Score_2019,xlab='Observed',ylab='Predicted',col=2); abline(a=0,b=1,col=4,lwd=2)
# Just the genomic part
gHat<-t(GTJ)%*%fm_j$ETA[[1]]$b
plot(y = gHat, x =juvs$Score_2019, xlab='Phenotype',
     ylab='Predicted Genomic Value',col=2)
#3# Godness of fit and related statistics
fm_j$fit
fm_j$varE # compare to var(y)
#4# Trace plots
# Residual variance
varE_j<-scan('varE.dat')
plot(varE_j,type='o',col=2,cex=.5,ylab=expression(var[e]));
abline(h=fm_j$varE,col=4,lwd=2);
abline(v=fm_j$burnIn/fm_j$thin,col=4)
# lambda (regularization parameter of the Bayesian Lasso)
lambda<-scan('ETA_3_lambda.dat')
plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda));
abline(h=fm_j$ETA[[1]]$lambda,col=4,lwd=2);
abline(v=fm_j$burnIn/fm_j$thin,col=4)

#Heritability estimate
# From variance companats 
varU=fm_j$ETA[[1]]$varB
varE=fm_j$varE
h2_j=varU/(varU+varE)
plot(h2_j,type='o',cex=.5,col=4);abline(h=c(h20,mean(h2[-c(1:200)])),lty=2,col=c(1,2),lwd=2)

#Leave one out training and prediction
gebvj <- c()
for (j in 61:nrow(juvs)){
  ETA <- list(list(X=t(GTJ[,-j]), model='BayesB'))
  fm <- BGLR(y = juvs$Score_2019[-j], ETA = ETA, nIter=20000, burnIn=2000, response_type = "ordinal")
  gebvj[j] <- t(GTJ[,j]) %*% fm$ETA[[1]]$b
  write(gebvj, file = "gebv.j.gblr.loo.2")
}
summary(lm(gebvj~juvs$Score_2019[1:length(gebvj)]))

########################################################
#ADULTS

#Leave one out training and prediction
gebva <- c()
for (a in 25:nrow(adults)){
  ETA <- list(list(X=t(GTA[,-a]), model='BayesB'))
  fm <- BGLR(y = adults$PercentScore_2019[-a], ETA = ETA, nIter=20000, burnIn=2000)
  gebva[a] <- t(GTA[,a]) %*% fm$ETA[[1]]$b
}
summary(lm(gebva~adults$PercentScore_2019[1:length(gebva)]))

ETA_A<-list(list(X=t(GTA), model='BayesB'))

fm_a <- BGLR(y = adults$PercentScore_2019, ETA = ETA_A, nIter=20000, burnIn=2000)
#Heritability estimate
# From variance companats 
varU=fm_a$ETA[[1]]$varB
varE=fm_a$varE
h2_a=varU/(varU+varE)
plot(h2_a,type='o',cex=.5,col=4);abline(h=c(h20,mean(h2[-c(1:200)])),lty=2,col=c(1,2),lwd=2)
# From variance a 
B=fm_j$ETA[[1]]$b
h2_new=rep(NA,nrow(B))
varU_new=h2_new
varE_new=h2_new
for(i in 1:length(h2_new)){
  u=X%*%B[i,]	
  varU_new[i]=var(u)
  varE_new[i]=var(y-u)
  h2_new[i]=varU_new[i]/(varU_new[i]+varE_new[i])
}
plot(h2_new,type='o',cex=.5,col=4);abline(h=c(h2_0,mean(h2_new)),lty=2,col=c(1,2),lwd=2)
mean(varU)
mean(varU+varE)
mean(varU_new+varE_new)

#1# Estimated Marker Effects & posterior SDs
bHat<- fm_a$ETA[[1]]$b
SD.bHat<- fm_a$ETA[[1]]$SD.b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')
#2# Predictions
# Total prediction
yHat<-fm_a$yHat
tmp<-range(c(y,yHat))
plot(y = yHat, x =adults$PercentScore_2019,xlab='Observed',ylab='Predicted',col=2); abline(a=0,b=1,col=4,lwd=2)
# Just the genomic part
gHat<-t(GTA)%*%fm_a$ETA[[1]]$b
plot(y = gHat, x =adults$PercentScore_2019, xlab='Phenotype',
     ylab='Predicted Genomic Value',col=2); abline(a=0,b=1,col=4,lwd=2)
#3# Godness of fit and related statistics
fm_a$fit
fm_a$varE # compare to var(y)
#4# Trace plots
list.files()
# Residual variance
varE<-scan('varE.dat')
plot(varE,type='o',col=2,cex=.5,ylab=expression(var[e]));
abline(h=fm_a$varE,col=4,lwd=2);
abline(v=fm_a$burnIn/fm_a$thin,col=4)
# lambda (regularization parameter of the Bayesian Lasso)
lambda<-scan('ETA_3_lambda.dat')
plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda));
abline(h=fm_a$ETA[[1]]$lambda,col=4,lwd=2);
abline(v=fm_a$burnIn/fm_a$thin,col=4)

#Use none GEBV snps
VCF <- read.vcfR("C:/Users/carey/Documents/MardonPark/Data/maf4.pass3.miss25.snps.only.LD.vcf")
#Convert to gt matrix
GT2 <- convert_GT(extract.gt(VCF, element = "GT"))
rm(VCF)
#Remove rows with NA
class(GT2) <- "numeric"
#Fix names
split_id <-
  sapply(row.names(GT2), function(x)
    str_split(x, pattern = "_"))
rownames(GT2) <- paste(sapply(split_id, "[[", 1), sapply(split_id, "[[", 2), sep = "_")
split_id_col <-
  sapply(colnames(GT2), function(x)
    str_split(x, pattern = "_"))
colnames(GT2) <- sapply(split_id_col, "[[", 1)

#Remove rows with NA
#GT <- GT[which(!(is.na(apply(GT, 1, sum)))),]
GT2[is.na(GT2)] <- as.numeric(0)
#Lets divide up the genetic data (a matrix of -1, 0, 1)
GT2 <- GT2 - 1
#Seperate adults and juveniles
GTA2 <- GT2[,which(colnames(GT2) %in% adults$Label)]
GTJ2 <- GT2[,which(colnames(GT2) %in% juvs$Label)]

#Order the df and matrix
GTA2 <- GTA2[, order(colnames(GTA))]
GTJ2 <- GTJ2[, order(colnames(GTJ))]

#Draw random sample - repeat 20 times
h2_j2 <- list()
for (i in 1:20){
subset <- sample(nrow(GTJ2),size=7985,replace=FALSE)
ETA_J2<-list(list(X=t(GTJ2[subset,]), model='BayesB'))
fm_j2 <- BGLR(y = juvs$Score_2019, ETA = ETA_J2, nIter=20000, burnIn=2000, response_type = "ordinal")
varU=fm_j2$ETA[[1]]$varB
varE=fm_j2$varE
h2_j2[i]=mean(varU/(varU+varE))
}


#Draw random sample - repeat 20 times
h2_a2 <- list()
for (i in 1:20){
  subset <- sample(nrow(GTA2),size=7985,replace=FALSE)
  ETA_A2<-list(list(X=t(GTA2[subset,]), model='BayesB'))
  fm_a2 <- BGLR(y = adults$PercentScore_2019, ETA = ETA_A2, nIter=20000, burnIn=2000)
  varU=fm_j2$ETA[[1]]$varB
  varE=fm_j2$varE
  h2_a2[i]=mean(varU/(varU+varE))
}
mean(unlist(h2_a2))


