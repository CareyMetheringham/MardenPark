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
VCF <- read.vcfR("selected_snps.maf01.vcf")
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

VCF2 <- read.vcfR("maf01.pass3.miss25.snps.only.LD.vcf",
                  nrows = 100000)
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
pheno <- read.csv("phenotypes.csv")
#Use only invividuals with genotype
pheno <- pheno[which(pheno$Label %in% colnames(GT2)),]
adults <- pheno[which(pheno$Type =="Adult" & !is.na(pheno$PercentScore_2019)),]
adults <- adults[order(adults$Label),]
#GTA <- GTF[,which(colnames(GTF) %in% adults$Label)]
juvs <- pheno[which(pheno$Type =="Juvenile" & !is.na(pheno$Score_2019)),]
juvs <- juvs[order(juvs$Label),]

#Define the number of sites to use
site_num <- round(nrow(GT))/4

# Get a list of effect sizes and find the sites with the largest absolute effect sizes
EES <- read.csv("effect_sizes.csv", sep="")
large_effect <- EES[order((EES$EES.MIA^2), decreasing = T),][1:site_num,]
large_effect_snps <- paste(large_effect$CONTIG, large_effect$SNP, sep = "_")

#JUVENILES ONLY

#Pull those sites from the GT table
GTF <- GT[which(rownames(GT) %in% large_effect_snps), which(colnames(GT) %in% juvs$Label)]
GTF[is.na(GTF)] <- as.numeric(0)

snp_effect_tableA <- data.frame(rep(0, nrow(GTF)*2))
snp_effect_tableJ <- data.frame(rep(0, nrow(GTF)*2))
for (i in 1:100){

#Pull an equal number of sites at random from the non GEBV sites
GTR <- GT2[sample(nrow(GT2), nrow(GTF)), which(colnames(GT2) %in% juvs$Label)]
GTR[is.na(GTR)] <- as.numeric(0)

#Combine the two sets of SNPs
GT3 <- rbind(GTR, GTF)
#Order the GT matrix by sample
GT3 <- GT3[order(rownames(GT3)), order(colnames(GT3))]
class(GT3) <- "numeric"

#Create the ETA object
ETAJ <-list(MRK=list(X=t(GT3), model='BayesB'))
#Fit the model
fmj <- BGLR(y = juvs$Score_2019, ETA = ETAJ, nIter=50000, burnIn=2000, verbose = F)

snp_effect_tableJ[,i] <- fmj$ETA[[1]]$b

#ADULTS ONLY

#Pull those sites from the GT table
GTFA <- GT[which(rownames(GT) %in% large_effect_snps), which(colnames(GT) %in% adults$Label)]
GTFA[is.na(GTFA)] <- as.numeric(0)

#Pull an equal number of sites at random from the non GEBV sites
GTRA <- GT2[sample(nrow(GT2), nrow(GTF)), which(colnames(GT2) %in% adults$Label)]
GTRA[is.na(GTRA)] <- as.numeric(0)

#Combine the two sets of SNPs
GT3A <- rbind(GTRA, GTFA)
#Order the GT matrix by sample
GT3A <- GT3A[order(rownames(GT3A)), order(colnames(GT3A))]
class(GT3A) <- "numeric"

#Create the ETA object
ETAA <-list(MRK=list(X=t(GT3A), model='BayesB'))
#Fit the model
fma <- BGLR(y = adults$PercentScore_2019, ETA = ETAA, nIter=50000, burnIn=2000, verbose = F)

snp_effect_tableA[,i] <- fma$ETA[[1]]$b

}

write.csv(snp_effect_tableA, file = "snp_effect_tableA.csv")
write.csv(snp_effect_tableJ, file = "snp_effect_tableJ.csv")
