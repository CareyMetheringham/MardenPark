library(rrBLUP)
library(Hmisc)
library(cowplot)

#Import csv files with data on the GEBV sites and unlinked sites
gebv_sites <- read.csv(file = "gebv_sites.csv")
unlinked_sites <- read.csv(file = "unlinked_sites.csv")
#Pull out the adults and juveniles
adult_unlinked <- unlinked_sites[unlinked_sites$Age=="Adult",]
juv_unlinked <- unlinked_sites[unlinked_sites$Age=="Juv",]
#Extract the genotype matrices
GTA <- adult_unlinked[, 4:ncol(adult_unlinked)]
GTJ <- juv_unlinked[, 4:ncol(juv_unlinked)]

#Select 5,000 sites at random
gt_sample <- sample(colnames(GTA), 5000, replace = F)
GTA <- GTA[, which(colnames(GTA) %in% gt_sample)]
GTJ <- GTJ[, which(colnames(GTJ) %in% gt_sample)]

#Train the model in adults
A_from_unlinked <-
  mixed.solve(
    adult_unlinked[,2],
    Z = GTA,
    K = NULL,
    SE = FALSE,
    return.Hinv = FALSE
  )
#Predict GEBV in both adults and juveniles
EST_A_from_unlinked <- as.vector(as.matrix(GTA) %*% A_from_unlinked$u)
lm(adult_unlinked$GEBV ~ EST_A_from_unlinked)$coefficients
EST_J_from_unlinkedA <- as.vector(as.matrix(GTJ) %*% A_from_unlinked$u)
lm(juv_unlinked$GEBV ~ EST_J_from_unlinkedA)$coefficients

#Train the model in juveniles
J_from_unlinked  <-
  mixed.solve(
    juv_unlinked[,2],
    Z = GTJ,
    K = NULL,
    SE = FALSE,
    return.Hinv = FALSE
  )
#Predict GEBV in both adults and juveniles
EST_J_from_unlinked <- as.vector(as.matrix(GTJ) %*% J_from_unlinked$u)
lm(juv_unlinked$GEBV ~ EST_J_from_unlinked)$coefficients
EST_A_from_unlinkedJ <- as.vector(as.matrix(GTA) %*% J_from_unlinked$u)
lm(adult_unlinked$GEBV ~ EST_A_from_unlinkedJ)$coefficients

#Leave one out models
loo_est_A <- c()
for (i in 1:nrow(adult_unlinked)){
  BLUPA <-
    mixed.solve(
      adult_unlinked[-i,2],
      Z = GTA[-i,],
      K = NULL,
      SE = FALSE,
      return.Hinv = FALSE
    )
  loo_est_A[i] <- as.vector(as.matrix(GTA[i,]) %*% BLUPA$u)
}

#Use LOO to validate the results
#Leave one out for Juveniles
loo_est_J <- c()
for (i in 1:nrow(juv_unlinked)){
  BLUPJ <-
    mixed.solve(
      juv_unlinked[-i,2],
      Z = GTJ[-i,],
      K = NULL,
      SE = FALSE,
      return.Hinv = FALSE
    )
  loo_est_J[i] <- as.vector(as.matrix(GTJ[i,]) %*% BLUPJ$u)
}

#Format output into table
juv_gebv_table <- data.frame(juv_unlinked$GEBV,EST_J_from_unlinkedA,EST_J_from_unlinked,loo_est_J)
adult_gebv_table <- data.frame(adult_unlinked$GEBV, EST_A_from_unlinked, EST_A_from_unlinkedJ, loo_est_A)
colnames(juv_gebv_table) <- colnames(adult_gebv_table) <- c("Previous", "PredA", "PredJ", "PredLOO")

#Check for significant difference in intercept
#Check the deviation from the 1:1 line
m0 <-  lm(juv_gebv_table$Previous ~ 0, offset = juv_gebv_table$PredA)
m1 <- lm(juv_gebv_table$Previous ~ juv_gebv_table$PredA)
anova(m0, m1)

#Plot output
p1 <- ggplot(data = juv_gebv_table, aes(x=PredLOO,y=Previous))+
  geom_point(alpha = 0.5)+
  geom_abline(slope = 1,intercept = 0, colour="red", linetype = "dashed")+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()+
  xlab("GEBV Predicted from Related Juveniles")+
  ylab("Juvenile GEBV")+
  ggtitle("A")+
  ylim(-0.6,1.5)+
  xlim(-0.6,0.7)

p2 <- ggplot(data = juv_gebv_table, aes(x=PredA,y=Previous))+
  geom_point(alpha = 0.5)+
  geom_abline(slope = 1,intercept = 0, colour="red", linetype = "dashed")+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()+
  xlab("GEBV Predicted from Related Adults")+
  ylab("Juvenile GEBV")+
  ggtitle("B")+
  ylim(-0.6,1.5)+
  xlim(-0.6,0.7)


plots <- ggdraw() +
  draw_plot(p1, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(p2, x = .5, y = 0, width = .5, height = 1) +

tiff("mp_ex_fig3.tiff", units="mm", width=180, height=100, res=300)
plots
dev.off()
