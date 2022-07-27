library(rrBLUP)
library(ggplot2)
library(cowplot)

#Detecting shifts in breeding values between adults and juveniles 
#Accounting for ancestry of juveniles
#Based on methods and code by Richard A. Nichols

#Import effect sizes - Supplementary Data 2
es_file <- read.csv("MP_effects_MIA_and_MAA.csv")
es <- with(es_file, EES.MIA - EES.MAA)

#Import csv files with data on the GEBV sites and unlinked sites
#- Supplementary Data 3
gebv_sites <- read.csv(file = "gebv_sites.csv")
#Extract the genotype matrices
GTAgebv <- gebv_sites[gebv_sites$Age=="Adult", 4:ncol(gebv_sites)]
GTJgebv <- gebv_sites[gebv_sites$Age=="Juv", 4:ncol(gebv_sites)]

#Calculate GEBV scores 
new_gebv <- (as.matrix(gebv_sites[,-(1:3)])%*%es) + 1
#Check that they are the same as listed
gebv_A <- new_gebv[gebv_sites$Age=="Adult"]
gebv_J <- new_gebv[gebv_sites$Age=="Juv"]

#Plot the raw difference in gebv between adults and juveniles
gebv_A_and_J <- data.frame(Age = gebv_sites$Age, GEBV = new_gebv)
Xlabs <- c("Adults (n=128)", "Juveniles (n = 452)")
fig2A <- ggplot(data = gebv_A_and_J, aes(x = Age, y = GEBV))+
  geom_boxplot(fill = c("#D55E00", "skyblue")) +
  theme_minimal() + xlab("") +
  theme(axis.text=element_text(size=12))+
  ggtitle("A")+
  ylab("GEBV")+
  scale_x_discrete(labels= Xlabs)

#Unlinked sites
unlinked_sites <- read.csv(file = "unlinked_sites.csv")
#Pull out the adults and juveniles
adult_unlinked <- unlinked_sites[unlinked_sites$Age=="Adult",]
juv_unlinked <- unlinked_sites[unlinked_sites$Age=="Juv",]
#Extract the genotype matrices
GTA <- adult_unlinked[, 4:ncol(adult_unlinked)]
GTJ <- juv_unlinked[, 4:ncol(juv_unlinked)]

#Train the model in juveniles
J_from_unlinked  <-
  mixed.solve(
    gebv_J,
    Z = GTJ,
    K = NULL,
    SE = FALSE,
    return.Hinv = FALSE
  )
#Predict GEBV in juveniles
#Overfitting of juveniles from juveniles
EST_J_from_unlinked <- as.vector(as.matrix(GTJ) %*% J_from_unlinked$u)+ rep(J_from_unlinked$beta, length(adult_unlinked[,2]))
#Fit linear model
lm(gebv_J ~ EST_J_from_unlinked)
plot(EST_J_from_unlinked ~ gebv_J)
abline(0,1, col = "red")

#Train the model in adults
A_from_unlinked <-
  mixed.solve(
    gebv_A,
    Z = GTA,
    K = NULL,
    SE = FALSE,
    return.Hinv = FALSE
  )
#Predict GEBV in juveniles
EST_J_from_unlinkedA <- as.vector(as.matrix(GTJ) %*% A_from_unlinked$u) + rep(A_from_unlinked$beta, length(gebv_J))
lm(gebv_J ~ EST_J_from_unlinkedA)$coefficients
#Check for significant difference in intercept

#Method by R.A.N:
#Quantify the shift between adults and juveniles
mod <- lm(juv_unlinked$GEBV ~ 1, offset = EST_J_from_unlinkedA)
summary(mod)
shift <- coefficients(mod)/sd(gebv_J)
shift

#Observed vs predicted
plot(gebv_J~ EST_J_from_unlinkedA)
abline(0,1)
abline(coefficients(mod), 1, col  = 'blue')

#Plot the predictions
juvenile_gebv <- data.frame(GEBV = gebv_J, Predicted = EST_J_from_unlinkedA)
fig2B <- ggplot(data = juvenile_gebv, aes(y = GEBV, x = Predicted))+
  geom_point(alpha = 0.5)+
  theme_minimal()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  geom_abline(slope = 1, intercept = coefficients(mod)[1], colour = "blue")+
  xlab("Juvenile GEBV scores predicted from ancestry")+
  ggtitle("B")

#Plot 2A and 2B together
combined_plot_2 <- ggdraw() +
  draw_plot(fig2A, x = 0.1, y = 0.5, width = 0.8, height = 0.45) +
  draw_plot(fig2B, x = 0.10, y = 0.05, width = 0.8, height = 0.5)
combined_plot_2

tiff("../Figures/marden_park_fig2.tiff", units="mm", width=180, height=200, res=300)
combined_plot_2 
dev.off()
png("../Figures/marden_park_fig2.png", units="mm", width=180, height=200, res=300)
combined_plot_2 
dev.off()
jpeg("../Figures/marden_park_fig2.jpeg", units="mm", width=180, height=200, res=300)
combined_plot_2 
dev.off()

#Compare with shift from neutral expectations (allowing for relatedness)
#to crude difference in mean GEBV score
rawDiff <- mean(gebv_J) - mean(gebv_A)
rawDiff/coefficients(mod)

#What truncation selection is required to give same mean shift?
shiftcalc <- function(cutoff, sd=1)
{integrate(function (x) dnorm(x, sd=sd)*x,
           lower = cutoff,
           upper= 11)$value / pnorm(cutoff, sd=sd, lower.tail = FALSE)
}

#solved by trial and error for a shift of 0.1535
shiftcalc(-1.78,sd=sqrt(2.5))*.4
pnorm(-1.78, sd = sqrt(2.5))

#HERITABILITY

#Assuming 0.4 was heritability 
#If heritability is dropped to 0.2
shiftcalc(-0.86,sd=sqrt(2.5))*.2
pnorm(-0.86, sd = sqrt(2.5))
#This same shift would be equivilent to truncation selection of 29%?

#If heritability is even lower at 0.1
shiftcalc(0.41,sd=sqrt(2.5))*.1
pnorm(0.41, sd = sqrt(2.5))
#Then that would imply 60% truncation!

#SANITY CHECKS ETC

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