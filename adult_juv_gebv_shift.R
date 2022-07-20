library(rrBLUP)

#Detecting shifts in breeding values between adults and juveniles 
#Accounding for relatedness 
#Based on methods and code by Richard A. Nichols

#Import effect sizes
es_file <- read.csv("MP_effects_MIA_and_MAA.csv")
es <- with(es_file, EES.MIA - EES.MAA)

#Import csv files with data on the GEBV sites and unlinked sites
gebv_sites <- read.csv(file = "gebv_sites.csv")
#Extract the genotype matrices
GTAgebv <- gebv_sites[gebv_sites$Age=="Adult", 4:ncol(gebv_sites)]
GTJgebv <- gebv_sites[gebv_sites$Age=="Juv", 4:ncol(gebv_sites)]

#Calculate GEBV scores 
new_gebv <- (as.matrix(gebv_sites[,-(1:3)])%*%es) + 1
#Check that they are the same as listed
plot(new_gebv, gebv_sites$GEBV)

#Unlinked sites
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

#Train the model in juveniles
J_from_unlinked  <-
  mixed.solve(
    juv_unlinked[,2],
    Z = GTJ,
    K = NULL,
    SE = FALSE,
    return.Hinv = FALSE
  )
#Predict GEBV in juveniles
#Overfitting of juveniles from juveniles
EST_J_from_unlinked <- as.vector(as.matrix(GTJ) %*% J_from_unlinked$u)+ rep(J_from_unlinked$beta, length(adult_unlinked[,2]))
#Fit linear model
lm(juv_unlinked$GEBV ~ EST_J_from_unlinked)
plot(EST_J_from_unlinked ~ juv_unlinked$GEBV)
abline(0,1, col = "red")

#Train the model in adults
A_from_unlinked <-
  mixed.solve(
    adult_unlinked[,2],
    Z = GTA,
    K = NULL,
    SE = FALSE,
    return.Hinv = FALSE
  )
#Predict GEBV in juveniles
EST_J_from_unlinkedA <- as.vector(as.matrix(GTJ) %*% A_from_unlinked$u) + rep(A_from_unlinked$beta, length(juv_unlinked[,2]))
lm(juv_unlinked$GEBV ~ EST_J_from_unlinkedA)$coefficients
#Check for significant difference in intercept

#Method by R.A.N:
#Check the deviation from the 1:1 line
mod <- lm(juv_unlinked$GEBV ~ 1, offset = EST_J_from_unlinkedA)
summary(mod)
shift <- coefficients(mod)/sd(juv_unlinked$GEBV)
shift

#Observed vs predicted
plot(EST_J_from_unlinkedA ~ juv_unlinked$GEBV)
abline(0,1)
abline(coefficients(mod), 1, col  = 'blue')


#Compare with shift from neutral expectations (allowing for relatedness)
#to crude difference in mean GEBV score
rawDiff <- mean(juv_gebv_table$Previous) - mean(adult_gebv_table$Previous)
rawDiff/coefficients(mod)

#What truncation selection is required to give same mean shift?
shiftcalc <- function(cutoff, sd=1)
{integrate(function (x) dnorm(x, sd=sd)*x,
           lower = cutoff,
           upper= 11)$value / pnorm(cutoff, sd=sd, lower.tail = FALSE)
}

# shiftcalc(-1.435)
# pnorm(-1.435)
# 
# # now do it for 0.422/0.4 = 1.06
# shiftcalc(-1.78,sd=sqrt(2.5))*.4
# pnorm(-1.78, sd = sqrt(2.5))

# replicate this estimate of selection on raw GEBV score with a random normal sample
randNorm <- rnorm(100000)
# by trial and error truncation of values below 0.13 gives a shift of 0.06
mean(randNorm[randNorm > -0.677])


# now add environmental variation with a heritability of 0.4 (i.e. Ve = 1.5 Va)
randPheno <- randNorm + rnorm(100000, sd = sqrt(1.5))

mean(randNorm[randPheno > -0.33])
pnorm(-0.33, sd = sqrt(2.5))

# Simulate selection on GEBV
# normal distribution to simulate the variation in predicted BV
# Observed sd of juveniles
sdj <- sqrt(0.1425774)
sdj
predBV <- rnorm(700,sd = 0.06)
hist(predBV)









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





#Plot output
p_4B <- ggplot(data = juv_gebv_table, aes(x=PredA,y=Previous))+
  geom_point(alpha = 0.5)+
  geom_abline(slope = 1,intercept = 0, colour="red", linetype = "dashed")+
  geom_smooth(method = "lm", se = FALSE)+
  theme_minimal()+
  xlab('Juveniles\' GEBV scores predicted from their ancestry')+
  ylab("Juvenile GEBV")+
  ggtitle("B")


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
  xlab('Juveniles\' GEBV scores predicted from their ancestry')+
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
