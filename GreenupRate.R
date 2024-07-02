#Model R.A.Nichols
#Plotting C.L.Metheringham

library(brms)
library(viridis)
library(ggplot2)
library(cowplot)

# set working directory to the current location
# replace "." with the path to the location of "gebv_model_df.csv"
# if it is elsewhere
setwd(".")

# read in the greening data 
dat1 <- read.csv("~/Documents/Ash/gebv_model_df.csv")

# add small amount to each date so first value (0) can be logged
dat1$time <- dat1$time + 0.00001

# Set random seeds
set.seed(42)

# Set broad priors
prior1 <- c(prior(normal(2.0, 0.03), nlpar = "s"), # note slope = exp(s)
            prior(normal(0.338, 0.006), nlpar = "l"),
            prior(normal(0.46, 0.03), nlpar = "u"),
            prior(normal(0.9,0.05), nlpar = "e")
)

# fit a single curve to the data
fit1 <- brm(
  bf(gcc ~ l + (u - l) / (1 + exp(-exp(s) * (log(time) - log(e)))),
     s ~ 1, l ~ 1, u ~ 1, e ~ 1,
     nl = TRUE),
  data = dat1, family = gaussian(),
  prior = prior1,
  control = list(adapt_delta = 0.9),
  chains = 10,
  cores = 10
)

# now fit a model with rate of greening differing between trees as a
# function of gebv, and a random effect giving the June score for each tree
prior2 <- c(prior(normal(2, 0.05), nlpar = "s"), # note slope = exp(s)
            prior(normal(0.2, 0.2), coef = "gebv", nlpar = "s"),
            
            prior(normal(0.338, 0.006), lb = 0.01, nlpar = "l"),
            
            prior(normal(0.43, 0.03), ub = 1, nlpar = "u"),
            prior(student_t(3, 0.03, 0.04), class = "sd", group = "treeID", nlpar = "u"),
            
            prior(normal(0.8,0.05), nlpar = "e")
)

fit2 <- brm(
  bf(gcc ~ l + (u - l) / (1 + exp(-exp(s) * (log(time) - log(e)))),
     s ~ 1 + gebv, 
     l ~ 1 , 
     u ~ 1 + (1|treeID),  
     e ~ 1 , 
     nl = TRUE),
  data = dat1, family = gaussian(),
  prior = prior2,
  control = list(adapt_delta = 0.95),
  chains = 10,
  cores = 10,
)

# choose small value to add to each date so the different
# trees can be visualized clearly 
# (jittered values used for display only, not in the analysis)
ntrees <- length(unique(dat1$treeID))
jitter <- runif(nrow(dat1),
                min = -sd(dat1$time)/20,
                max = sd(dat1$time)/20
)

# ggplot
# Define the color palette
cols <- viridis(ntrees)

# Prepare data
treenames <- unique(dat1$treeID)
indiv_gebv <- rep(0, ntrees)
for (i in 1:ntrees) indiv_gebv[i] <- mean(dat1$gebv[dat1$treeID == treenames[i]])

xvals <- seq(0,1,length.out = 1000)
lines_data <- data.frame()
for (i in 1:ntrees) {
    fdat <- data.frame(time = xvals,
                       treeID = rep(treenames[i], length(xvals)),
                       gebv = rep(indiv_gebv[i], length(xvals)))
    fitted_values <- predict(fit2, newdata = fdat)[,1]
    lines_data <- rbind(lines_data, data.frame(time = xvals, fitted = fitted_values, treeID = treenames[i], col = cols[i]))
}

# Create plot
p1 <- ggplot(dat1, aes(x = time + jitter, y = gcc)) +
    geom_point(aes(color = as.factor(treeID)), size = 2) +
    scale_color_manual(values = cols) +
    geom_line(data = lines_data, aes(x = time, y = fitted, color = as.factor(treeID)), alpha = 1) +
    labs(x = "Date (jittered)", y = "Greening score (gcc)",
         title = "Rate of spring greening with GEBV (yellow ... purple)") +
    theme_minimal()+
    theme(legend.position = "none")
# Add fitted line to the first plot
xvals <- seq(0, 1, length.out = 1000)
fdat <- data.frame(time = xvals)
fitted_values <- predict(fit1, newdata = fdat)[,1]
p1 <- p1 + geom_line(data = data.frame(time = xvals, fitted = fitted_values),
                     aes(x = time, y = fitted))

# Calculate quantiles
quantiles <- quantile(dat1$gebv, probs = c(0.25, 0.50, 0.75))

# Split data into quantile groups
group_1 <- dat1[dat1$gebv <= quantiles[1], ]
group_2 <- dat1[dat1$gebv >= quantiles[1] & dat1$gebv < quantiles[2], ]
group_3 <- dat1[dat1$gebv >= quantiles[2] & dat1$gebv < quantiles[3], ]
group_4 <- dat1[dat1$gebv >= quantiles[3], ]

# Calculate mean GEBV for each group
mean_gebv_1 <- mean(group_1$gebv)
mean_gebv_2 <- mean(group_2$gebv)
mean_gebv_3 <- mean(group_3$gebv)
mean_gebv_4 <- mean(group_4$gebv)

# Select a treeID from each group to use for mean line prediction
treeID_1 <- unique(group_1$treeID)[1]
treeID_2 <- unique(group_2$treeID)[1]
treeID_3 <- unique(group_3$treeID)[1]
treeID_4 <- unique(group_4$treeID)[1]

# Create data for mean lines
xvals <- seq(0, 1, length.out = 100)
mean_fdat_1 <- data.frame(time = xvals, treeID = rep(treeID_1, 100), gebv = rep(mean_gebv_1, 100))
mean_fdat_2 <- data.frame(time = xvals, treeID = rep(treeID_2, 100), gebv = rep(mean_gebv_2, 100))
mean_fdat_3 <- data.frame(time = xvals, treeID = rep(treeID_3, 100), gebv = rep(mean_gebv_3, 100))
mean_fdat_4 <- data.frame(time = xvals, treeID = rep(treeID_4, 100), gebv = rep(mean_gebv_4, 100))

# Predictions for mean lines
mean_fdat_1$pred <- predict(fit2, newdata = mean_fdat_1)[, 1]
mean_fdat_2$pred <- predict(fit2, newdata = mean_fdat_2)[, 1]
mean_fdat_3$pred <- predict(fit2, newdata = mean_fdat_3)[, 1]
mean_fdat_4$pred <- predict(fit2, newdata = mean_fdat_4)[, 1]

# Combine mean line data
mean_lines_data <- rbind(mean_fdat_1, mean_fdat_2, mean_fdat_3, mean_fdat_4)

# Create the plot
p2 <- ggplot(dat1, aes(x = time + jitter, y = gcc)) +
    geom_point(aes(color = as.factor(treeID)), size = 2) +
    scale_color_manual(values = cols) +
    geom_line(data = mean_lines_data[mean_lines_data$gebv == mean_gebv_1,], aes(x = time, y = pred), color = "#E5E419FF", lwd = 1, lty = 2) +
    geom_line(data = mean_lines_data[mean_lines_data$gebv == mean_gebv_2,], aes(x = time, y = pred), color = "#53C569FF", lwd = 1, lty = 2) +
    geom_line(data = mean_lines_data[mean_lines_data$gebv == mean_gebv_3,], aes(x = time, y = pred), color = "#365D8DFF", lwd = 1, lty = 2) +
    geom_line(data = mean_lines_data[mean_lines_data$gebv == mean_gebv_4,], aes(x = time, y = pred), color = "#440154FF", lwd = 1, lty = 2) +
    annotate("text", x = 0.93, y = 0.42, label = "Q1", colour = "#E5E419FF", size = 4)+
    annotate("text", x = 0.91, y = 0.403, label = "Q2", colour = "#53C569FF", size = 4)+
    annotate("text", x = 0.88, y = 0.38, label = "Q3", colour = "#365D8DFF", size = 4)+
    annotate("text", x = 0.95, y = 0.45, label = "Q4", colour = "#440154FF", size = 4)+
    labs(x = "Date (jittered)", y = "Greening score (gcc)",
         title = "Rate of spring greening per GEBV quantile (yellow ... purple)") +
    theme_minimal() +
    theme(legend.position = "none")

# Combine the two plots
combined_plot <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 1)

#Save output to file
tiff("~/Documents/Ash/Figures/greenup.tiff", units="mm", width=180, height=200, res=300)
combined_plot
dev.off()
png("~/Documents/Ash/Figures/greenup.png", units="mm", width=180, height=200, res=300)
combined_plot
dev.off()
jpeg("~/Documents/Ash/Figures/greenup.jpg", units="mm", width=180, height=200, res=300)
combined_plot
dev.off()
