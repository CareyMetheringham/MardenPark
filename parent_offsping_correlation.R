library(ggplot2)
library(weights)
library(cowplot)

#Import data included in Supplementary Table 3
parent_data<- read.csv(file.choose())
colnames(parent_data) <- c("Parent", "DBH", "Percent_2019", "GEBV", "Num_Offsping", "Offspring_Score")

# Test for correlations
cor.test(parent_data$GEBV, parent_data$Num_Offspring)
cor.test(parent_data$Percent_2019, parent_data$Num_Offspring)
# parent phenotype vs juv scores
wtd.cor(parent_data$Offspring_Score, parent_data$Percent_2019, weight = parent_data$Num_Offspring)
# GEBV vs juv scores
cor.test(parent_data$GEBV, parent_data$Offspring_Score)
wtd.cor(parent_data$Offspring_Score, parent_data$GEBV, weight = parent_data$Num_Offspring)

#Plot results
canopy_plot <- ggplot(data = parent_data, aes(x = Percent_2019, y = Offspring_Score, size = Num_Offspring))+
  geom_point(show.legend = FALSE) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(y = "Mean Offspring Score", x = "Canopy Cover of Parent Tree (%)", size = "Number of Offspring")+
theme(legend.position="none")+
  ggtitle("A")

gebv_plot <- ggplot(data = parent_data, aes(x = GEBV, y = Offspring_Score, size = Num_Offspring))+
  geom_point(show.legend = FALSE) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(y = "Mean Offspring Score", x = "GEBV of Parent Tree", size = "Number of Offspring")+
  theme(legend.position="bottom")+
  ggtitle("B")

both_plots <- ggdraw()+
  draw_plot(canopy_plot, x = 0.05, y = 0.55, height = 0.4, width = 0.9)+
  draw_plot(gebv_plot, x = 0.05, y = 0, height = 0.5, width = 0.9)

both_plots

#Save plots to file
tiff("mp_fig4.tiff", units="mm", width=180, height=180, res=300)
both_plots
dev.off()
png("mp_fig4.png", units="mm", width=180, height=180, res=300)
both_plots
dev.off()
jpeg("mp_fig4.jpeg", units="mm", width=180, height=180, res=300)
both_plots
dev.off()