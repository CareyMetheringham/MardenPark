#load packages

#Import effect sizes - Supplementary Data 2
es_file <- read.csv("effect_sizes.csv")
es <- with(es_file, EES.MIA - EES.MAA)
  
#Import csv files with data on the GEBV sites and unlinked sites
#- Supplementary Data 3
gebv_sites <- read.csv(file = "gebv_sites.csv")
#Extract the genotype matrices
GTAgebv <- gebv_sites[gebv_sites$Age=="Adult", 4:ncol(gebv_sites)]
GTJgebv <- gebv_sites[gebv_sites$Age=="Juv", 4:ncol(gebv_sites)]
  
#Calculate GEBV scores 
new_gebv <- (as.matrix(gebv_sites[,-(1:3)])%*%es) + 1

# Load in phenotype data
gebv_and_ph <- read.csv("phenotypes.csv")
gebv_and_ph_j <- gebv_and_ph[which(gebv_and_ph$Type == "Juvenile"), ]
gebv_and_ph_a <- gebv_and_ph[which(gebv_and_ph$Type == "Adult"), ]
# 
qnorm(gebv_and_ph_j$GEBV)
qqline(gebv_and_ph_j$GEBV)
qqnorm(gebv_and_ph_a$GEBV)
qqline(gebv_and_ph_a$GEBV)
# 
qqnorm(gebv_and_ph_j$GEBV[which(gebv_and_ph_j$Height > 0.5)])
qqline(gebv_and_ph_j$GEBV[which(gebv_and_ph_j$Height > 0.5)])
# ```
# 
 ```{r}
wilcox.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV)
t.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV)
mean(gebv_and_ph_j$GEBV)
mean(gebv_and_ph_a$GEBV)
 ```

 ```{r}
 t.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV[which(gebv_and_ph_j$Height < 0.5)])
 wilcox.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV[which(gebv_and_ph_j$Height < 0.5)])
 t.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV[which(gebv_and_ph_j$Height >= 0.5)])
 wilcox.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV[which(gebv_and_ph_j$Height >= 0.5)])
 t.test(gebv_and_ph_j$GEBV[which(gebv_and_ph_j$Height >= 0.5)], gebv_and_ph_j$GEBV[which(gebv_and_ph_j$Height < 0.5)])
 ```
# 
# 
# Shape data for boxplots
# ```{r}
# h_cutoff <- 0.5
# 
# gebv_and_ph_j1 <-
#   gebv_and_ph_j[which(gebv_and_ph_j$Height..m. < h_cutoff), ]
# gebv_and_ph_j2 <-
#   gebv_and_ph_j[which(gebv_and_ph_j$Height..m. >= h_cutoff), ]
# gebv_for_plots <-
#   data.frame(
#     c(
#       rep("Adults (n = 133)", nrow(gebv_and_ph_a)),
#       rep("Juveniles (n = 442)", nrow(gebv_and_ph_j))
#       #rep("Juveniles < 0.5m (n = 136)", nrow(gebv_and_ph_j1)),
#       #rep("Juveniles >= 0.5m (n = 291)", nrow(gebv_and_ph_j2))
#     ),
#     c(
#       gebv_and_ph_a$GEBV,
#       gebv_and_ph_j$GEBV
#       #gebv_and_ph_j1$GEBV,
#       #gebv_and_ph_j2$GEBV
#     )
#   )
# colnames(gebv_for_plots) <- c("Cohort", "GEBV")
# ```
# ```{r}
# gebv_boxplot <- ggplot(gebv_for_plots, aes(x = Cohort, y = GEBV)) +
#   geom_boxplot(fill = c("#D55E00", "skyblue")) +
#   theme_minimal() + xlab("") +
#   theme(axis.text=element_text(size=12))
# 
# gebv_boxplot
# ```
# 
# Plot Score by location
# 
# ```{r}
# library(ggpubr)
# 
# adult_loc_2019 <- ggplot(data = gebv_and_ph_a, aes(x = Long, y = Lat, size = DBH, colour = PercentScore_2019))+
#   geom_point(show.legend = FALSE) +
#   theme_minimal() +
#   theme(legend.position="bottom", legend.box = "vertical") +
#   scale_color_gradient(low = "red", high = "orange") +
#   labs(colour = "% Canopy Cover", size = "DBH (cm)") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   ylim(51.27117, 51.27011) +
#   xlim(-0.04020264, -0.038833) +
#   guides(colour = guide_colourbar(order = 1),
#          size = guide_legend(order = 2))
# 
# juv_loc_2019 <- ggplot(data = gebv_and_ph_j, aes(x = Long, y = Lat, size = Height, colour = Score_2019))+
#   geom_point(show.legend = FALSE) +
#   theme_minimal() +
#   theme(legend.position="bottom", legend.box = "vertical") +
#   labs(colour = "Health Score", size = "Height (m)") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   ylim(51.27117, 51.27011) +
#   xlim(-0.04020264, -0.038833) +
#   guides(colour = guide_colourbar(order = 1),
#          size = guide_legend(order = 2))
# 
# adult_loc_2021 <- ggplot(data = gebv_and_ph_a, aes(x = Long, y = Lat, size = DBH, colour = PercentScore_2021))+
#   geom_point(xlab="GEBV") +
#   theme_minimal() +
#   theme(legend.position="bottom", legend.box = "vertical") +
#   scale_color_gradient(low = "red", high = "orange") +
#   labs(colour = "% Canopy Cover", size = "DBH (cm)") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   ylim(51.27117, 51.27011) +
#   xlim(-0.04020264, -0.038833) +
#   guides(colour = guide_colourbar(order = 1),
#          size = guide_legend(order = 2))
# 
# juv_loc_2021 <- ggplot(data = gebv_and_ph_j, aes(x = Long, y = Lat, size = Height, colour = Score_2021))+
#   geom_point() +
#   theme_minimal() +
#   theme(legend.position="bottom", legend.box = "vertical") +
#   labs(colour = "Health Score", size = "Height (m)") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   ylim(51.27117, 51.27011) +
#   xlim(-0.04020264, -0.038833) +
#   guides(colour = guide_colourbar(order = 1),
#          size = guide_legend(order = 2))
# 
# ```
# ```{r}
# ggarrange(adult_loc_2019, juv_loc_2019, adult_loc_2021, juv_loc_2021,
#           labels = c("Adults : 2019", "Juveniles : 2019", "Adults : 2021", "Juveniles : 2021"),
#           ncol = 2, nrow = 2,
#           heights = c(1.5, 2.2))
# ```
# ```{r}
adult_loc <- ggplot(data = gebv_and_ph_a, aes(x = Long, y = Lat, size = DBH, colour = GEBV))+
   geom_point(xlab="GEBV") +
   theme_minimal() +
   theme(legend.position="bottom", legend.box = "vertical") +
   scale_color_gradient(low = "#fee6ce", high = "#e63f0d") +
   labs(colour = "GEBV", size = "DBH (cm)") +
   ylim(51.27117, 51.27011) +
   xlim(-0.04020264, -0.038833) +
   guides(colour = guide_colourbar(order = 1),
          size = guide_legend(order = 2))
 
 juv_loc <- ggplot(data = gebv_and_ph_j, aes(x = Long, y = Lat, size = Height, colour = GEBV))+
   geom_point() +
   theme_minimal() +
   theme(legend.position="bottom", legend.box = "vertical") +
   labs(colour = "GEBV", size = "Height (m)") +
   ylim(51.27117, 51.27011) +
   xlim(-0.04020264, -0.038833) +
   guides(colour = guide_colourbar(order = 1),
          size = guide_legend(order = 2))
 
 ggarrange(adult_loc, juv_loc, 
           labels = c("Adults", "Juveniles"),
           ncol = 2, nrow = 1)
# ```
# 
# Exclude the children of tree 51
# ```{r}
children_of_51 <- related$id[which(related$dam == "S51R" | related$sire == "S51R")]
not_51_a <- gebv_and_ph_a[!which(gebv_and_ph_a$Label %in% children_of_51 & gebv_and_ph_a$Label != "S51R"),]
not_51_j <- gebv_and_ph_j[!which(gebv_and_ph_j$Label %in% children_of_51),]
wilcox.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV)
t.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV)
# ```
# 
# ```{r}
 juv_gebv <- ggplot(data = gebv_and_ph_j, aes(x = GEBV-1, y = Score_2019, group = Score_2019))+
   geom_violin(fill="skyblue")+
   theme_minimal()+
   ylab(label = "Health Score")+
   xlab("GEBV")
 
 adult_gebv <- ggplot(data = gebv_and_ph_a, aes(x = GEBV-1, y = PercentScore_2019))+
   geom_point()+
   theme_minimal()+
   ylab(label = "Canopy Cover (%)")+
   xlab("GEBV")
 
 fig2 <- ggarrange(adult_gebv, juv_gebv, 
           labels = c("Adults", "Juveniles"),
           ncol = 2, nrow = 1)+
          xlab("GEBV")
 
 fig2
 ```
 ```{r}
 tiff("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_ex_data_fig3.tiff", units="mm", width=180, height=100, res=300)
 fig2
 dev.off()
 png("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_ex_data_fig3.png", units="mm", width=180, height=100, res=300)
 fig2
 dev.off()
 jpeg("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_ex_data_fig3.jpeg", units="mm", width=180, height=100, res=300)
 fig2
 dev.off()
 ```
