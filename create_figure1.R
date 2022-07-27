#Create Figure 1 for the Marden Park Ash Paper
#Phenotypic estimates of ash dieback damage as a percentage of trees scored each year.
#A) Percentage canopy coverage estimates of adult trees and B) Health scores
#(where 5 is the most healthy) of juvenile trees as a percentage of the total trees
#scored for each age group in 2019 and 2021. 

library(ggplot2)
library(cowplot)

#Import the data from Supplementary Data 1
pheno_data <- read.csv(file.choose())

#Pull out adults and juveniles
adult_data <- pheno_data[which(pheno_data$Type == "Adult"),]
juv_data <- pheno_data[which(pheno_data$Type == "Juvenile"),]

#Split percentage canopy cover into 5 categories
adult_cat_score_2019 <-cut(adult_data$PercentScore_2019,
                           breaks=c(0, 1, 25, 50, 75, 100),
                           labels=c("0","0<x<=25","25<x<=50","50<x<=75",">75"),
                           include.lowest = T)
adult_cat_score_2021 <-cut(adult_data$PercentScore_2021,
                           breaks=c(0, 1, 25, 50, 75, 100),
                           labels=c("0","0<x<=25","25<x<=50","50<x<=75",">75"),
                           include.lowest = T)

#Calculate percentages of each years population - adults
percent_2019 <- table(adult_cat_score_2019)/sum(!(is.na(adult_data$PercentScore_2019)))*100
percent_2021 <- table(adult_cat_score_2021)/sum(!(is.na(adult_data$PercentScore_2021)))*100
#Calculate percentages of each years population - juveniles
jpercent_2019 <- table(juv_data$Score_2019)/sum(!(is.na(juv_data$Score_2019)))*100
jpercent_2021 <- table(juv_data$Score_2021)/sum(!(is.na(juv_data$Score_2021)))*100

#Create data frame for adult scores
adult_scores <- data.frame(Year = c(rep(2019,5), rep(2021,5)),
                           Score = rep(c("0","0<x<=25","25<x<=50","50<x<=75",">75"),2),
                           Count = c(table(adult_cat_score_2019), table(adult_cat_score_2021)),
                           Percentage = c(percent_2019, percent_2021))

#Create data frame for juvenile scores 
juv_scores <- data.frame(Year = c(rep(2019,5), rep(2021,5)),
                          Score = as.factor(rep(c(1,2,3,4,5),2)),
                          Count = c(table(juv_data$Score_2019), table(juv_data$Score_2021)),
                          Percentage = c(jpercent_2019, jpercent_2021))

#Plot output as stacked barchart
juv_plot_stacked <- ggplot(data = juv_scores, aes(fill=Score, y=Percentage, x=Year)) +
    geom_bar(position="stack", stat="identity")+
  ylab("Percentage of Trees Scored")+
  xlab("Year")+
  labs(title = "Juveniles")+
  theme_minimal()+
  theme(legend.position="bottom", legend.title = element_blank())+
  scale_fill_brewer(palette = "Blues")+
  scale_x_discrete(limit = c(2019, 2021))

adult_plot_stacked <- ggplot(data = adult_scores, aes(fill=Score, y=Percentage, x=Year)) +
    geom_bar(position="stack", stat="identity")+
  ylab("Percentage of Trees Scored")+
  xlab("Year")+
  labs(title = "Adults")+
  theme_minimal()+
  theme(legend.position="bottom", legend.title = element_blank())+
  scale_fill_brewer(palette = "Oranges")+
  scale_x_discrete(limit = c(2019, 2021))

combined_plot_stacked <- ggdraw() +
  draw_plot(juv_plot_stacked, x = 0.0, y = 0.05, width = 0.4, height = 0.9) +
  draw_plot(adult_plot_stacked, x = 0.45, y = 0.05, width = 0.4, height = 0.9)

combined_plot_stacked

#Save output to file
tiff("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_fig1.tiff", units="mm", width=180, height=100, res=300)
combined_plot_stacked
dev.off()
png("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_fig1.png", units="mm", width=180, height=100, res=300)
combined_plot_stacked
dev.off()
jpeg("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_fig1.jpeg", units="mm", width=180, height=100, res=300)
combined_plot_stacked
dev.off()