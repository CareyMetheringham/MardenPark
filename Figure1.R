#Create Figure 1 for the Marden Park Ash Paper
#Phenotypic estimates of ash dieback damage as a percentage of trees scored each year.
#A) Percentage canopy coverage estimates of adult trees and B) Health scores
#(where 5 is the most healthy) of juvenile trees as a percentage of the total trees
#scored for each age group in 2019 and 2021. 

library(ggplot2)
library(cowplot)
library(readxl)

#Import the data from Supplementary Data 1
pheno_data_unfiltered <- read_excel("~/Documents/Ash/Data_S1_phenotyping.xlsx")

replicates <- c("S31R", "S54R","S67R", "S76R", "S82R", "S223R", "S332Q", "S410Q", "S415R", "S571Q",
                "S655R", "S868R", "S885R", "S895R", "S901Q", "S925Q", "S928Q", "S939Q", "S954R", "S987Q",
                "SA", "SB", "SE", "SH", "SJ", "S394C")
pheno_data <- pheno_data_unfiltered[which(!(pheno_data_unfiltered$Sample %in% replicates)),]

#Pull out adults and juveniles
adult_data <- pheno_data[which(pheno_data$Type == "Adult"),]
juv_data <- pheno_data[which(pheno_data$Type == "Juvenile"),]

#Split percentage canopy cover into 5 categories
adult_cat_score_2019 <-cut(as.numeric(adult_data$PercentScore_2019),
                           breaks=c(0, 1, 25, 50, 75, 100),
                           labels=c("0","0<x<=25","25<x<=50","50<x<=75",">75"),
                           include.lowest = T)
adult_cat_score_2021 <-cut(as.numeric(adult_data$PercentScore_2021),
                           breaks=c(0, 1, 25, 50, 75, 100),
                           labels=c("0","0<x<=25","25<x<=50","50<x<=75",">75"),
                           include.lowest = T)

#Juveniles in 2019 
juv_data_2019 <- juv_data[which(juv_data$Score_2019 %in% c("1","2","3","4","5")),]
jpercent_2019 <- table(juv_data_2019$Score_2019)/sum(!(is.na(juv_data_2019$Score_2019)))*100

#Juveniles in 2021 
juv_data_2021 <- juv_data[which(juv_data$Score_2021 %in% c("1","2","3","4","5")),]
jpercent_2021 <- table(juv_data_2021$Score_2021)/sum(!(is.na(juv_data_2021$Score_2021)))*100

#Calculate percentages of each years population - adults
percent_2019 <- table(adult_cat_score_2019)/sum(!(is.na(adult_cat_score_2019)))*100
percent_2021 <- table(adult_cat_score_2021)/sum(!(is.na(adult_cat_score_2021)))*100

#Create data frame for adult scores
adult_scores <- data.frame(Year = c(rep(2019,5), rep(2021,5)),
                           Score = as.factor(rep(c("0","0<x<=25","25<x<=50","50<x<=75",">75"),2)),
                           Count = c(table(adult_cat_score_2019), table(adult_cat_score_2021)),
                           Percentage = c(percent_2019, percent_2021))

#Create data frame for juvenile scores 
juv_scores <- data.frame(Year = c(rep(2019,5), rep(2021,5)),
                         Score = as.factor(rep(c(1,2,3,4,5),2)),
                         Count = c(table(juv_data_2019$Score_2019), table(juv_data_2021$Score_2021)),
                         Percentage = c(jpercent_2019, jpercent_2021))

#Plot output as stacked barchart
juv_plot_stacked <- ggplot(data = juv_scores, aes(fill=Score, y=Count, x=Year)) +
    geom_bar(position="stack", stat="identity")+
    ylab("Trees Scored")+
    xlab("Year")+
    labs(title = "A: Juveniles")+
    theme_minimal()+
    theme(legend.position="bottom", legend.title = element_blank())+
    scale_fill_brewer(palette = "Blues")+
    scale_x_discrete(limit = c(2019, 2021))

adult_scores$Score <- factor(adult_scores$Score, levels = c("0","0<x<=25","25<x<=50","50<x<=75",">75"))
adult_plot_stacked <- ggplot(data = adult_scores, aes(fill=as.factor(Score), y=Count, x=Year)) +
    geom_bar(position="stack", stat="identity")+
    ylab("")+
    xlab("Year")+
    labs(title = "B: Adults")+
    theme_minimal()+
    theme(legend.position="bottom", legend.title = element_blank())+
    scale_fill_brewer(palette = "Oranges")+
    scale_x_discrete(limit = c(2019, 2021))

combined_plot_stacked_1 <- ggdraw() +
    draw_plot(juv_plot_stacked, x = 0.0, y = 0.05, width = 0.5, height = 0.9) +
    draw_plot(adult_plot_stacked, x = 0.45, y = 0.05, width = 0.5, height = 0.9)

combined_plot_stacked_1

#Save output to file
tiff("~/Documents/Ash/Figures/fig1_supp.tiff", units="mm", width=180, height=120, res=300)
combined_plot_stacked_1
dev.off()
png("~/Documents/Ash/Figures/fig1_supp.png", units="mm", width=180, height=120, res=300)
combined_plot_stacked_1
dev.off()
jpeg("~/Documents/Ash/Figures/fig1_supp.jpg", units="mm", width=180, height=120, res=300)
combined_plot_stacked_1
dev.off()

#Plot separate years
juv_plot_stacked_2019 <- ggplot(data = juv_scores[juv_scores$Year=="2019",], aes(fill=Score, y=Percentage, x = Year)) +
    geom_bar(position="stack", stat="identity")+
    ylab("Percentage of Trees Scored")+
    xlab("")+
    labs(title = "A: Juveniles - 2019")+
    theme_minimal()+
    theme(legend.position="bottom", legend.title = element_blank())+
    scale_fill_brewer(palette = "Blues")+
    scale_x_discrete(limit = c(2019))

adult_scores$Score <- factor(adult_scores$Score, levels = c("0","0<x<=25","25<x<=50","50<x<=75",">75"))
adult_plot_stacked_2019 <- ggplot(data = adult_scores[adult_scores$Year=="2019",], aes(fill=as.factor(Score), y=Percentage, x=Year)) +
    geom_bar(position="stack", stat="identity")+
    ylab("")+
    xlab("")+
    labs(title = "B: Adults - 2019")+
    theme_minimal()+
    theme(legend.position="bottom", legend.title = element_blank())+
    scale_fill_brewer(palette = "Oranges")+
    scale_x_discrete(limit = c(2019, 2021))

combined_plot_stacked_2 <- ggdraw() +
    draw_plot(juv_plot_stacked_2019, x = 0.0, y = 0.05, width = 0.5, height = 0.9) +
    draw_plot(adult_plot_stacked_2019, x = 0.45, y = 0.05, width = 0.5, height = 0.9)

combined_plot_stacked_2

#Save output to file
tiff("~/Documents/Ash/Figures/fig1.tiff", units="mm", width=180, height=120, res=300)
combined_plot_stacked_2
dev.off()
png("~/Documents/Ash/Figures/fig1.png", units="mm", width=180, height=120, res=300)
combined_plot_stacked_2
dev.off()
jpeg("~/Documents/Ash/Figures/fig1.jpg", units="mm", width=180, height=120, res=300)
combined_plot_stacked_2
dev.off()
