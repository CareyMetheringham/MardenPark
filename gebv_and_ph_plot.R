library(ggplot2)
library(ggpubr)

#Import dataframe
gebv_and_ph <- read.csv(file = file.choose()) 

#Get adults and juveniles
gebv_and_ph_j <- gebv_and_ph[which(gebv_and_ph$Type == "Juvenile"), ]
gebv_and_ph_a <- gebv_and_ph[which(gebv_and_ph$Type == "Adult"), ]

qqnorm(gebv_and_ph_j$GEBV)
qqline(gebv_and_ph_j$GEBV)
qqnorm(gebv_and_ph_a$GEBV)
qqline(gebv_and_ph_a$GEBV)

wilcox.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV)
t.test(gebv_and_ph_a$GEBV, gebv_and_ph_j$GEBV)
mean(gebv_and_ph_j$GEBV)
mean(gebv_and_ph_a$GEBV)

juv_gebv <- ggplot(data = gebv_and_ph_j, aes(x = GEBV, y = Score_2019, group = Score_2019))+
  geom_violin(fill="skyblue")+
  theme_minimal()+
  ylab(label = "Health Score")+
  xlab("GEBV")

adult_gebv <- ggplot(data = gebv_and_ph_a, aes(x = GEBV, y = PercentScore_2019))+
  geom_point()+
  theme_minimal()+
  ylab(label = "Canopy Cover (%)")+
  xlab("GEBV")

fig2 <- ggarrange(adult_gebv, juv_gebv,
          labels = c("Adults", "Juveniles"),
          ncol = 2, nrow = 1)+
          xlab("GEBV")

fig2

tiff("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_ex_data_fig3.tiff", units="mm", width=180, height=100, res=300)
fig2
dev.off()
png("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_ex_data_fig3.png", units="mm", width=180, height=100, res=300)
fig2
dev.off()
jpeg("/Users/carey/University/Marden_Park/Figures/Final_Figures/mp_ex_data_fig3.jpeg", units="mm", width=180, height=100, res=300)
fig2
dev.off()