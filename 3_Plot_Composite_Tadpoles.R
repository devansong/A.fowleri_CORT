########################################################################
############March 24 2023, Code written by Anne Devan-Song##############
######Plot composite XGBOOst results WITH trends #############################
###ANFO paper ########################
############Bend, OR, USA###############################################
########################################################################

rm(list=ls())
graphics.off()
library(tidyverse)
setwd("~/Dropbox/ANFO_2023/Data")
library(ggplot2) 
library(ggpubr)
library(data.table)
library(zoo)
library(dplyr)
library(viridis)
library(colorspace)

tad <-read.csv("tadpolefinal.csv")
tad <-subset(tad, Treatment !="C")
tad$D16_TotalLength <- as.numeric(as.character(tad$D16_TotalLength))

#Scale BCI so it's between 0 and 1
tad$DamBCI_Scale <- ((tad$DamBCI - min(tad$DamBCI)) / (max(tad$DamBCI) - min(tad$DamBCI)))+1  # Scale to 0/1
hist(tad$DamBCI_Scale) 
tad$SireBCI_Scale <- ((tad$SireBCI - min(tad$SireBCI)) / (max(tad$SireBCI) - min(tad$SireBCI)))+1  # Scale to 0/1
hist(tad$SireBCI_Scale) 

tad$dellen <- (tad$delta)/(tad$D16_TotalLength)
tad$auclen <- (tad$AUC)/(tad$D16_TotalLength)
tad$cortlen <- (tad$CORT_1)/(tad$D16_TotalLength)

tad$delfat <- (tad$SireDelta/tad$SireSVL)
tad$delmot <- (tad$DamDelta/tad$DamSVL)

tad$aucfat <- (tad$SireAUC/tad$SireSVL)
tad$aucmot <- (tad$DamAUC/tad$DamSVL)
tad$basfat <- (tad$SireBaselineCORT/tad$SireSVL)
tad$basmot <- (tad$DamBaselineCORT/tad$DamSVL)
tad$ferfat <- (tad$SireFert/tad$SireSVL)
tad$fermot <- (tad$DamFert/tad$DamSVL)
tad$FatherBC <- tad$SireSVL/tad$SireMass
tad$MotherBC <- tad$DamSVL/tad$DamMass




pD <- ggplot(data = tad,aes(x = basfat, y = cortlen))+
  geom_point(alpha=0.4, size=4, color="#003f5c")+
  xlab("Sire Baseline CORT (ng/mL)/SVL (mm)")+ 
  labs(y = paste0("Tadpole Baseline", 
                  "\n", 
                  "CORT (ng/mL)/Length (mm)"))+ 
  #ylab("Tadpole Baseline CORT (ng/mL)/Length (mm)")+
  #scale_color_manual(values = c("#444e86"))+
  theme(strip.background =  element_rect(fill = NA, colour = NA), 
        panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(), 
        legend.position="top")+
  geom_smooth(method = "lm", color="black", size=0.5, alpha=0.5)+
  ggtitle("B")
pD

pF <- ggplot(data = tad,aes(x = basfat, y = dellen))+
  geom_point(alpha=0.4, size=4, color="#003f5c")+
  xlab("Sire Baseline CORT (ng/mL)/SVL (mm)")+ 
  ylab("Tadpole ΔCORT (ng/mL)/Length (mm)")+
  ylab(bquote("Tadpole ΔCORT (ng/mL)/Length (mm)"))+
  #scale_color_manual(values = c("#444e86"))+
  theme(strip.background =  element_rect(fill = NA, colour = NA), 
        panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(), 
        legend.position="top")+
  geom_smooth(method = "lm", color="black", size=0.5, alpha=0.5)+
  ggtitle("D")
pF


cortcorr <- read.csv("xgboost_results_cortlen.csv")
cortcorr$group <- "Size-corrected CORT (R2 = 2%)"
cortcorr <- cortcorr[order(cortcorr$Importance, decreasing = TRUE), ]
cortcorr$Feat <- c("Sire Base. CORT", 
                   "Dam ΔCORT",
                   "Dam Base. CORT", 
              "Sire ΔCORT", 
              "Sire BCI", 
              "Treatment")

delcorr <- read.csv("xgboost_results_dellen.csv")
delcorr$group <- "Size-corrected Delta (R2 = 3%)"
delcorr <- delcorr[order(delcorr$Importance, decreasing = TRUE), ]
delcorr$Feat <- c("Sire Base. CORT", 
                  "Dam ΔCORT", 
                  "Sire ΔCORT", 
                  "Treatment",
                  "Dam Base. CORT")


df <- rbind (cortcorr, delcorr)
df <- df[order(df$Importance, decreasing = TRUE), ]


#df$Feature <- c("Father CORT/SVL", "Father CORT/SVL", "Father Mass",
#                "Father BCI",
#                "Mother SVL",
#                "Mother Delta/SVL",
#                "Mother Delta",
#                "Mother CORT",
#                "Mother Delta",
#                "Father SVL",
#                "Mother Mass",
#                "Mother SVL")

cortcorr <- subset(df, group == "Size-corrected CORT (R2 = 2%)")
delcorr <- subset(df, group == "Size-corrected Delta (R2 = 3%)")




#plot2 <- ggarrange(pB, pD, pF, ncol=1, nrow=3)
#plot2

#ggarrange(plot, plot2)


pC <- ggplot(cortcorr, aes(reorder(Feat, Importance), 
                      Importance, fill=factor(Importance))) + 
  geom_bar(stat='identity', position=position_dodge(), color="black", size=0.2,show.legend=FALSE) + 
  scale_fill_discrete_sequential(palette = "Viridis")+
  #geom_errorbar(aes(ymin=Importance-SE, ymax=Importance+SE), width=0.7,
  #              position=position_dodge(.9))+
  coord_flip()+
  theme_classic()+
  xlab("")+
  labs(y = paste0("Relative Importance in Predicting", 
                  "\n", 
                  "Tadpole Baseline CORT/Length"))+   
  ggtitle("A")+
  annotate(geom="text", x=2.5, y=0.5, label="2.3%",
                         color="maroon")
pC

pE <- ggplot(delcorr, aes(reorder(Feat, Importance), 
                      Importance, fill=factor(Importance))) + 
  geom_bar(stat='identity', position=position_dodge(), color="black", size=0.2,show.legend=FALSE) + 
  scale_fill_discrete_sequential(palette = "Viridis")+
  #geom_errorbar(aes(ymin=Importance-SE, ymax=Importance+SE), width=0.7,
  #              position=position_dodge(.9))+
  coord_flip()+
  theme_classic()+
  xlab("")+
  labs(y = paste0("Relative Importance in Predicting", 
                  "\n", 
                  "Tadpole ΔCORT/Length"))+  ggtitle("C")+
  annotate(geom="text", x=2, y=0.5, label="3.7%",
           color="maroon")
pE

ggarrange(pC, pD, pE, pF, ncol = 2, nrow = 3)


setwd("~/Dropbox/ANFO_2023/Figures")
png("Tadpole_Composite_XGBOOST_AND_scatterplots.png", units="in", width=7.5, height=6, res=300) #this sets dimensions of image
ggarrange(pC, pD, pE, pF, ncol = 2, nrow = 2)
dev.off()





####Extra
#plot <- ggplot(df, aes(reorder(Feature, Average.Importance), 
#                       Average.Importance, fill=factor(Average.Importance))) + 
  geom_bar(stat='identity', position=position_dodge(), color="black", size=0.2,show.legend=FALSE) + 
  scale_fill_discrete_sequential(palette = "Viridis")+
  geom_errorbar(aes(ymin=Average.Importance-SE, ymax=Average.Importance+SE), width=0.7,
                position=position_dodge(.9))+
  coord_flip()+
  theme_classic()+
  labs(y = "Average Importance", x ="Feature")+
  ggtitle("")+
  facet_wrap(~group, scale="free", ncol = 1)
plot

#ng/ml
