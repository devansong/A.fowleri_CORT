########################################################################
############March 2 2023, Code written by Anne Devan-Song###############
######Exploratory data for ANFO stress manuscript!######################
############Bend, OR, USA###############################################
########################################################################


rm(list=ls())
graphics.off()

setwd("~/Dropbox/ANFO_2023/Data")
library(ggplot2)
library(ggpubr)
library(ggcorrplot)
tad <-read.csv("tadpolefinal.csv")
tad <-subset(tad, Treatment !="C")
tad$D16_TotalLength <- as.numeric(as.character(tad$D16_TotalLength))


ggplot(data=tad, aes(y=delta/D16_TotalLength, x=))

#tad <- subset(tad, delta<0.73)
#Plot treatment vs delta
plotdelt <- ggplot(data = tad,aes(x = Treatment, y = delta/D16_TotalLength, fill = Treatment))+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA, show.legend=F) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=1.4, dodge.width = .75, color="black",alpha=0.7
                               ,show.legend = F)+
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.5, alpha = 0,show.legend = F)+
  scale_fill_manual(values = c("#003f5c", "#ffa600"))+
  #ylab("TadpoleΔCORT/Length")+ 
  ylab(bquote("Tadpole Δ"[CORT]~"(ng/mL)/Length (mm)"))+
  xlab("Parental Group")+ 
  scale_x_discrete(labels=c("A" = "Control", "B" = "Treatment"))+
  theme(strip.background =  element_rect(fill = NA, colour = NA), 
        panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank())+
  ggtitle("C")


plotdelt


tadtot <-read.csv("tadpolewithoutCORT.csv")
tadtot <-subset(tadtot, Treatment !="C")

plotlen <- ggplot(data = tadtot,aes(x = Treatment, y = Total.Length..mm., fill = Treatment))+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA, show.legend=F) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=1.4, dodge.width = .75, color="black",alpha=0.7
                               ,show.legend = F)+
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.5, alpha = 0,show.legend = F)+
  scale_fill_manual(values = c("#003f5c", "#ffa600"))+
  ylab("Tadpole Length (mm)")+ 
  xlab("Parental Group")+ 
  scale_x_discrete(labels=c("A" = "Control", "B" = "Treatment"))+
  theme(strip.background =  element_rect(fill = NA, colour = NA), 
        panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank())+
  ggtitle("A")


plotlen


plotcor <- ggplot(data = tad,aes(x = Treatment, y = CORT_1/D16_TotalLength, fill = Treatment))+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA, show.legend=F) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=1.4, dodge.width = .75, color="black",alpha=0.7
                               ,show.legend = F)+
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=0.5, alpha = 0,show.legend = F)+
  scale_fill_manual(values = c("#003f5c", "#ffa600"))+
  ylab("Tadpole Baseline CORT (ng/mL)/Length (mm)")+ 
  xlab("Parental Group")+ 
  scale_x_discrete(labels=c("A" = "Control", "B" = "Treatment"))+
  theme(strip.background =  element_rect(fill = NA, colour = NA), 
        panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank())+
  ggtitle("B")


plotcor

ggarrange(plotlen, plotcor, plotdelt)

setwd("~/Dropbox/ANFO_2023/Figures")
png("Tadpole_Len_Delt.png", units="in", width=8, height=4, res=300)
ggarrange(plotlen, plotcor, plotdelt, ncol=3)
dev.off()

tad$dellen <- (tad$delta)/(tad$D16_TotalLength)
tad$SireΔCORT <- (tad$SireDelta/tad$SireSVL)
tad$DamΔCORT <- (tad$DamDelta/tad$DamSVL)
tad$SireBaselineCORT <- (tad$SireBaselineCORT/tad$SireSVL)
tad$DamBaselineCORT <- (tad$DamBaselineCORT/tad$DamSVL)
tad$'SireFertility/SVL'<- (tad$SireFertility/tad$SireSVL)
tad$'DamFertility/SVL' <- (tad$DamFertility/tad$DamSVL)
###check for colinearity in covariates
# Compute a correlation matrix

subdf <- tad[, c("D16_TotalLength",
                 "SireBCI", 
                 "SireSVL", 
                 "SireMass",
                 "SireΔCORT", 
                  "SireFertility/SVL", 
                 "SireBaselineCORT",
                 "DamBCI",
                 "DamSVL", 
                 "DamMass",
                 "DamΔCORT",
                 "DamFertility/SVL", 
                 "DamBaselineCORT"
                 )]

for(i in c(1:ncol(subdf))){
  subdf[,i] <- as.numeric(subdf[,i]) #make all data points numeric
}

colnames(subdf) <- c("Tadpole Length", 
                     "Sire BCI", 
                     "Sire SVL", 
                     "Sire Mass",
                     "Sire ΔCORT/SVL", 
                     "Sire Fertility/SVL", 
                     "Sire Base. CORT/SVL",
                     "Dam BCI",
                     "Dam SVL", 
                     "Dam Mass",
                     "Dam ΔCORT",
                     "Dam Fertility/SVL", 
                     "Dam Base. CORT/SVL")

df <-na.omit(subdf)
corr <- round(cor(df), 1)

# Visualize
ggcorrplot(corr, p.mat = cor_pmat(df),
           hc.order = TRUE, type = "lower",
           color = c("#FC4E07", "white", "#00AFBB"),
           outline.col = "white", lab = TRUE)

setwd("~/Dropbox/ANFO_2023/Figures")
png("Cov_corr.png", units="in", width=5.7, height=5.7, res=300)
ggcorrplot(corr, p.mat = cor_pmat(df),
           hc.order = TRUE, type = "lower",
           color = c("#FC4E07", "white", "#00AFBB"),
           outline.col = "white", lab = TRUE)
dev.off()



