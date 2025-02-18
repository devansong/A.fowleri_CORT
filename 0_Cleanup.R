########################################################################
############March 2 2023, Code written by Anne Devan-Song###############
######Analysis for for ANFO CORT manuscript#############################
############Bend, OR, USA###############################################
########################################################################

#input: AdultsALL.csv, Male.csv, Female.csv, Tadpoles.csv, Tadpoles_complete_data.csv
#output: adultfinal.csv, tadpolefinal.csv, tadpolecwithourtCORT.csv 

rm(list=ls()) #clear working space
graphics.off() #clear plots 

setwd("~/Dropbox/ANFO_2023/Data") #set working directory to where data are stored 
#load relevant libraries
library(ggplot2) 
library(ggpubr)
library(data.table)
library(zoo)
library(dplyr)

adu <-read.csv("AdultsALL.csv") #load adult data
str(adu) #see the structure of this dataset

#Convert FertRate to another column Fert, which is in proportions, and numerical format
adu$Fert <- gsub('%','',adu$FertRate) #remove the % sign
adu$Fert <- as.numeric(adu$Fert) 
adu$Fert <- adu$Fert/100

#Make sure Sex has only 2 factors
adu$Sex <- gsub(' ','',adu$Sex) #removes all spaces
adu$Sex <- as.factor(adu$Sex)

#Calculate adult delta (change in CORT
adu$delta <- adu$CORT_2-adu$CORT_1
adu <- subset(adu, select=-c(CORT_3))


###############################################################################################
###############################################################################################
#Convert adult data to long form, so each CORT measurement has its own row
#long should have 106 x 3 (318) rows
long <- melt(setDT(adu), id.vars = c("AdultID",
                                     "Sex", 
                                     "ClutchID",
                                     "Treatment",
                                     "BCI",
                                     "ClutchSize",
                                     "FertRate",
                                     "X", 
                                     "Data_Complete",
                                     "Use",
                                     "Notes",
                                     "delta",
                                     "Fert"), 
             variable.name = "CORT")
#Calculate AUC for each toad
AUC <- function(x, fs)
  setNames(as.data.frame(
    lapply(fs, function(f) sum(diff(x$value)*rollmean(x[,f], 2)))), 
    fs)
aduAUC <- long %>% 
  group_by(AdultID) %>%
  arrange(value) %>%
  do(AUC(., grep("value", names(.), value=T)))
colnames(aduAUC) <- c("AdultID", "AUC")

#merge newly created AUC dataframe with adult dataframe
#new should now have 106 observations 
new <- merge(adu, aduAUC, by="AdultID", all=FALSE)
###############################################################################################
###############################################################################################


new<- adu


#Create a 'bred/not bred' column 
bred <- subset(new, ClutchID !="NA") #subset all toads with clutch ID (i.e., bred)
unbred <- new[is.na(new$ClutchID),] #subset all toads with no clutch ID
bred$Bred <- "Yes"
unbred$Bred <- "No"
withbred <- rbind(bred, unbred) 

#add SVL and Mass 
fem <- read.csv("Female.csv")
mal <- read.csv("Male.csv")

fem <- fem[, c("ID", "SVL", "Mass_Initial")]
colnames(fem) <- c("AdultID", "SVL", "Mass")
mal <- mal[, c("ID", "SVL", "Mass_Initial")]
colnames(mal) <- c("AdultID", "SVL", "Mass")
svlmass <- rbind(mal, fem)
mergsize <- merge(withbred, svlmass, by="AdultID", all=TRUE)
#mergsize$BC <- (mergsize$SVL)/(mergsize$Mass)

#create new csv in data folder, that will now be used for analysis
write.csv(mergsize, file="adultfinal.csv", row.names=FALSE)


########################################################################
########################################################################
########################################################################
########################################################################
#Now, prep Tadpole data 

tad <- read.csv("Tadpoles.csv") #read in tadpole data


#calculate delta for tadpoles
tad$delta <- tad$CORT_3-tad$CORT_1
tad <- subset(tad, select=-c(CORT_2))


########################################################################
########################################################################
#Calculate AUC for tadpoles, same code as adults
tadlong <- melt(setDT(tad), id.vars = c("TadpoleID", "ClutchID", "Treatment",
                                        "D16_TotalLength","D16_TailLength","D16_BodyWidth",
                                        "Data_Compete", "delta"), variable.name = "CORT")
tadlongform <- tadlong %>% 
  group_by(TadpoleID) %>%
  arrange(value) %>%
  do(AUC(., grep("value", names(.), value=T)))
colnames(tadlongform) <- c("TadpoleID", "AUC")
newtad <- merge(tad, tadlongform, by="TadpoleID", all=FALSE)

#extract enclosure ID from tadpole ID 
########################################################################
########################################################################



newtad$enclosure <- sub(".*-", "", newtad$TadpoleID)  
newtad$enclosure <-substr(newtad$enclosure,1,1)
########################################################################
########################################################################
####Match tadpoles to their parents data
adult <- read.csv("adultfinal.csv")

#Create a father only dataframe
father <- subset(adult, Sex =="M" & ClutchID!="")
father <- father[, c(
  "AdultID","ClutchID","BCI","Fert", "delta", "SVL", "Mass", "CORT_1"
)]
colnames(father) <- c("SireID", "ClutchID", "SireBCI", "SireFertility", "SireDelta",
                      "SireSVL", "SireMass", "SireBaselineCORT")

#Create a mother only dataframe
mother <- subset(adult, Sex =="F" & ClutchID!="")
mother <- mother[, c(
  "AdultID","ClutchID","BCI","ClutchSize", "delta","SVL", "Mass", "CORT_1"  
)]
colnames(mother) <- c("DamID", "ClutchID", "DamBCI", "DamFertility", "DamDelta", 
                      "DamSVL", "DamMass", "DamBaselineCORT")

#Merge tadpole data with father and mother

newtad <- tad
tadfath <- merge(newtad, father, by="ClutchID", all=FALSE)
tadpar <-merge(tadfath, mother, by="ClutchID", all=FALSE)


##################################################################################
#For tadpole delta, CORT_1 and AUC, calculate Z-values to identify outliers 
delSD <- sd(tadpar$delta, na.rm = TRUE)
delSE <- delSD/(sqrt(107))
delmean <- mean(tadpar$delta, na.rm=TRUE) #0.1882056
tadpar$Zdelt <- (tadpar$delta-delmean)/(delSD)

cortSD <- sd(tadpar$CORT_1, na.rm = TRUE)
cortSE <- cortSD/(sqrt(107))
cortmean <- mean(tadpar$CORT_1, na.rm=TRUE) #0.1882056
tadpar$Zcort <- (tadpar$CORT_1-cortmean)/(cortSD)

tadpar$D16_TotalLength <- as.numeric(tadpar$D16_TotalLength)

tadpar$dellen <- (tadpar$delta)/(tadpar$D16_TotalLength)
tadpar$cortlen <- (tadpar$CORT_1)/(tadpar$D16_TotalLength)
tadpar$delfat <- (tadpar$SireDelta/tadpar$SireSVL)
tadpar$delmot <- (tadpar$DamDelta/tadpar$DamSVL)
tadpar$basfat <- (tadpar$SireBaselineCORT/tadpar$SireSVL)
tadpar$basmot <- (tadpar$DamBaselineCORT/tadpar$DamSVL)
tadpar$ferfat <- (tadpar$SireFertility/tadpar$SireSVL)
tadpar$fermot <- (tadpar$DamFertility/tadpar$DamSVL)

str(tadpar)
#create final tadpole csv, to be used in future analyses 
write.csv(tadpar, file="tadpolefinal.csv", row.names=FALSE)


#Merge parent with complete data for tadpoles 

tadcom <- read.csv("Tadpoles_complete_data.csv")
tadcom <- tadcom[c(1:281),] #remove empty rows

colnames(tadcom)[colnames(tadcom) == "Clutch.ID"] ="ClutchID"
tadcomfa <- merge(tadcom, father, by ="ClutchID", all=FALSE)
tadcompar <- merge(tadcomfa, mother, by ="ClutchID", all=FALSE)


tad <-read.csv("tadpolewithoutCORT.csv")
tad <-subset(tad, Treatment !="C")

#tad$dellen <- (tad$delta)/(tad$D16_TotalLength)
tadcompar $delfat <- (tadcompar$SireDelta/tadcompar$SireSVL)
tadcompar $delmot <- (tadcompar$DamDelta/tadcompar$DamSVL)
tadcompar$basfat <- (tadcompar$SireBaselineCORT/tadcompar$SireSVL)
tadcompar$basmot <- (tadcompar$DamBaselineCORT/tadcompar$DamSVL)
tadcompar$ferfat <- (tadcompar$SireFertility/tadcompar$SireSVL)
tadcompar$fermot <- (tadcompar$DamFertility/tadcompar$DamSVL)


write.csv(tadcompar, file = "tadpolewithoutCORT.csv", row.names=FALSE)



