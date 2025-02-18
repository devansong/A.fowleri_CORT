########################################################################
############March 2 2023, Code written by Anne Devan-Song###############
######Analysis for for ANFO CORT manuscript#############################
############Bend, OR, USA###############################################
########################################################################

rm(list=ls())
graphics.off()

setwd("~/Dropbox/ANFO_2023/Data")
library(ggplot2)
library(ggpubr)
library(broom)
library(tidyr)
library(geepack)
library(dplyr)
adu <-read.csv("adultfinal.csv")

adu <- subset(adu, Treatment !="C") #remove toads in treatment C
TrA <- subset(adu, Treatment =="A")
TrB <- subset(adu, Treatment =="B")
TrA$Treatment<- "Control"
TrB$Treatment<- "Treatment"
adu <- rbind(TrA, TrB)


########################################################################
#Plots to look at differences between treatments
########################################################################

Sex.labs <- c("Female", "Male")
names(Sex.labs) <- c("F", "M")

plotdelt <- ggplot(data = adu,aes(x = Treatment, y = (delta/SVL), fill = Treatment))+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=0.5, col="grey",
              show.legend=F) +
  scale_fill_manual(values = c("#444e86", "#dd5182"))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=0.7
                               ,show.legend = F)+
  geom_boxplot(notch = FALSE,  outlier.size = -1, color="black",lwd=0.5, alpha = 0,show.legend = F)+
  xlab("")+ 
  #ylab("Adult ΔCORT/SVL")+
  ylab(bquote("Adult Δ"[CORT]~"(ng/mL)/SVL (mm)"))+
  #ylab(bquote("Δ"[CORT]~"hi"))+
  theme(panel.background =  element_rect(fill = "white"),
        panel.border =      element_rect(fill = NA, colour="black")
       )+
  ggtitle("A")+
  facet_wrap(~Sex, labeller = labeller(Sex = Sex.labs))


plotdelt

fem <- subset(adu, Sex == "F")
mal <- subset(adu, Sex == "M")


plotF <- ggplot(data = fem,aes(x = Treatment, y = (ClutchSize/SVL), fill=Treatment))+
  #geom_point()+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=0.5, col="grey", show.legend=F) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=0.7
                               ,show.legend = F)+
  xlab("")+ 
  ylab("Clutch Size/Female SVL (mm)")+
  scale_fill_manual(values = c("#444e86", "#dd5182"))+
  theme(strip.background =  element_rect(fill = NA, colour = NA), 
        panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank())+
  ggtitle("B")

plotF

plotM <- ggplot(data = mal,aes(x = Treatment, y = Fert, fill = Treatment))+
  #geom_point()+
  scale_fill_manual(values = c("#444e86", "#dd5182"))+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=0.5, col="grey", show.legend=F) +
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=0.7
                               ,show.legend = F)+
  xlab("")+ 
  ylab("Percent Clutch Fertilized (%)")+
  theme(strip.background =  element_rect(fill = NA, colour = NA), 
        panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank())+
  ggtitle("C")

fert <- ggarrange(plotF, plotM)

setwd("~/Dropbox/ANFO_2023/Figures")
png("Adult_fert.png", units="in", width=6, height=6, res=300)
ggarrange(plotdelt, fert, ncol=1)
dev.off()

########################################################################
########################################################################

adu$delcor <- adu$delta/adu$SVL
adu$cortcor <- adu$CORT_1/adu$SVL
adu$Bred <- as.factor(adu$Bred)
adu$SVL <- as.numeric(adu$SVL)
#Question: Is bred/not bred correlated with stress and/or treatment? 

newadult <- adu[, c("Bred", "delcor", "SVL")]
newadult <- na.omit(newadult)


model <- glm(Bred ~., data = newadult, 
             family = binomial)

# Predict the probability (p) of diabetes positivity
probabilities <- predict(model, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
head(predicted.classes)


predictors <- c("delcor", "SVL")
# Bind the logit and tidying the data for plot
mydata <- newadult[, c("SVL", "delcor")]
mydata <- mydata %>%
  mutate(logit = log(probabilities/(1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)


ggplot(mydata, aes(logit, predictor.value))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() + 
  facet_wrap(~predictors, scales = "free_y")

plot(model, which = 4, id.n = 3)

model.data <- augment(model) %>% 
  mutate(index = 1:n()) 

model.data %>% top_n(3, .cooksd)

car::vif(model)

#delcor      SVL 
#1.000968 1.000968 
mylogit <- glm(Bred ~ delcor + Treatment + Sex+ SVL, data = adu, family = "binomial")
summary(mylogit)
###################################
new <- adu[, c("Bred", 
               "Treatment", 
               "delcor", 
               "Sex", 
               "SVL")]

new$Bred <- as.numeric(new$Bred) -1
new$Treatment <- as.factor(new$Treatment)
new$Treatment <- as.numeric(new$Treatment) + 0
new$Sex <- as.factor(new$Sex)
new$Sex <- as.numeric(new$Sex) + 0

new <- subset(new, delcor > -200)

gee_model <- geeglm(Bred ~ delcor + Treatment + Sex + SVL, id = 1:nrow(new), data = new, family = "binomial", corstr = "exchangeable")
summary(gee_model)
#################################
##NUTTIN 


logreg <- summary(mylogit) 
ggplot(model.data, aes(index, .std.resid)) + 
  geom_point(aes(color = Bred), alpha = .5) +
  theme_bw()

logre <- as.data.frame(logreg$coefficients)
setwd("~/Dropbox/ANFO_2023/Data")
write.csv(logre, file="LogisticRegResults.csv")

#Answer: No, nothing is significantly correlated with whether a toad bred or not
#Linear regression: Is a toad's delta related to its treatment, sex, or BCI? 
#1. Generalized Estimating Equations (GEE)
#GEEs extend generalized linear models to account for correlated observations within clusters. They don’t require explicit #group numbers if you only know that there are correlated observations. 
##########################
new2 <- adu[, c("delcor", 
               "Treatment", 
               "Sex", 
               "SVL", 
               "cortcor")]

new2$Treatment <- as.factor(new2$Treatment)
new2$Treatment <- as.numeric(new2$Treatment) + 0
new2$Sex <- as.factor(new2$Sex)
new2$Sex <- as.numeric(new2$Sex) + 0

new2 <- subset(new2, delcor > -200)

gee_model <- geeglm(log(delcor)~ Treatment + Sex + SVL + log(cortcor), id = 1:nrow(new2), data = new2, family = "gaussian", corstr = "exchangeable")
summary(gee_model)
############################


myreg <- lm(log(delcor)~ Treatment + Sex + SVL + log(cortcor), data=adu)
summ <- summary(myreg)
regre <- as.data.frame(summ$coefficients)
write.csv(regre, file="RegResults.csv")

plot(delcor~SVL, data=adu)
plot(cortcor~SVL, data=adu)

myreg$df
model_residuals = myreg$residuals
hist(model_residuals)
qqnorm(model_residuals)
# Plot the Q-Q line
qqline(model_residuals)

summary(myreg)
plot(myreg, which = 1)

library(polr)





###############NEW DO XGBOOST 
library(xgboost)
library(caret)

library(fastDummies)
library(DiagrammeR)
library(Ckmeans.1d.dp) # for xgb.ggplot.importance
library(viridis)
library(colorspace)

#tad <-read.csv("tadpolefinal.csv")
#tad <-subset(tad, Treatment !="C")
str(adu)

adu <- subset(adu, Use == "Yes")

df <- adu[, c("Sex", 
              "Treatment", 
              "BCI",
              "Bred",
              "SVL",
              "delcor", 
              "cortcor",
              "AdultID"
)]

df$Treatment <- as.numeric(as.factor(df$Treatment)) + 0
df$Sex <- as.numeric(as.factor(df$Sex)) + 0
str(df)
df$Bred <- as.numeric(df$Bred)

column_to_check <- "delcor"

# Remove rows where the specified column has NA values
df <- df[complete.cases(df[[column_to_check]]), ]


#df <- subset(df, dellen < 10)
#for(i in c(1:ncol(df))){
#  df[,i] <- as.numeric(df[,i]) #make all data points numeric
#}

X_matrix <- as.matrix(df %>% select(-delcor, -cortcor))
Y_vector <- df$delcor
#clutches <- df$ClutchID

set.seed(123) # For reproducibility
train_index <- createDataPartition(Y_vector, p = 0.7, list = FALSE)
X_train <- X_matrix[train_index, ]
Y_train <- Y_vector[train_index]
X_test <- X_matrix[-train_index, ]
Y_test <- Y_vector[-train_index]


dtrain <- xgb.DMatrix(data = X_train, label = Y_train)
dtest <- xgb.DMatrix(data = X_test, label = Y_test)

params <- list(
  objective = "reg:squarederror", # Regression task
  eval_metric = "rmse" # Root Mean Squared Error
)

model <- xgb.train(params = params, 
                   data = dtrain, 
                   nrounds = 100, # Number of boosting rounds
                   verbose = 1) # Print training log


predictions <- predict(model, X_test)

# Calculate RMSE
rmse <- sqrt(mean((predictions - Y_test)^2))
print(paste("RMSE: ", rmse))
summary(lm(Y_test~predictions))
importance_matrix <- xgb.importance(feature_names = colnames(X_train), model = model)

# Print importance matrix
print(importance_matrix)

# Plot feature importance
xgb.plot.importance(importance_matrix)



########################################################################
#REPEAT FOR Bred/NotBred but use binary 
########################################################################

X_matrix <- as.matrix(df %>% select(-Bred, -cortcor, -AdultID))
Y_vector <- df$Bred - 1 
#clutches <- df$ClutchID

set.seed(123) # For reproducibility
train_index <- createDataPartition(Y_vector, p = 0.7, list = FALSE)
X_train <- X_matrix[train_index, ]
Y_train <- Y_vector[train_index]
X_test <- X_matrix[-train_index, ]
Y_test <- Y_vector[-train_index]


dtrain <- xgb.DMatrix(data = X_train, label = Y_train)
dtest <- xgb.DMatrix(data = X_test, label = Y_test)

params <- list(
  objective = "binary:logistic", # Binary classification task
  eval_metric = "logloss" # Logarithmic loss
)

model <- xgb.train(params = params, 
                   data = dtrain, 
                   nrounds = 100, # Number of boosting rounds
                   verbose = 1) # Print training log

predictions_prob <- predict(model, X_test)
predictions <- ifelse(predictions_prob > 0.5, 1, 0)

# Confusion Matrix and Accuracy
conf_matrix <- confusionMatrix(factor(predictions), factor(Y_test))
print(conf_matrix)

accuracy <- conf_matrix$overall['Accuracy']
print(paste("Accuracy: ", accuracy))

summary(glm(Y_test~predictions, family="binomial"))

# Feature importance
importance_matrix <- xgb.importance(feature_names = colnames(X_train), model = model)
print(importance_matrix)
xgb.plot.importance(importance_matrix)





########################################################################
########################################################################
#plot relationships between various things
p1 <- ggplot(data = adu,aes(x = CORT_1, y = (delta/SVL), col=Treatment))+
  geom_point(alpha=0.4, size=4)+
  xlab("Adult Baseline CORT (ng/mL)/SVL (mm)")+ 
  #ylab("ΔCORT/SVL")+
  ylab(bquote("Adult Δ"[CORT]~"(ng/mL)/SVL (mm)"))+
  scale_color_manual(values = c("#444e86", "#dd5182"))+
  theme(strip.background =  element_rect(fill = NA, colour = NA), 
        panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill = NA, colour="black"), 
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(), 
        legend.position="top")+
  geom_smooth(method = "lm", color="black", size=0.5, alpha=0.5)+
  ggtitle("")
p1

setwd("~/Dropbox/ANFO_2023/Figures")
png("Adult_baseline_Delta.png", units="in", width=3.5, height=3.5, res=300)
p1
dev.off()


plot(adu$delcor, adu$cortcor)


################POWER ANALYSIS############################

install.packages("pwr") # Install pwr package if not already installed
library(pwr)            # Load pwr package
library(lme4)






