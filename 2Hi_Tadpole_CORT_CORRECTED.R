########################################################################
############March 2 2023, Code written by Anne Devan-Song###############
######XGBOOST for tadpole data##########################################
############Bend, OR, USA###############################################
########################################################################



rm(list=ls())
graphics.off()

setwd("~/Dropbox/ANFO_2023/Data")
library(xgboost)
library(caret)

library(fastDummies)
library(DiagrammeR)
library(Ckmeans.1d.dp) # for xgb.ggplot.importance
library(viridis)
library(colorspace)

tad <-read.csv("tadpolefinal.csv")
tad <-subset(tad, Treatment !="C")
tad$D16_TotalLength <- as.numeric(as.character(tad$D16_TotalLength))

tad$dellen <- (tad$delta)/(tad$D16_TotalLength)
tad$cortlen <- (tad$CORT_1)/(tad$D16_TotalLength)
tad$delfat <- (tad$SireDelta/tad$SireSVL)
tad$delmot <- (tad$DamDelta/tad$DamSVL)
tad$basfat <- (tad$SireBaselineCORT/tad$SireSVL)
tad$basmot <- (tad$DamBaselineCORT/tad$DamSVL)
tad$ferfat <- (tad$SireFertility/tad$SireSVL)
tad$fermot <- (tad$DamFertility/tad$DamSVL)

#tad <- subset(tad, Zdelt<3 & Zdelt>-3)
#Consider adding parent baseline CORT

df <- tad[, c("CORT_1", 
              "D16_TotalLength",
              "Treatment", 
              "ClutchID",
              "delfat",
              "delmot",
              "basfat",
              "basmot",
              "ferfat",
              "fermot",
              "SireBCI", 
              "SireSVL", 
              "SireMass",
              "DamBCI", 
              "DamSVL", 
              "DamMass"
)]


df$Treatment <- as.numeric(as.factor(df$Treatment)) + 0
df$ClutchID <- as.numeric(as.factor(df$ClutchID)) + 0

df <- subset(df, CORT_1 < 10)


for(i in c(1:ncol(df))){
  df[,i] <- as.numeric(df[,i]) #make all data points numeric
}

X <- as.matrix(df %>% select(-CORT_1, -ClutchID))
Y <- df$CORT_1
clutches <- df$ClutchID

cv_control <- trainControl(
  method = "cv",
  number = 9,  # Number of folds
  index = createFolds(clutches, k = 9),  # Create folds based on clutches
  savePredictions = "final",  # Save predictions for later use
  summaryFunction = defaultSummary  # Use default summary function
)

xgb_params <- list(
  objective = "reg:squarederror", # Regression task
  eval_metric = "rmse",           # Root Mean Squared Error
  max_depth = 6,                  # Example parameter
  eta = 0.3,                      # Learning rate
  nrounds = 100                   # Number of boosting rounds
)

model <- train(
  x = X,
  y = Y,
  method = "xgbTree",
  trControl = cv_control,
  tuneGrid = expand.grid(
    nrounds = xgb_params$nrounds,
    max_depth = xgb_params$max_depth,
    eta = xgb_params$eta,
    gamma = 0,  # Default value
    colsample_bytree = 1,  # Default value
    min_child_weight = 1,  # Default value
    subsample = 1  # Default value
  ),
  metric = "RMSE"  # Metric for evaluation
)
print(model)

final_model <- model$finalModel

importance_matrix <- xgb.importance(feature_names = colnames(X), model = final_model)
print(importance_matrix)
xgb.plot.importance(importance_matrix)

# Print the best model performance
print("Best RMSE:")
print(model$results[which.min(model$results$RMSE), ])

# Print the final model summary
print(model$finalModel)


predictions <- model$pred$pred
actuals <- model$pred$obs

summary(lm(predictions~actuals))
# Calculate performance metrics
rmse <- sqrt(mean((predictions - actuals)^2))
mae <- mean(abs(predictions - actuals))
r_squared <- cor(predictions, actuals)^2 
r_squared
summary(lm(predictions~actuals))

# Print performance metrics
print(paste("RMSE:", round(rmse, 3)))
print(paste("MAE:", round(mae, 3)))
print(paste("R-squared:", round(r_squared, 3)))

summary(lm(df$cortlen~df$basfat + df$basmot + df$ClutchID)) + 
             #df$delmot + df$delfat + 
             df$ClutchID)


summary(lm(tad$D16_TotalLength~tad$basfat + tad$basmot + tad$delmot + tad$delfat + tad$SireBCI + tad$ClutchID))

library(lme4)
library(lmerTest) # Adds p-values and significance tests
library(dplyr)


model <- lmer(
  CORT_1 ~ D16_TotalLength + SireBCI + DamBCI + basfat +basmot+ (1 | ClutchID),
  data = tad
)

model <- lmer(
  cortlen ~ SireBCI + DamBCI + basfat + basmot+ (1 | ClutchID),
  data = tad
)

summary(model)

summary(lm(df$cortlen~ df$SireBCI + df$DamBCI + df$basfat + df$basmot + df$ClutchID))

# Print the model summary
summary(model)




