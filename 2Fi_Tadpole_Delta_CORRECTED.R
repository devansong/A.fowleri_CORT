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
tad <- subset(tad, dellen < 1000000000)
str(tad)

df <- tad[, c("dellen", 
              "Treatment", 
              "ClutchID",
              "delfat",
              "delmot",
              "basfat",
              "basmot",
              "SireBCI", 
              "DamBCI"
)]

df <- subset(df, dellen<10000000000)
df$Treatment <- as.numeric(as.factor(df$Treatment)) + 0
df$ClutchID <- as.numeric(as.factor(df$ClutchID)) + 0
str(df)

column_to_check <- "dellen"

# Remove rows where the specified column has NA values
df <- df[complete.cases(df[[column_to_check]]), ]


#df <- subset(df, dellen < 10)
#for(i in c(1:ncol(df))){
#  df[,i] <- as.numeric(df[,i]) #make all data points numeric
#}

X <- as.matrix(df %>% select(-dellen, -ClutchID))
Y <- df$dellen
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
  max_depth = 3,                  # Example parameter
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

write.csv(importance_matrix, file="xgboost_results_dellen.csv", row.names=FALSE)
# Print the best model performance
print("Best RMSE:")
print(model$results[which.min(model$results$RMSE), ])

# Print the final model summary
print(model$finalModel)


predictions <- model$pred$pred
actuals <- model$pred$obs

plot(actuals, predictions)


# Calculate performance metrics
#rmse <- sqrt(mean((predictions - actuals)^2))
#mae <- mean(abs(predictions - actuals))
#r_squared <- cor(predictions, actuals)^2 
#r_squared
summary(lm(predictions~actuals))

# Print performance metrics
#print(paste("RMSE:", round(rmse, 3)))
#print(paste("MAE:", round(mae, 3)))
#print(paste("R-squared:", round(r_squared, 3)))








######################Repeat for Cort#####


df <- tad[, c("cortlen", 
              "Treatment", 
              "ClutchID",
              "delfat",
              "delmot",
              "basfat",
              "basmot",
              "SireBCI", 
              "DamBCI"
)]

df <- subset(df, cortlen <1000000000)
df$Treatment <- as.numeric(as.factor(df$Treatment)) + 0
df$ClutchID <- as.numeric(as.factor(df$ClutchID)) + 0
str(df)

column_to_check <- "cortlen"

# Remove rows where the specified column has NA values
df <- df[complete.cases(df[[column_to_check]]), ]


#df <- subset(df, dellen < 10)
#for(i in c(1:ncol(df))){
#  df[,i] <- as.numeric(df[,i]) #make all data points numeric
#}

X <- as.matrix(df %>% select(-cortlen, -ClutchID))
Y <- df$cortlen
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
  max_depth = 3,                  # Example parameter
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

write.csv(importance_matrix, file="xgboost_results_cortlen.csv", row.names=FALSE)
# Print the best model performance
print("Best RMSE:")
print(model$results[which.min(model$results$RMSE), ])

# Print the final model summary
print(model$finalModel)


predictions <- model$pred$pred
actuals <- model$pred$obs

plot(actuals, predictions)
summary(lm(predictions~actuals))

# Calculate performance metrics
#rmse <- sqrt(mean((predictions - actuals)^2))
#mae <- mean(abs(predictions - actuals))
#r_squared <- cor(predictions, actuals)^2 
#r_squared
#summary(lm(predictions~actuals))

# Print performance metrics
#print(paste("RMSE:", round(rmse, 3)))
#print(paste("MAE:", round(mae, 3)))
#print(paste("R-squared:", round(r_squared, 3)))


####POWER ANalysis 


library(simr)
library(lme4)
str(df)
num_groups <- 9
df <- tad[, c("dellen", 
              "cortlen",
              "Treatment", 
              "ClutchID",
              "delfat",
              "delmot",
              "basfat",
              "basmot",
              "SireBCI", 
              "DamBCI"
)]

str(df)
df$Treatment <- as.numeric(as.factor(df$Treatment)) + 0
df$ClutchID <- as.numeric(as.factor(df$ClutchID)) + 0

num_clutches <- 9

#df <- subset(df, cortlen > 0.000000000000000001)
summary(lmer(dellen ~ basfat + basmot + Treatment + (1 | ClutchID), data = df))

model <- lmer(dellen ~ basfat + basmot + Treatment + (1 | ClutchID), data = df)
sim_treat <- powerSim(model, nsim=100, test = fcompare(y~basfat))
sim_treat

model_ext_class <- extend(model, along="ClutchID", n=18)
model_ext_class


sim_treat_class <- powerSim(model_ext_class, nsim=100, test = fcompare(y~basfat))
sim_treat_class

#NEED SAMPLE SIZE OF 18!!




# Function to perform power simulation and check results
perform_power_simulation <- function(model, num_clutches) {
  # Extend the model to simulate different numbers of clutches
  extended_model <- extend(model, along = "ClutchID", n = num_clutches)
  
  # Perform power simulation
  power_result <- powerSim(extended_model, nsim = 100)  # Adjust nsim as needed
  
  # Return the detailed result for debugging
  return(power_result)
}

power_result$x

# Define target power level
target_power <- 0.8

# Start with the current number of clutches
num_clutches <- 9

# Set a maximum number of clutches to avoid infinite loop
max_clutches <- 200

repeat {
  power_result <- perform_power_simulation(modd, num_clutches)
  
  # Print detailed power results for inspection
  print(paste("Number of clutches:", num_clutches))
  print(power_result)
  
  # Extract power estimate
  power_estimate <- power_result$x
  
  if (power_estimate >= target_power || num_clutches >= max_clutches) break
  
  # Increment the number of clutches
  num_clutches <- num_clutches + 1
}

# Print the result
print(paste("Number of clutches needed for", target_power * 100, "% power is:", power_result$x))


##############REPEAT FOR cortlen##########
model <- lmer(cortlen ~ basfat + basmot + Treatment + (1 | ClutchID), data = df)

summary(model)

sim_treat <- powerSim(model, nsim=100, test = fcompare(y~basfat))
sim_treat

model_ext_class <- extend(model, along="ClutchID", n=31)
model_ext_class


sim_treat_class <- powerSim(model_ext_class, nsim=100, test = fcompare(y~basfat))
sim_treat_class












extended_model <- extend(modd, along = "ClutchID", n = num_clutches)
power_result <- powerSim(extended_model, nsim = 100)
power_result$x

# Function to perform power simulation and check results
perform_power_simulation <- function(model, num_clutches) {
  # Extend the model to simulate different numbers of clutches
  extended_model <- extend(model, along = "ClutchID", n = num_clutches)
  
  # Perform power simulation
  power_result <- powerSim(extended_model, nsim = 100)  # Adjust nsim as needed
  
  # Return the detailed result for debugging
  return(power_result)
}

# Define target power level
target_power <- 0.8

# Start with the current number of clutches
num_clutches <- 9

# Set a maximum number of clutches to avoid infinite loop
max_clutches <- 200

repeat {
  power_result <- perform_power_simulation(modd, num_clutches)
  
  # Print detailed power results for inspection
  print(paste("Number of clutches:", num_clutches))
  print(power_result)
  
  # Extract power estimate
  power_estimate <- power_result$x
  
  if (power_estimate >= target_power || num_clutches >= max_clutches) break
  
  # Increment the number of clutches
  num_clutches <- num_clutches + 1
}

# Print the result
print(paste("Number of clutches needed for", target_power * 100, "% power is:", power_result$x))





