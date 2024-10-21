###################################################################
####################### Multiple Regression ####################### 
###################################################################

###################################################################
# Import Packages
###################################################################

library(readxl)
library(car)
library(dplyr)
library(fitdistrplus)
library(corrplot)
library(fastDummies)
library(reshape2)
library(quantreg)
library(ggplot2)

###################################################################
# Import Data set - The World Values Survey Wave 7
###################################################################

WVS = read.csv("WVS_clean.csv")

###################################################################
# Omit Na values
###################################################################
#WVS = read.csv("WVS_Cross-National_Wave_7_csv_v6_.csv")
#sum(is.na(WVS))
#WVS = na.omit(WVS)

###############################################################################
# Filter 10 most + least religious countries by average DISBELIEF index
###############################################################################

sum(which(is.na(WVS)))
average_religiosity <- WVS %>%
  group_by(B_COUNTRY_ALPHA) %>%
  summarise(mean_religiosity = mean(DISBELIEF, na.rm = TRUE)) %>%
  arrange(desc(mean_religiosity)) %>%
  top_n(11, mean_religiosity)
most_secular_countries <- average_religiosity$B_COUNTRY_ALPHA
most_secular_countries <- setdiff(most_secular_countries, "AND")
WVS_secular <- WVS %>% filter(B_COUNTRY_ALPHA %in% most_secular_countries)

average_religiosity <- WVS %>%
  group_by(B_COUNTRY_ALPHA) %>%
  summarise(mean_religiosity = mean(DISBELIEF, na.rm = TRUE)) %>%
  arrange(mean_religiosity) %>%  # Sort in ascending order to find the most secular
  top_n(-10, mean_religiosity)  # Negative to get the lowest values
most_religeous_countries <- average_religiosity$B_COUNTRY_ALPHA
WVS_relig <- WVS %>% filter(B_COUNTRY_ALPHA %in% most_religeous_countries)

WVS_20 <- bind_rows(
  mutate(WVS_secular, ReligCountry = "Secular"),
  mutate(WVS_relig, ReligCountry = "Religious")
)

###############################################################################
# Create a Dummy Variable - ReligCountry 
### indicates if the individual record is from religious country
###############################################################################

WVS_20$ReligCountry = factor(WVS_20$ReligCountry)
head(WVS_20)

###############################################################################
# Remove values that used to calculate RESEMAVAL value
###############################################################################

WVS_20 = WVS_20[!names(WVS_20) %in% c("VOICE","CHOICE","AUTONOMY","EQUALITY")]

###############################################################################
# RESEMAVAL - normality check
### Plot histogram
### Create qqplot
### Computes descriptive parameters of the empirical distribution
###############################################################################
par(mar = c(4, 4, 2, 2))
hist(WVS_20$RESEMAVAL)
qqnorm(WVS_20$RESEMAVAL)
qqline(WVS_20$RESEMAVAL)
d = descdist(WVS_20$RESEMAVAL, discrete = F)

###############################################################################
# Covariance & Correlation 
###############################################################################

numeric_cols = sapply(WVS_20,is.numeric)
cor_WVS = cor(WVS_20[,unlist(numeric_cols)] ,method = "pearson")
corrplot(cor_WVS, method = "number", type = "upper" )

###############################################################################
# Removing High correlated independent variables 
###############################################################################

WVS_20 = WVS_20[!names(WVS_20) %in% c("SACSECVAL", "G_TOWNSIZE", "H_URBRURAL")]

###############################################################################
# Covariance & Correlation (Highly correlated independent variables removed)
###############################################################################

numeric_cols = sapply(WVS_20,is.numeric)
cor_WVS = cor(WVS_20[,unlist(numeric_cols)] ,method = "pearson")
corrplot(cor_WVS, method = "number", type = "upper" )

###############################################################################
# Function Build - 
### iteratively builds a linear regression model 
### by selecting significant terms based on their p-values 
### while considering interaction effects with the variable ReligCountry
###############################################################################

add_significant_terms <- function(data, terms, enforced_terms, num_shuffles = 100) {
  # Initialize variables to store the best model and its MSE
  best_model <- NULL
  best_summary <- NULL
  best_mse <- Inf
  
  # Perform the shuffling and model fitting multiple times
  for (i in 1:num_shuffles) {
    # Shuffle the terms
    shuffled_terms <- sample(terms)
    
    # Start with a formula including only the enforced terms and country fixed effects
    formula <- as.formula(
      paste("RESEMAVAL ~(", paste(enforced_terms, collapse = " + "), ")* ReligCountry"))
    #formula <- as.formula(
      #paste("RESEMAVAL ~", paste(enforced_terms, collapse = " + ")))
    
    # Fit the initial model
    regModel <- lm(formula, data)
    summary_reg <- summary(regModel)
    
    # Get initial p-values
    p_values <- summary_reg$coefficients[, 4]
    
    # Create a vector to store significant terms
    significant_terms <- enforced_terms
    
    # Iteratively add terms
    for (term in shuffled_terms) {
      # Update the formula by adding the current term
      new_formula <- as.formula(
        paste("RESEMAVAL ~(", paste(c(significant_terms, term), collapse = " + "), ")* ReligCountry"))
      #new_formula <- as.formula(
        #paste("RESEMAVAL ~", paste(c(significant_terms, term), collapse = " + ")))
      
      # Fit the model with the updated formula
      new_model <- lm(new_formula, data)
      new_summary <- summary(new_model)
      
      # Get p-values
      new_p_values <- new_summary$coefficients[, 4]
      
      # Check if adding the term causes any p-value to exceed 0.05
      if (all(new_p_values <= 0.05, na.rm = TRUE)) {
        # If not, add the term to the list of significant terms
        significant_terms <- c(significant_terms, term)
        
        # Update the model and summary
        regModel <- new_model
        summary_reg <- new_summary
        p_values <- new_p_values
      }
    }
    
    # Calculate the Mean Squared Error (MSE)
    mse <- mean(resid(regModel)^2)
    
    # Update the best model if the current one has a lower MSE
    if (mse < best_mse) {
      best_model <- regModel
      best_summary <- summary_reg
      best_mse <- mse
      aic <- AIC(regModel)
    }
  }
  
  # Return the best model and its summary
  return(list(Regression = best_model, Summary = best_summary, MSE = best_mse, AIC = aic))
}

terms <- c("DEFIANCE", 
           "SCEPTICISM", "I((10*(SCEPTICISM -0.5))^2)",
           "I((10*(DEFIANCE -0.5))^2)",
           "AGE","RELATIVISM", "I((10*(RELATIVISM -0.5))^2)",
           "I((AGE)^2)", "I(AGE*EDUCATION)", "EDUCATION", "INCOME_LEVEL")


terms_of_interest <- c("I((10*(DISBELIEF-0.5))^2)", "DISBELIEF")

# Generate all combinations of 2 or 3 terms
combinations <- c(
  combn(terms_of_interest, 2, simplify = FALSE),
  combn(terms_of_interest, 1, simplify = FALSE)
)

# Initialize variables to store the best result
best_result <- NULL
best_mse <- Inf
best_aic <- Inf

# Iterate over each combination of enforced terms
for (enforced_terms in combinations) {
  result <- add_significant_terms(data = WVS_20, terms = terms, enforced_terms = enforced_terms, num_shuffles = 50)
  
  print(paste("Enforced terms:", paste(enforced_terms, collapse = ", ")))
  print(result$Summary)
  print(result$MSE)
  print(result$AIC)
  
  # Update the best result based on MSE
  if (!is.na(result$MSE) && result$MSE < best_mse) {
    best_result <- result
    best_mse <- result$MSE
    best_aic <- result$AIC
  }
}

# Print the best result
print("Best Result:")
print(best_result$Summary)
print(best_result$MSE)
print(best_result$AIC)

formula <- RESEMAVAL ~ (DISBELIEF+ I((10 * (DISBELIEF - 0.5))^2)+RELATIVISM+ I((10 * (RELATIVISM - 0.5))^2)+ SCEPTICISM+I((10 * (SCEPTICISM - 0.5))^2)+ AGE+ I(AGE * EDUCATION)+ DEFIANCE)*ReligCountry
ols_model <- lm(formula, data = WVS_20)
summary(ols_model)

###############################################################################
# Residuals Homoscedasticity Check
###############################################################################

plot(ols_model$residuals~ols_model$fitted.values, col='#071952')

###############################################################################
# Quantile Regression
###############################################################################

tau_values <- c(0.1, 0.25, 0.5, 0.75, 0.9)
tau_values <- c(0.1, 0.5, 0.9)

qr_models <- lapply(tau_values, function(tau) rq(formula, data = WVS_20, tau = tau))

# Display the summary of each quantile regression model
qr_summaries <- lapply(qr_models, function(model) summary(model, se = "boot", R = 100, bmethod = "xy"))

#qr_summaries <- lapply(qr_models, function(model) summary(model, se = "boot"))
names(qr_summaries) <- paste0("tau_", tau_values)
qr_summaries

ols_coef <- coef(ols_model)
ols_coef_df <- as.data.frame(t(ols_coef))
ols_coef_df$tau <- "OLS"

qr_coefs <- lapply(qr_models, coef)
qr_coef_df <- do.call(rbind, lapply(1:length(qr_coefs), function(i) {
  data.frame(t(qr_coefs[[i]]), tau = tau_values[i])
}))
colnames(qr_coef_df) <- c(names(ols_coef), "tau")

# Combine OLS and quantile regression coefficients for comparison
coef_comparison <- rbind(ols_coef_df, qr_coef_df)

# Melt the data for plotting
coef_melt <- melt(coef_comparison, id.vars = "tau")

# Plot the coefficients
ggplot(coef_melt, aes(x = tau, y = value, color = variable, group = variable)) +
  geom_line() +
  labs(title = "OLS vs Quantile Regression Coefficients", x = "Quantile/OLS", y = "Coefficient") +
  theme_minimal()


# Calculate residuals for the OLS model
ols_residuals <- resid(ols_model)

# Calculate residuals for the median regression model
qr_residuals <- resid(qr_models[[3]])

# Plot residuals for OLS
ggplot(WVS_20, aes(x = predict(ols_model), y = ols_residuals)) +
  geom_point() +
  labs(title = "OLS Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_minimal()


# Plot residuals for quantile regression
ggplot(WVS_20, aes(x = predict(qr_models[[3]]), y = qr_residuals)) +
  geom_point() +
  labs(title = "Quantile Regression Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_minimal()


# Predict values for OLS and quantile regression
ols_predictions <- predict(ols_model)
qr_predictions <- sapply(qr_models, predict)

# Combine predictions with the original data
WVS_20_pred <- cbind(WVS_20, ols_predictions, qr_predictions)
colnames(WVS_20_pred)[(ncol(WVS_20)+1):(ncol(WVS_20)+length(tau_values)+1)] <- c("OLS", paste0("Q", tau_values))

# Melt the data for plotting
WVS_20_melt <- melt(WVS_20_pred, id.vars = colnames(WVS_20))

# Plot predictions
ggplot(WVS_20_melt, aes(x = RESEMAVAL, y = value, color = variable)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#071952") +
  scale_color_manual(values = c("OLS" = "#B31312", "Q0.1" = "#ECA869", "Q0.5" = "#B08BBB", "Q0.9" = "#AED2FF")) +
  labs(title = "OLS vs Quantile Regression Predictions", x = "Actual Values", y = "Predicted Values", color = "Model") +
  theme_minimal()

# Explanation:
# This plot compares the actual values of the response variable (medv) with the predicted values from 
# both the OLS and quantile regression models. It helps to visualize the performance of the models.

# Perform Wald tests for quantile regression coefficients
qr_tests <- lapply(qr_models, function(model) {
  summary(model,se = "boot", R = 100, bmethod = "xy")$coefficients
})

# Combine the coefficients and p-values into a data frame
qr_pvals_df <- do.call(rbind, lapply(1:length(qr_tests), function(i) {
  coefs <- qr_tests[[i]]
  data.frame(tau = tau_values[i], variable = rownames(coefs), coef = coefs[, 1], p_value = coefs[, 4])
}))

# Plot the p-values
ggplot(qr_pvals_df, aes(x = tau, y = p_value, color = variable, group = variable)) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(title = "P-values of Quantile Regression Coefficients", x = "Quantile", y = "P-value") +
  theme_minimal()

###################################################################
#################### Censored Data Regression #####################
###################################################################

###################################################################
# Import Packages
###################################################################
#install.packages("VGAM")
#install.packages("gridExtra") # Run this if the package is not installed
library(VGAM)
library(ggplot2)
library(gridExtra)

update.packages("lmtest")
update.packages("VGAM")
################################################################################
# Define a left-censoring threshold (e.g., RESEMAVAL < 0.1 will be censored)
### performs left censoring on the RESEMAVAL variable in the WVS_20 data frame
################################################################################

left_censor_threshold <- 0.2
WVS_20$RESEMAVAL_censored <- ifelse(WVS_20$RESEMAVAL < left_censor_threshold, left_censor_threshold, WVS_20$RESEMAVAL)

#################################################################################
# Tobit model Fit
##################################################################################

tobit_formula <- RESEMAVAL_censored ~ (DISBELIEF + I((10 * (DISBELIEF - 0.5))^2) + RELATIVISM + 
                                         I((10 * (RELATIVISM - 0.5))^2) + SCEPTICISM + 
                                         I((10 * (SCEPTICISM - 0.5))^2) + AGE + 
                                         I(AGE * EDUCATION) + DEFIANCE) * ReligCountry

censored_ols_model <- lm(tobit_formula, data = WVS_20)
summary(censored_ols_model)

tobit_model <- vglm(formula=tobit_formula, tobit(Lower = left_censor_threshold), data = WVS_20)
summary(tobit_model)

################################################################################
################################################################################
# OLS VS. Tobit
################################################################################
### AIC Comparison
### Residuals Comparison
### Prediction Comparison
################################################################################
################################################################################

ols_aic <- AIC(censored_ols_model)
tobit_aic <- AIC(tobit_model)

ols_residuals <- resid(censored_ols_model)
tobit_residuals_vec <- resid(tobit_model)
tobit_residuals <- tobit_residuals_vec[,"mu"]
tobit_sdErrors <- tobit_residuals_vec[,"loglink(sd)"]
#par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))
plot(ols_residuals ~ predict(censored_ols_model), main = "Censored OLS Residuals", xlab = "Fitted Values", ylab = "Residuals",col="#071952")
plot(tobit_residuals ~ predict(tobit_model), main = "Tobit Residuals", xlab = "Fitted Values", ylab = "Residuals",col="#1F6E8C")

censored_ols_predictions <- predict(censored_ols_model)
tobit_predictions <- predict(tobit_model)

WVS_20$censored_OLS <- censored_ols_predictions
WVS_20$Tobit <- tobit_predictions

ggplot(WVS_20, aes(x = RESEMAVAL, y = censored_OLS)) +
  geom_point(aes(color = "OLS"), alpha = 0.75) +
  geom_point(aes(x = RESEMAVAL, y =Tobit[,"mu"], color = "Tobit",), alpha = 0.25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#071952") +
  scale_color_manual(values = c("OLS" = "#8E05C2", "Tobit" = "#4E9F3D")) +
  labs(title = "OLS vs Tobit Predictions", x = "Actual Values", y = "Predicted Values", color = "Model") +
  theme_minimal()

ols_residuals_density <- ggplot(data.frame(Residuals = ols_residuals), aes(x = Residuals)) +
  geom_density(fill = "#8E05C2", alpha = 0.5,) +
  labs(title = "Censored OLS Residuals Distribution", x = "Residuals") +
  theme_minimal()

tobit_residuals_density <- ggplot(data.frame(Residuals = tobit_residuals), aes(x = Residuals)) +
  geom_density(fill = "#4E9F3D", alpha = 0.5) +
  labs(title = "Tobit Residuals Distribution", x = "Residuals") +
  theme_minimal()

grid.arrange(ols_residuals_density, tobit_residuals_density, ncol = 2)

ols_data <- data.frame(Residuals = ols_residuals, Model = "Censored OLS")
tobit_data <- data.frame(Residuals = tobit_residuals, Model = "Tobit")
#combined_data <- rbind(ols_data, tobit_data)
ggplot(rbind(ols_data, tobit_data), aes(x = Residuals, fill = Model)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Censored OLS" = "#8E05C2", "Tobit" = "#4E9F3D")) +
  labs(title = "Residuals Distribution: Censored OLS vs. Tobit", 
       x = "Residuals", 
       y = "Density") +
  theme_minimal()

# Arrange the residuals density plots
# residuals_density_plot <- grid.arrange(ols_residuals_density, tobit_residuals_density, nrow = 1)

# 3. Prediction Error Plot
WVS_20$OLS_Error <- WVS_20$censored_OLS - WVS_20$RESEMAVAL
WVS_20$Tobit_Error <- WVS_20$Tobit[,"mu"] - WVS_20$RESEMAVAL


qq_ols <- ggplot(data.frame(Residuals = ols_residuals), aes(sample = Residuals)) +
  stat_qq(data = data.frame(Residuals = ols_residuals), aes(sample = Residuals), color = "#8E05C2") +
  stat_qq_line(data = data.frame(Residuals = ols_residuals), aes(sample = Residuals), color = "#8E05C2") +
  labs(title = "QQ Plot of Censored OLS Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

qq_tobit <- ggplot(data.frame(Residuals = tobit_residuals), aes(sample = Residuals)) +
  stat_qq(data = data.frame(Residuals = tobit_residuals), aes(sample = Residuals), color = "#4E9F3D") +
  stat_qq_line(data = data.frame(Residuals = tobit_residuals), aes(sample = Residuals), color = "#4E9F3D") +
  labs(title = "QQ Plot of Tobit Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

grid.arrange(qq_ols, qq_tobit, ncol = 2)

#ggplot() +
#  stat_qq(data = data.frame(Residuals = ols_residuals), aes(sample = Residuals), color = "#8E05C2") +
#  stat_qq_line(data = data.frame(Residuals = ols_residuals), aes(sample = Residuals), color = "#8E05C2") +
#  stat_qq(data = data.frame(Residuals = tobit_residuals), aes(sample = Residuals), color = "#4E9F3D") +
#  stat_qq_line(data = data.frame(Residuals = tobit_residuals), aes(sample = Residuals), color = "#4E9F3D") +
#  labs(title = "QQ Plot of Residuals: Censored OLS vs. Tobit", 
#       x = "Theoretical Quantiles", 
#       y = "Sample Quantiles") +
#  theme_minimal()


ols_homoscedasticity <- ggplot(data.frame(Fitted = predict(censored_ols_model), Residuals = ols_residuals), 
                               aes(x = Fitted, y = Residuals)) +
  geom_point(color = "#8E05C2", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Censored OLS Residuals vs. Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

tobit_homoscedasticity <- ggplot(data.frame(Fitted = predict(tobit_model)[,"mu"], Residuals = tobit_residuals), 
                                 aes(x = Fitted, y = Residuals)) +
  geom_point(color = "#4E9F3D", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Tobit Residuals vs. Fitted Values", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

grid.arrange(ols_homoscedasticity, tobit_homoscedasticity, ncol = 2)


cat("Models Comparison:\n")
cat("OLS Model:\n")
summary(censored_ols_model)
cat("Tobit Model - Left censoring (< 0.2):\n")
summary(tobit_model)
cat("=========================================================================")
