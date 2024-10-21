
###################################################################
#################### Censored Data Regression #####################
###################################################################

###################################################################
# Import Packages
###################################################################
#install.packages("VGAM")
#install.packages("gridExtra") # Run this if the package is not installed
library(readxl)
library(car)
library(dplyr)
library(fitdistrplus)
library(corrplot)
library(fastDummies)
library(reshape2)
library(quantreg)
library(VGAM)
library(ggplot2)
library(gridExtra)
###################################################################
# Import Data set - The World Values Survey Wave 7
###################################################################

WVS_20 = read.csv("wvs20.csv")
head(WVS_20)

###############################################################################
# Covariance & Correlation 
###############################################################################

numeric_cols = sapply(WVS_20,is.numeric)
cor_WVS = cor(WVS_20[,unlist(numeric_cols)] ,method = "pearson")
corrplot(cor_WVS, method = "number", type = "upper" )

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
# Removing High correlated independent variables 
###############################################################################

WVS_20 = WVS_20[!names(WVS_20) %in% c("SACSECVAL", "G_TOWNSIZE", "H_URBRURAL")]

###############################################################################
# Covariance & Correlation (Highly correlated independent variables removed)
###############################################################################

numeric_cols = sapply(WVS_20,is.numeric)
cor_WVS = cor(WVS_20[,unlist(numeric_cols)] ,method = "pearson")
corrplot(cor_WVS, method = "number", type = "upper" )


################################################################################
# Define a left-censoring threshold (e.g., RESEMAVAL < 0.1 will be censored)
### performs left censoring on the RESEMAVAL variable in the WVS_20 data frame
################################################################################

left_censor_threshold1 <- 0.16
left_censor_threshold2 <- 0.32
right_censor_threshold1 <- 0.7
right_censor_threshold2 <- 0.8

left_censor_threshold <- left_censor_threshold2
right_censor_threshold <- right_censor_threshold2

WVS_20$RESEMAVAL_censored_left <- ifelse(WVS_20$RESEMAVAL < left_censor_threshold, left_censor_threshold, WVS_20$RESEMAVAL)
WVS_20$RESEMAVAL_censored_right <- ifelse(WVS_20$RESEMAVAL > right_censor_threshold, right_censor_threshold, WVS_20$RESEMAVAL)
#################################################################################
# Tobit model Fit
##################################################################################

tobit_formula_left <- RESEMAVAL_censored_left ~ (DISBELIEF + I((10 * (DISBELIEF - 0.5))^2) + RELATIVISM + 
                                         I((10 * (RELATIVISM - 0.5))^2) + SCEPTICISM + 
                                         I((10 * (SCEPTICISM - 0.5))^2) + AGE + 
                                         I(AGE * EDUCATION) + DEFIANCE) * ReligCountry
tobit_formula_right <- RESEMAVAL_censored_right ~ (DISBELIEF + I((10 * (DISBELIEF - 0.5))^2) + RELATIVISM + 
                                                    I((10 * (RELATIVISM - 0.5))^2) + SCEPTICISM + 
                                                    I((10 * (SCEPTICISM - 0.5))^2) + AGE + 
                                                    I(AGE * EDUCATION) + DEFIANCE) * ReligCountry
tobit_formula <- tobit_formula_left

censored_ols_model <- lm(tobit_formula, data = WVS_20)
summary(censored_ols_model)

#tobit_model_right <- vglm(formula=tobit_formula, tobit(Upper = right_censor_threshold), data = WVS_20)
tobit_model_left <- vglm(formula=tobit_formula, tobit(Lower = left_censor_threshold), data = WVS_20)
tobit_model <- tobit_model_left
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

#plot(ols_residuals ~ predict(censored_ols_model), main = "Censored OLS Residuals", xlab = "Fitted Values", ylab = "Residuals",col="#071952")
#plot(tobit_residuals ~ predict(tobit_model), main = "Tobit Residuals", xlab = "Fitted Values", ylab = "Residuals",col="#1F6E8C")

censored_ols_predictions <- predict(censored_ols_model)
tobit_predictions <- predict(tobit_model)
WVS_20$censored_OLS_pred <- censored_ols_predictions
WVS_20$tobit_pred <- tobit_predictions[,"mu"]



ggplot(WVS_20, aes(x = RESEMAVAL, y = censored_OLS_pred)) +
  geom_point(aes(color = "OLS"), alpha = 0.75) +
  geom_point(aes(x = RESEMAVAL, tobit_pred, color = "Tobit",), alpha = 0.25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#071952") +
  scale_color_manual(values = c("OLS" = "#6DA9E4", "Tobit" = "#FBA834")) +
  labs(title = "Predictions VS Actual Values - Tobit & Censored Data OLS", x = "Actual Values", y = "Predicted Values", color = "Model") +
  theme_minimal()



ols_residuals_density <- ggplot(data.frame(Residuals = ols_residuals), aes(x = Residuals)) +
  geom_density(fill = "#6DA9E4", alpha = 0.5,) +
  labs(title = "Censored OLS Residuals Distribution", x = "Residuals") +
  theme_minimal()

tobit_residuals_density <- ggplot(data.frame(Residuals = tobit_residuals), aes(x = Residuals)) +
  geom_density(fill = "#FBA834", alpha = 0.5) +
  labs(title = "Tobit Residuals Distribution", x = "Residuals") +
  theme_minimal()

grid.arrange(ols_residuals_density, tobit_residuals_density, ncol = 2)

ols_data <- data.frame(Residuals = ols_residuals, Model = "Censored OLS")
tobit_data <- data.frame(Residuals = tobit_residuals, Model = "Tobit")

ggplot(rbind(ols_data, tobit_data), aes(x = Residuals, fill = Model)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Censored OLS" = "#6DA9E4", "Tobit" = "#FBA834")) +
  labs(title = "Residuals Distribution: Censored OLS vs. Tobit", 
       x = "Residuals", 
       y = "Density") +
  theme_minimal()


qq_ols <- ggplot(data.frame(Residuals = ols_residuals), aes(sample = Residuals)) +
  stat_qq(data = data.frame(Residuals = ols_residuals), aes(sample = Residuals), color = "#6DA9E4") +
  stat_qq_line(data = data.frame(Residuals = ols_residuals), aes(sample = Residuals), color = "#6DA9E4") +
  labs(title = "QQ Plot of Censored OLS Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

qq_tobit <- ggplot(data.frame(Residuals = tobit_residuals), aes(sample = Residuals)) +
  stat_qq(data = data.frame(Residuals = tobit_residuals), aes(sample = Residuals), color = "#FBA834") +
  stat_qq_line(data = data.frame(Residuals = tobit_residuals), aes(sample = Residuals), color = "#FBA834") +
  labs(title = "QQ Plot of Tobit Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

grid.arrange(qq_ols, qq_tobit, ncol = 2)


# 3. Prediction Error Plot
WVS_20$OLS_Error <- WVS_20$RESEMAVAL - WVS_20$censored_OLS_pred
WVS_20$Tobit_Error <- WVS_20$RESEMAVAL - WVS_20$tobit_pred

ggplot(WVS_20, aes(x = censored_OLS_pred, y = OLS_Error)) +
  geom_point(aes(color = "OLS"), alpha = 0.75) +
  geom_point(aes(x = tobit_pred, y=Tobit_Error, color = "Tobit",), alpha = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", color = "#071952") +
  scale_color_manual(values = c("OLS" = "#6DA9E4", "Tobit" = "#FBA834")) +
  labs(title = "OLS vs Tobit Homoscedasticity", x = "Predicted Values", y = "Standart Error", color = "Model") +
  theme_minimal()

ols_homoscedasticity <- ggplot(data.frame(Fitted = WVS_20$censored_OLS_pred, Residuals = WVS_20$OLS_Error), 
                               aes(x = Fitted, y = Residuals)) +
  geom_point(color = "#6DA9E4",fill='#AED2FF', alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed",color = "#071952") +
  labs(title = "Censored OLS Homoscedasticity", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

# Create Tobit homoscedasticity plot
tobit_homoscedasticity <- ggplot(data.frame(Fitted = WVS_20$tobit_pred, Residuals = WVS_20$Tobit_Error), 
                                 aes(x = Fitted, y = Residuals)) +
  geom_point(color = "#FBA834",fill='#FFE6C7', alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed",color = "#071952") +
  labs(title = "Tobit Homoscedasticity", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

# Display both plots side by side
grid.arrange(ols_homoscedasticity, tobit_homoscedasticity, ncol = 2)


#cat("Models Comparison:\n")
#cat("OLS Model:\n")
#summary(censored_ols_model)
#cat("Tobit Model - Left censoring (< 0.2):\n")
#summary(tobit_model)
#cat("=========================================================================")

