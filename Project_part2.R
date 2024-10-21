# Load necessary libraries
library(plm)
library(ggplot2)
library(lmtest)
########################
###### Panel Data ######
########################
ADHD = read.csv("lab_treatment_data_30.csv")
hist(log(ADHD$hrv_mean))
hist(ADHD$rec_num)
hist(ADHD$med)
hist(ADHD$user)
hist(ADHD$rec_num)

ADHD <- subset(ADHD, rec_num %in% c(30,31,32)) # 481 records
ADHD <- subset(ADHD, med %in% c(0,1)) # 481 records

ADHD$hrv_mean <- log(ADHD$hrv_mean)

ADHD <- subset(ADHD, rec_num %in% c(26,27,28,29, 30,31,32)) # 811 records

num_rows <- nrow(ADHD)

# Print the result
print(num_rows)




numeric_cols = sapply(ADHD,is.numeric)
cor_ADHD = cor(ADHD[,unlist(numeric_cols)] ,method = "pearson")
corrplot(cor_ADHD, method = "number", type = "upper" )

# Removing High correlated independent variables 
ADHD <- ADHD %>%
  mutate(acc_std = (x_std + y_std + z_std) / 3)

ADHD <- ADHD %>%
  mutate(gyro_std = (gx_std + gy_std + gz_std) / 3)

ADHD = ADHD[!names(ADHD) %in% c("gx_std", "gy_std", "gz_std","x_std", "y_std", "x_mean", "y_mean")]
ADHD = ADHD[!names(ADHD) %in% c("rms", "acc_std", "gyro_std")]
ADHD = ADHD[!names(ADHD) %in% c("acc_std", "x_std", "y_std")]
ADHD = ADHD[!names(ADHD) %in% c("g_rms", "rms")]
#ADHD = ADHD[!names(ADHD) %in% c("g_rms", "rms", "acc_std")]
ADHD = ADHD[!names(ADHD) %in% c("time_bin_start", "hrv_std")]
ADHD = ADHD[!names(ADHD) %in% c("gyro_std")]
table(index(ADHD), useNA = "ifany")


panel_data = pdata.frame(ADHD, index = c("user", "med"))

fixed_effects_model <- plm(hrv_mean ~ x_mean + y_mean + z_mean + gx_mean + gy_mean + gz_mean + gx_zcr + gy_zcr + gz_zcr + g_rms, data = panel_data, model = "within")
transformed_parameters <- c(
  "x_mean",
  "y_mean",
  "z_mean",
  "gx_mean",
  "gy_mean",
  "gz_mean",
  "gx_zcr",
  "gy_zcr",
  "gz_zcr",
  "g_rms",
  
  "I(x_mean^2)",
  "I(y_mean^2)",
  "I(z_mean^2)",
  "I(gx_mean^2)",
  "I(gy_mean^2)",
  "I(gz_mean^2)",
  "I(gx_zcr^2)",
  "I(gy_zcr^2)",
  "I(gz_zcr^2)",
  "I(g_rms^2)",
  
  "I(log(x_mean + 1e-6))",
  "I(log(y_mean + 1e-6))",
  "I(log(z_mean + 1e-6))",
  "I(log(gx_mean + 1e-6))",
  "I(log(gy_mean + 1e-6))",
  "I(log(gz_mean + 1e-6))",
  "I(log(gx_zcr + 1e-6))",
  "I(log(gy_zcr + 1e-6))",
  "I(log(gz_zcr + 1e-6))",
  "I(log(g_rms + 1e-6))",
  
  "I(sqrt(abs(x_mean)))",
  "I(sqrt(abs(y_mean)))",
  "I(sqrt(abs(z_mean)))",
  "I(sqrt(abs(gx_mean)))",
  "I(sqrt(abs(gy_mean)))",
  "I(sqrt(abs(gz_mean)))",
  "I(sqrt(abs(gx_zcr)))",
  "I(sqrt(abs(gy_zcr)))",
  "I(sqrt(abs(gz_zcr)))",
  "I(sqrt(abs(g_rms)))",
  
  "I(1 / (x_mean + 1e-6))",
  "I(1 / (y_mean + 1e-6))",
  "I(1 / (z_mean + 1e-6))",
  "I(1 / (gx_mean + 1e-6))",
  "I(1 / (gy_mean + 1e-6))",
  "I(1 / (gz_mean + 1e-6))",
  "I(1 / (gx_zcr + 1e-6))",
  "I(1 / (gy_zcr + 1e-6))",
  "I(1 / (gz_zcr + 1e-6))",
  "I(1 / (g_rms + 1e-6))",
  
  "gyro_std",
  "acc_std",
  "gx_std",
  "gy_std",
  "gz_std",
  "x_std",
  "y_std",
  "z_std",
  
  "I(gyro_std^2)",
  "I(acc_std^2)",
  "I(gx_std^2)",
  "I(gy_std^2)",
  "I(gz_std^2)",
  "I(x_std^2)",
  "I(y_std^2)",
  "I(z_std^2)",
  
  "I(log(gyro_std + 1e-6))",
  "I(log(acc_std + 1e-6))",
  "I(log(gx_std + 1e-6))",
  "I(log(gy_std + 1e-6))",
  "I(log(gz_std + 1e-6))",
  "I(log(x_std + 1e-6))",
  "I(log(y_std + 1e-6))",
  "I(log(z_std + 1e-6))",
  
  "I(sqrt(abs(gyro_std)))",
  "I(sqrt(abs(acc_std)))",
  "I(sqrt(abs(gx_std)))",
  "I(sqrt(abs(gy_std)))",
  "I(sqrt(abs(gz_std)))",
  "I(sqrt(abs(x_std)))",
  "I(sqrt(abs(y_std)))",
  "I(sqrt(abs(z_std)))",
  
  "I(1 / (gyro_std + 1e-6))",
  "I(1 / (acc_std + 1e-6))",
  "I(1 / (gx_std + 1e-6))",
  "I(1 / (gy_std + 1e-6))",
  "I(1 / (gz_std + 1e-6))",
  "I(1 / (x_std + 1e-6))",
  "I(1 / (y_std + 1e-6))",
  "I(1 / (z_std + 1e-6))"
)


# Loop through each transformed parameter and fit the fixed effects model
results <- vector("list", length(transformed_parameters))

# Loop through each transformed parameter and fit the fixed effects model
for (i in seq_along(transformed_parameters)) {
  param <- transformed_parameters[i]
  formula <- as.formula(paste("hrv_mean ~", param))
  fd_model <- plm(formula, data = panel_data, model = "fd")
  
  # Store results
  results[[i]] <- list(param = param, p_value = summary(fd_model)$coefficients[, "Pr(>|t|)"][2])
}

# Extract p-values from results
p_values <- sapply(results, function(x) x$p_value)

# Find the indices of the top 5 models based on p-values
top_5_indices <- order(p_values)[1:80]

# Print the top 5 results
for (i in top_5_indices) {
  result <- results[[i]]
  cat("Parameter:", result$param, "\n")
  cat("P-value:", result$p_value, "\n\n")
}




# Define your terms
terms <- c(
  "I(log(gz_mean + 1e-6))",
  "z_mean",
  "I(sqrt(abs(gx_mean)))",
  "I(log(gz_zcr + 1e-6))",
  "I(1 / (gy_mean + 1e-6))",
  "I(sqrt(abs(acc_std)))",
  "rec_num"
)

# Function to get combinations of a specified size
fit_plm <- function(formula, data) {
  plm(formula, data = data, model = "fd")
}

# Get all combinations from size 2 and above
all_combinations <- lapply(2:length(terms), function(size) combn(terms, size, simplify = FALSE))

best_mse <- Inf
# Fit models for each combination
for (comb_list in all_combinations) {
  for (comb in comb_list) {
    formula_str <- paste("hrv_mean ~", paste(comb, collapse = " + "))
    formula <- as.formula(formula_str)
    result <- plm(formula, data = panel_data, model = "fd")
    aa <- summary(result)$coefficients[-1, "Pr(>|t|)"]
    p_values <- as.numeric(sub(".*\\s+", "", aa))
    significance_level <- 0.05
    
    # Check if all p-values are less than the significance level
    all_significant <- all(p_values < significance_level)
    
    # Print or store results as needed
    mse <- mean(resid(result)^2)
    if (mse < best_mse && all_significant) {
      best_model <- result
      best_summary <- summary(result)
      best_mse <- mse
    }
  }
}



fd_model <- plm(hrv_mean~ I(log(gz_zcr + 1e-6)) +I(1 / (gy_mean + 1e-6))+I(sqrt(abs(gyro_std)))  +I(sqrt(abs(z_std))), data = panel_data, model = "fd")
summary(fd_model)
mse <- mean(resid(fd_model)^2)


base_formula <- formula()

for (param in transformed_parameters) {
  formula_str <- paste("hrv_mean ~ rec_num + I(1 / (gy_mean + 1e-6)) + I(sqrt(abs(acc_std)))", "+", param)
  formula <- as.formula(formula_str)
  model_result <- plm(formula, data = panel_data, model = "fd")
  mse <- mean(resid(result)^2)
  aa <- summary(model_result)$coefficients[-1, "Pr(>|t|)"]
  p_values <- as.numeric(sub(".*\\s+", "", aa))
  significance_level <- 0.05
  
  # Check if all p-values are less than the significance level
  all_significant <- all(p_values < significance_level)
  
  if (all_significant) {
    best_model <- model_result
    best_summary <- summary(model_result)
    best_mse <-  mean(resid(model_result)^2)}
}

best_mse <- Inf
for (param in transformed_parameters) {
  formula_str <- paste("hrv_mean ~ rec_num + I(1 / (gy_mean + 1e-6)) + I(sqrt(abs(acc_std)))", "+", param)
  formula <- as.formula(formula_str)
  model_result <- plm(formula, data = ADHD,index = c("user", "med"), 
                      model = "random")
  mse <- mean(resid(model_result)^2)
  # Check if all p-values are less than the significance level
  all_significant <- all(p_values < significance_level)
  
  if (mse < best_mse) {
    best_model <- model_result
    best_summary <- summary(model_result)
    best_mse <-  mean(resid(model_result)^2)}
}

aa <- summary(fd_model)$coefficients[-1, "Pr(>|t|)"]
p_values <- as.numeric(sub(".*\\s+", "", aa))

# Set the significance level
significance_level <- 0.05

# Check if all p-values are less than the significance level
all_significant <- all(p_values < significance_level)
#####################################
# Estimate a first Difference model #
#####################################

fd_model = plm(hrv_mean ~  I(1/(gy_mean + 1e-06)) + I(sqrt(abs(acc_std))) + rec_num, data = panel_data, model = "fd")
summary(fd_model)


re_model <- plm(hrv_mean ~ I(1/(gy_mean + 1e-06)) + I(sqrt(abs(acc_std))) + rec_num, 
                data = ADHD, 
                index = c("user", "med"), 
                model = "random")

# Summary of the model
summary(re_model)

## filter by 30-32

fd_model = plm(hrv_mean ~  I(x_std^2)  + rec_num +I(log(gy_std + 1e-6)) , data = panel_data, model = "fd")
summary(fd_model)


re_model <- plm(hrv_mean ~  I(x_std^2) +I(log(gy_std + 1e-6))+ I(log(z_mean + 1e-06)), 
                data = ADHD, 
                index = c("user", "med"), 
                model = "random")

# Summary of the model
summary(re_model)


# when filter from 26-32
re_model <- plm(hrv_mean ~  I(sqrt(abs(acc_std))) +rec_num+ gx_std, 
                data = ADHD, 
                index = c("user", "med"), 
                model = "random")

# Summary of the model
summary(re_model)



cat("Coefficients:\n")
print(coef(re_model))

# Conduct Hausman Test to compare with fixed effects model
hausman_test <- phtest(re_model, plm(hrv_mean ~  I(sqrt(abs(acc_std))) +rec_num+ gx_std, 
                                     data = ADHD, 
                                     index = c("user", "med"), 
                                     model = "within"))
print(hausman_test)


residuals <- resid(re_model)
plot(ADHD$med, residuals, main = "Scatter Plot of Residuals from Random Effects Model",
     xlab = "Year", ylab = "Residuals", pch = 20, col = "blue")
abline(h = 0, col = "red")






#FE for 30-32
formula <- hrv_mean ~ I(gy_std^2) + I(log(gx_mean + 1e-06)) + med

# Estimate the fixed effects model

fixed_effects_model <- plm(formula, data = ADHD, index = c("user", "med"), model = "within")
# Summary of the fixed effects model to view coefficients and statistical significance

summary(fixed_effects_model)

# Statistical Test for Fixed Effects
# Conduct the Hausman test
# First, estimate the random effects model for comparison

random_effects_model <- plm(hrv_mean ~ I(gy_std^2) + I(log(gx_mean + 1e-06)), data = ADHD, index = c("user", "med"), model = "random")

# Perform the Hausman test to compare fixed and random effects models

hausman_test <- phtest(fixed_effects_model, random_effects_model)

# Display the results of the Hausman test

hausman_test

fixed_effects <- fixef(fixed_effects_model)

# Convert the fixed effects to a data frame for plotting

fixed_effects_df <- data.frame(State = names(fixed_effects), Effect = as.numeric(fixed_effects))

# Plot the fixed effects by state to visualize the individual-specific effects

ggplot(fixed_effects_df, aes(x = reorder(State, Effect), y = Effect)) +
  geom_point() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Fixed Effects by State", x = "State", y = "Fixed Effect")


residuals <- resid(fixed_effects_model)
residuals_df <- data.frame(Residuals = residuals)
ggplot(residuals_df, aes(x = Residuals)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency")

for (i in seq_along(transformed_parameters)) {
  param <- transformed_parameters[i]
  formula <- as.formula(paste("hrv_mean ~ med +", param))
  fd_model <- plm(formula, data = ADHD, index = c("user", "med"), model = "within")
  
  # Store results
  results[[i]] <- list(param = param, p_value = summary(fd_model)$coefficients[, "Pr(>|t|)"][2])
}

# Extract p-values from results
p_values <- sapply(results, function(x) x$p_value)

# Find the indices of the top 5 models based on p-values
top_5_indices <- order(p_values)[1:80]

# Print the top 5 results
for (i in top_5_indices) {
  result <- results[[i]]
  cat("Parameter:", result$param, "\n")
  cat("P-value:", result$p_value, "\n\n")
}


best_mse <- Inf
for (param in transformed_parameters) {
  formula_str <- paste("hrv_mean ~ med + I(sqrt(abs(acc_std))) + rec_num + x_mean + I(y_mean^2) + I(log(z_mean + 1e-06)) + I(sqrt(abs(x_mean))) + z_mean+ I(1/(z_mean + 1e-06))", "+", param)
  formula <- as.formula(formula_str)
  model_result <- plm(formula, data = ADHD,index = c("user", "med"), 
                      model = "within")
  mse <- mean(resid(model_result)^2)
  aa <- summary(model_result)$coefficients[, "Pr(>|t|)"]
  p_values <- as.numeric(sub(".*\\s+", "", aa))
  
  # Set the significance level
  significance_level <- 0.05
  # Check if all p-values are less than the significance level
  all_significant <- all(p_values < significance_level)
  
  if (mse < best_mse && all_significant) {
    best_model <- model_result
    best_summary <- summary(model_result)
    best_mse <-  mean(resid(model_result)^2)}
}


#FE for 25-32
formula <- hrv_mean ~ med + I(sqrt(abs(acc_std))) + rec_num + x_mean + I(y_mean^2)

# Estimate the fixed effects model

fixed_effects_model <- plm(formula, data = ADHD, index = c("user", "med"), model = "within")
# Summary of the fixed effects model to view coefficients and statistical significance

summary(fixed_effects_model)

# Statistical Test for Fixed Effects
# Conduct the Hausman test
# First, estimate the random effects model for comparison

random_effects_model <- plm(hrv_mean ~ I(sqrt(abs(acc_std))) + rec_num + x_mean + I(y_mean^2), data = ADHD, index = c("user", "med"), model = "random")
summary(random_effects_model)

# Perform the Hausman test to compare fixed and random effects models

hausman_test <- phtest(fixed_effects_model, random_effects_model)

# Display the results of the Hausman test

hausman_test


formula <- hrv_mean ~ med + I(sqrt(abs(acc_std))) + rec_num + x_mean + I(y_mean^2) + I(log(z_mean + 1e-06)) + I(sqrt(abs(x_mean))) + z_mean+ I(1/(z_mean + 1e-06))

# Estimate the fixed effects model

fixed_effects_model <- plm(formula, data = ADHD, index = c("user", "med"), model = "within")

fixed_effects <- fixef(fixed_effects_model)

# Convert the fixed effects to a data frame for plotting

fixed_effects_df <- data.frame(State = names(fixed_effects), Effect = as.numeric(fixed_effects))

# Plot the fixed effects by state to visualize the individual-specific effects

ggplot(fixed_effects_df, aes(x = reorder(State, Effect), y = Effect)) +
  geom_point() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Fixed Effects by State", x = "State", y = "Fixed Effect")


residuals <- resid(fixed_effects_model)
residuals_df <- data.frame(Residuals = residuals)
ggplot(residuals_df, aes(x = Residuals)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Residuals", x = "Residuals", y = "Frequency")



