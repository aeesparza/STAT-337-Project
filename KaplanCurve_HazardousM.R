#Code for Kaplan Curves and Hazardous Model
# Load necessary libraries
if (!require(survival)) {
    install.packages("survival")
    library(survival)
}

# Load ARFF file
bone_marrow_data <- read.arff("bone-marrow.arff")

# View the first few rows and general structure of the dataset
head(bone_marrow_data)
str(bone_marrow_data)
summary(bone_marrow_data)

# Basic correlation matrix for numeric variables and a scatterplot matrix of the first four numeric variables
numeric_vars <- bone_marrow_data[, sapply(bone_marrow_data, is.numeric)]
correlations <- cor(numeric_vars)
print(correlations)
pairs(numeric_vars[, 1:4])  # Adjust indices according to specific needs

# Create survival objects from survival time and status
surv_obj <- Surv(time = bone_marrow_data$survival_time, event = bone_marrow_data$survival_status)

# Fit Kaplan-Meier survival curve and plot it
km_fit <- survfit(surv_obj ~ 1)  # Overall survival curve
plot(km_fit, main = "Kaplan-Meier Overall Survival Curve", xlab = "Time (days)", ylab = "Survival probability")

# Fit Cox proportional hazards model
cox_model <- coxph(surv_obj ~ CD34kgx10d6 + Recipientageint + Disease + Riskgroup + HLAmatch + aGvHDIIIIV + extcGvHD, data = bone_marrow_data)
summary(cox_model)

# Fit Kaplan-Meier curves stratified by key factors and plot them
km_fit_cd34 <- survfit(surv_obj ~ CD34kgx10d6, data = bone_marrow_data)
plot(km_fit_cd34, col = 1:6, lty = 1, lwd = 2, main = "Kaplan-Meier Curve by CD34+ Dosage",
     xlab = "Time (days)", ylab = "Survival Probability")
legend("topright", legend = names(km_fit_cd34$strata), col = 1:6, lty = 1, lwd = 2, title = "CD34+ Dosage")

# Kaplan-Meier for other key factors
plot_factors <- c("Recipientgender", "Stemcellsource", "Donorage", "Disease", "Riskgroup", "Txpostrelapse", "HLAmatch", "aGvHDIIIIV", "extcGvHD")
colors <- c("blue", "red", "green", "purple", "orange", "cyan")
for (factor in plot_factors) {
    km_fit_factor <- survfit(surv_obj ~ get(factor), data = bone_marrow_data)
    plot(km_fit_factor, col = colors, lty = 1, lwd = 2, main = paste("KM Curve for", factor),
         xlab = "Time (days)", ylab = "Survival Probability")
    legend("topright", legend = names(km_fit_factor$strata), col = colors, lty = 1, lwd = 2, title = factor)
}
# Assuming surv_obj and fitting have already been performed:
km_fit_cd34 <- survfit(surv_obj ~ CD34kgx10d6, data = bone_marrow_data)

# Plot the Kaplan-Meier survival curve for CD34+ cell dosage
plot(km_fit_cd34, main = "Kaplan-Meier Survival Curve for CD34+ Cells", xlab = "CD34+ Cells Dosage", ylab = "Survival Probability")

p_value_cd34 <- survdiff(surv_obj ~ CD34kgx10d6, data = bone_marrow_data)$p.value

# Adding the p-value to the plot
text(x = max(bone_marrow_data$CD34kgx10d6) * 0.8, y = 0.2, labels = paste("p =", format(p_value_cd34, digits = 3)), cex = 1.2, col = "red")
