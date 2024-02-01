#--------------------------------------------------
# Plots of bias associated with the variance
# estimators for appendix
# Author: Ann Marie Weideman
# Date: 1/17/24
#--------------------------------------------------

library(ggplot2)
library(rstudioapi)

# Set working directory to current location
script_path<-setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#-----------------------------
# Function Definitions
#-----------------------------

# Function to calculate percent bias for binary endpoint
binary_bias <- function(n, p1) {
  
  e1 = p1
  e2 = (p1 * (1 - p1)) / n + p1^2
  e3 = (1 / n^3) * (n * p1 + 3 * n * (n - 1) * p1^2 + n * (n - 1) * (n - 2) * p1^3)
  e4 = (1 / n^4) * (n * p1 + 7 * n * (n - 1) * p1^2 + 6 * n * (n - 1) * (n - 2) * p1^3+
                      n*(n-1)*(n-2)*(n-3)*p1^4)
  
  truth <- 1/n*(p1 * (1 - p1) * (1 - 2 * p1)^2)
  sd_truth <- sqrt(truth)
  estimate <- 1/n*(-4 * e4 + 8 * e3 - 5 * e2 + e1)
  sd_estimate <- sqrt(estimate)
  bias <- abs(sd_estimate - sd_truth)
  pct_bias <- bias / sd_truth * 100
  
  return(round(pct_bias, 2))
}

# Function to calculate percent bias for time-to-event dndpoint
TTE_bias <- function(n, p1) {

  sum1 <- sum(1 / (1:n) * choose(n, 1:n) * p1^(1:n) * (1 - p1)^(n - (1:n)))
  
  truth <- 1/n*(1 / p1 - 1)
  sd_truth <- exp(1.96*sqrt(truth))
  estimate <- 1/n*(n * sum1 - 1)
  sd_estimate <- exp(1.96*sqrt(estimate))
  bias <- sd_estimate - sd_truth
  pct_bias <- bias / sd_truth * 100
  
  return(round(pct_bias, 2))
}

#-----------------------------
# Data Preparation
#-----------------------------

p1 <- 0.1  # Pooled event rate
n_values <- seq(200, 1000, by = 100) # Sample size

# Binary Data
binary_data <- data.frame(
  n = n_values,
  percent_bias = sapply(n_values, binary_bias, p1 = p1),
  Endpoint = "Binary"
)

# Time-to-Event Data
TTE_data <- data.frame(
  n = n_values,
  percent_bias = sapply(n_values, TTE_bias, p1 = p1),
  Endpoint = "Time-to-Event"
)

# Combined Data
combined_data <- rbind(binary_data, TTE_data)

#-----------------------------
# Plotting
#-----------------------------

# Define colors for accessibility
colors <- c("Binary" = "#4d85d1", "Time-to-Event" = "#085b67")

# Create combined plot
p1 <- ggplot(combined_data, aes(x = n, y = percent_bias, color = Endpoint)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = percent_bias), vjust = -2.5, size = 3, hjust = 0.25) +
  scale_color_manual(values = colors) +
  labs(x = expression("Sample size at interim analysis (" * n["int"] * ")"), 
       y = expression("Bias of 95% CI half-width (%)"), color = "Endpoint Type") +
  scale_x_continuous(breaks = seq(200, 1000, by = 100), 
                     expand = c(0, 0), limits = c(175, 1025)) +
  scale_y_continuous(limits = c(0, max(combined_data$percent_bias) * 1.1)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "black"),
    strip.text.x = element_text(size = 12)
  )

print(p1)

# Save the plot
ggsave(file = paste0(script_path,"/output/combined_variance_estimator_bias.jpg"),
       plot = p1, width = 6.92, height = 5.90, dpi = 600)
