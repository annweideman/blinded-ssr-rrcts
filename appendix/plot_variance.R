#--------------------------------------------------
# Plots of estimated variance for appendix
# Author: Ann Marie Weideman
# Date: 1/18/24
#--------------------------------------------------

library(ggplot2)
library(cowplot)
library(scales)  # For scientific notation

# Parameters
p1 <- seq(0, 1, by = 0.01)
n <- 500

# Variance calculations
binary_variance <- 1/n * p1 * (1 - p1) * (1 - 2 * p1)^2
tte_variance <- 1/n * (1/p1 - 1)

# Data frames for ggplot
binary_data <- data.frame(p1 = p1, Variance = binary_variance)
tte_data <- data.frame(p1 = p1, Variance = tte_variance)

# Common theme settings for axis text
common_theme <- theme(
  legend.position = 'none',
  plot.margin = margin(5.5, 5.5, 5.5, 5.5),
  axis.text = element_text(size = 10),
  axis.title.x = element_text(size = 14),  
  axis.title.y = element_text(size = 14)
)

# Plot for binary variance
binary_plot <- ggplot(binary_data, aes(x = p1, y = Variance)) +
  geom_line(color = "#4d85d1", size = 1) +
  labs(x = '', y = '') +
  scale_y_continuous(labels = scientific) +
  theme_minimal() +
  common_theme

# Plot for time-to-event Variance
tte_plot <- ggplot(tte_data, aes(x = p1, y = Variance)) +
  geom_line(color = "#085b67", size = 1) +
  labs(x = expression("Pooled event rate " ~ (hat(p)[1])), y = '') +  # Corrected hat syntax
  theme_minimal() +
  common_theme

# Combine the plots
combined_plots <- plot_grid(binary_plot, tte_plot, ncol = 1, align = 'v', 
                            axis = 'l', rel_heights = c(1, 1))

# Add a common y-axis label (left-aligned and centered)
y_axis_label <- ggdraw() + 
  draw_label("Estimated variance", x = 0, y = 0.4, angle = 90, 
             hjust = 0, vjust = 3, size = 14) +
  theme_void() +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

# Dummy data for legend
legend_data <- data.frame(x = c(1, 1), y = c(1, 1), 
                          group = c("Binary", "Time-to-Event"))

# Create a combined legend
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = group)) +
  geom_line(size = 1) +
  scale_color_manual("", values = 
                       c("Binary" = "#4d85d1", "Time-to-Event" = "#085b67")) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12)
  )

# Combine the plots with the shared legend
final_plot <- plot_grid(
  plot_grid(
    y_axis_label,
    combined_plots, 
    ncol = 2, rel_widths = c(0.05, 1), align = 'v'
  ),
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.1) # Adjust the relative heights if needed
)

# View the plot
print(final_plot)

# Save the final plot with combined legend
ggsave(file = paste0(script_path,"/output/combined_estimated_variances.jpg"),
       plot = final_plot, width = 6.43, height = 5.90, dpi = 600)
