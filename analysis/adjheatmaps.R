# Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)
library(gridExtra)
library(cowplot)

# Load all the functions
filesrcs = list.files(path = "../functions", pattern="*.R$")
sapply(paste0("../functions/", filesrcs), source)

# Load results
load("../data/ADNI/FBS.rda")
load("../data/ADNI/ADNIlist.rda")

# Estimate networks and graphs
FCdata = est.FC(ADNIlist, CPDlist, 2)
FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)

# Append state information
FCdata$state = rep(1:3, times = 68)

# Group by 'class' and 'state' and calculate the mean for each edge
avg_df = FCdata %>%
  group_by(class, state) %>%
  summarise(across(starts_with("_V"), mean, na.rm = TRUE))

### 1 Plotting adjacency matrices ###
# Define your custom theme with no gridlines and axis labels every 10 nodes
custom_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),  # Rotate x-axis labels
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(), 
    axis.ticks = element_blank(),
    legend.position = "none"
  )

# Function to create a single heatmap plot without title and legend
create_heatmap_plot <- function(curr_row, num_nodes, midpoint, limit) {
  adj_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  adj_matrix[upper.tri(adj_matrix)] <- as.numeric(curr_row)
  adj_matrix <- adj_matrix + t(adj_matrix)
  diag(adj_matrix) <- 0
  adj_matrix_long <- melt(adj_matrix)
  
  ggplot(adj_matrix_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#0C6291", mid = "#FBFEF9", high = "#A63446", 
                         midpoint = midpoint, limit = limit, space = "Lab") +
    custom_theme +
    scale_x_continuous(breaks = seq(0, num_nodes, by = 10)) +
    scale_y_continuous(breaks = seq(0, num_nodes, by = 10))
}

# List to store plots
plot_list <- list()

# Loop through all rows in avg_df
midpoint = 0.2; limit = c(0, 1)
for(i in 1:nrow(avg_df)) {
  curr_row <- avg_df[i, -c(1, 2)]
  num_nodes <- (1 + sqrt(1 + 8 * length(curr_row))) / 2
  p <- create_heatmap_plot(curr_row, num_nodes, midpoint, limit)
  plot_list[[i]] <- p
}

# Create a dummy plot to extract the legend
dummy_plot <- ggplot(melt(matrix(runif(100), nrow=10)), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0C6291", mid = "#FBFEF9", high = "#A63446", 
                       midpoint = 0.2, limit = c(0, 1), space = "Lab", name = "Edge Weight") +
  theme_minimal()
legend <- get_legend(dummy_plot)

# Create row and column labels
row_labels <- c("CN", "eMCI")
col_labels <- c("Segment 1", "Segment 2", "Segment 3")

# Function to create labeled grobs
label_grob <- function(label, rot = 0) {
  textGrob(label, rot = rot, gp = gpar(fontsize = 15, fontface = "bold"))
}

# Arrange the plots in a grid
grid_plots <- arrangeGrob(
  arrangeGrob(
    textGrob(" "),  # Empty placeholder for top-left corner
    label_grob(col_labels[1]), label_grob(col_labels[2]), label_grob(col_labels[3]), 
    ncol = 4, 
    widths = unit(c(0.5, 7, 7, 7), "cm")  # Adjust the widths to make plots bigger
  ),
  arrangeGrob(
    label_grob(row_labels[1], rot = 90), plot_list[[1]], plot_list[[2]], plot_list[[3]], 
    ncol = 4, 
    widths = unit(c(0.5, 7, 7, 7), "cm")  # Adjust the widths to make plots bigger
  ),
  arrangeGrob(
    label_grob(row_labels[2], rot = 90), plot_list[[4]], plot_list[[5]], plot_list[[6]], 
    ncol = 4, 
    widths = unit(c(0.5, 7, 7, 7), "cm")  # Adjust the widths to make plots bigger
  ),
  nrow = 3,
  heights = unit(c(0.5, 7, 7), "cm")  # Adjust the heights to make plots bigger
)

# Combine the grid of plots with the legend
final_plot <- plot_grid(grid_plots, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# Save the final plot
ggsave("adjmats.png", plot = final_plot, width = 25, height = 15, units = "cm")
dev.off()

### Plot difference
avgdiff_df <- avg_df[1:3,] - avg_df[4:6,]
plot_list <- list()
midpoint = 0; limit = c(-0.3, 0.3)
for(i in 1:nrow(avgdiff_df)) {
  curr_row <- avgdiff_df[i, -c(1, 2)]
  num_nodes <- (1 + sqrt(1 + 8 * length(curr_row))) / 2
  p <- create_heatmap_plot(curr_row, num_nodes, midpoint, limit)
  plot_list[[i]] <- p
}

# Define grid of plots
grid_plots <- arrangeGrob(
  arrangeGrob(
    label_grob("Difference", rot = 90), plot_list[[1]], plot_list[[2]], plot_list[[3]], 
    ncol = 4, 
    widths = unit(c(0.5, 7, 7, 7), "cm")  # Adjust the widths to make plots bigger
  ),
  nrow = 1,
  heights = unit(c(7), "cm")  # Adjust the heights to make plots bigger
)

# Create a dummy plot to extract the legend
dummy_plot <- ggplot(melt(matrix(runif(100), nrow=10)), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0C6291", mid = "#FBFEF9", high = "#A63446", 
                       midpoint = 0, limit = c(-0.3, 0.3), space = "Lab", name = "Edge Weight") +
  theme_minimal()
legend <- get_legend(dummy_plot)

# Combine the grid of plots with the legend
final_plot <- plot_grid(grid_plots, legend, ncol = 2, rel_widths = c(0.9, 0.1))

# Save the final plot
ggsave("adjmatdiff.png", plot = final_plot, width = 25, height = 7, units = "cm")

### 2 Calculate squared differences ###
# Assuming avg_df contains columns 'class', 'state', and adjacency matrix entries
class_0 <- avg_df[avg_df$class == 0, ]
class_1 <- avg_df[avg_df$class == 1, ]

compute_squared_diff <- function(state, class_0, class_1) {
  # Extract the adjacency matrices for the given state
  adj_class_0 <- class_0[class_0$state == state, -c(1, 2)]
  adj_class_1 <- class_1[class_1$state == state, -c(1, 2)]
  
  # Compute the squared differences
  squared_diff <- (as.numeric(adj_class_0) - as.numeric(adj_class_1))^2
  
  # Convert to matrix form
  num_nodes <- (1 + sqrt(1 + 8 * length(squared_diff))) / 2
  adj_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
  adj_matrix[upper.tri(adj_matrix)] <- squared_diff
  adj_matrix <- adj_matrix + t(adj_matrix)
  diag(adj_matrix) <- 0
  
  return(adj_matrix)
}

# Compute squared differences for each state
squared_diff_state_1 <- compute_squared_diff(1, class_0, class_1)
squared_diff_state_2 <- compute_squared_diff(2, class_0, class_1)
squared_diff_state_3 <- compute_squared_diff(3, class_0, class_1)

# Function to sum the upper triangular elements
sum_upper_tri <- function(matrix) {
  return(sum(matrix[upper.tri(matrix)]))
}

# Sum the squared differences for the upper triangle
sum_squared_diff_state_1 <- sum_upper_tri(squared_diff_state_1)
sum_squared_diff_state_2 <- sum_upper_tri(squared_diff_state_2)
sum_squared_diff_state_3 <- sum_upper_tri(squared_diff_state_3)

print(sum_squared_diff_state_1)
print(sum_squared_diff_state_2)
print(sum_squared_diff_state_3)