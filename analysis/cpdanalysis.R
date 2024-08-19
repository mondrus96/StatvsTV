# This script is for analysis of change point information between CN/eMCI
library(ggplot2)

# Load change point data
load("../data/ADNI/FBS.rda")
load("../data/ADNI/ADNIlist.rda")
all_results = c()
for(i in 1:length(CPDlist)){
  # Select the change points
  cpts = CPDlist[[i]]$change_points
  cpts = cpts$T[order(cpts$stat_test)][1:2]
  rowvec = c(ADNIlist[[i]]$class, cpts, min(cpts))
  
  # Append to output
  all_results = rbind(all_results, rowvec)
}

# Plot all change points across CN/eMCI
CNcpts = c(all_results[all_results[,1] == 0, 2:3])
eMCIcpts = c(all_results[all_results[,1] == 1, 2:3])

# Prepare the data for CN changepoints
CN_table <- as.data.frame(table(CNcpts))
names(CN_table) <- c("Time", "Frequency")

# Prepare the data for eMCI changepoints
eMCI_table <- as.data.frame(table(eMCIcpts))
names(eMCI_table) <- c("Time", "Frequency")

# Convert Time to numeric for both datasets
CN_table$Time <- as.numeric(as.character(CN_table$Time))
eMCI_table$Time <- as.numeric(as.character(eMCI_table$Time))

# Plot CN changepoints as vertical lines
png("CNcpts.png")
ggplot(CN_table, aes(x = Time, y = Frequency)) +
  geom_segment(aes(xend = Time, yend = 0), color = "#0f7ab5") +
  labs(x = "Time", y = "Frequency", title = "CN Changepoints") +
  scale_x_continuous(breaks = seq(30, 100, by = 10)) +  # Adjust x-axis breaks
  theme_minimal() +
  ylim(0, 6) + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
dev.off()

# Plot eMCI changepoints as vertical lines
png("eMCIcpts.png")
ggplot(eMCI_table, aes(x = Time, y = Frequency)) +
  geom_segment(aes(xend = Time, yend = 0), color = "#0f7ab5") +
  labs(x = "Time", y = "Frequency", title = "eMCI Changepoints") +
  scale_x_continuous(breaks = seq(30, 100, by = 10)) +  # Adjust x-axis breaks
  theme_minimal() +
  ylim(0, 6) + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
dev.off()

# KS test result
paste0("CNcpts vs. eMCIcpts ALL: ", ks.test(CNcpts, eMCIcpts, exact = FALSE)$p.value) # overall different
firstCNcpts = c(all_results[all_results[,1] == 0, 4])
firsteMCIcpts = c(all_results[all_results[,1] == 1, 4])
t.test(firstCNcpts, firsteMCIcpts, alternative = "less") # first change point for CN < eMCI

# Save tables separately
CNcpts = all_results[all_results[,1] == 0, 2:3]
eMCIcpts = all_results[all_results[,1] == 1, 2:3]
rownames(CNcpts) = rownames(eMCIcpts) = NULL
colnames(CNcpts) = colnames(eMCIcpts) = c("1st cpt", "2nd cpt")
write.csv(CNcpts, "CNcpts.csv", row.names = FALSE)
write.csv(eMCIcpts, "eMCIcpts.csv", row.names = FALSE)