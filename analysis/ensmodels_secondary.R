# This script is for plotting combinations of predictions
library(ggplot2)
library(scales)

# Combine all the data
load("../classifiers/data2/FBS_cpDFC2.rda")
CVeval$Name = "cpd2"
cols_to_convert <- setdiff(colnames(CVeval), c("Model", "Name"))
CVeval[cols_to_convert] <- lapply(CVeval[cols_to_convert], as.numeric)
# Calculate SE for each of the quantities
CVeval$AccSE = sqrt((CVeval$Accuracy*(1-CVeval$Accuracy))/68)
CVeval$SenSE = sqrt((CVeval$Sensitivity*(1-CVeval$Sensitivity))/68)
CVeval$SpeSE = sqrt((CVeval$Specificity*(1-CVeval$Specificity))/68)
CVeval$F1SE = sqrt((CVeval$F1*(1-CVeval$F1))/68)
df = CVeval
for(i in 2:10){
  load(paste0("../classifiers/ensdata2/", i, "ens.rda"))
  CVeval[cols_to_convert] <- lapply(CVeval[cols_to_convert], as.numeric)
  CVeval$Name = paste0("top", i)
  # Calculate SE for each of the quantities
  CVeval$AccSE = sqrt((CVeval$Accuracy*(1-CVeval$Accuracy))/68)
  CVeval$SenSE = sqrt((CVeval$Sensitivity*(1-CVeval$Sensitivity))/68)
  CVeval$SpeSE = sqrt((CVeval$Specificity*(1-CVeval$Specificity))/68)
  CVeval$F1SE = sqrt((CVeval$F1*(1-CVeval$F1))/68)
  df = rbind(df, CVeval)
}
df = df[df$Model == "linsvm",]

# Define the custom theme
custom_theme <- theme_minimal() +
  theme(
    text = element_text(size = rel(5)),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# Plot Accuracy
png("EnsAccuracy_secondary.png", width = 10, height = 5, units = "in", res = 300, pointsize = 4)
df$Name <- factor(df$Name, levels = df$Name)
ggplot(data = df, aes(x = Name, y = Accuracy)) +
  geom_point(size = 3, color = "#0f7ab5") +
  geom_ribbon(aes(ymin = Accuracy - AccSE, ymax = Accuracy + AccSE, group = 1), fill = "#0f7ab5", alpha = 0.2) +
  geom_line(aes(group = 1), color = "#0f7ab5") +
  scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits = c(0.5, 1), oob = rescale_none) +
  custom_theme +
  xlab("")
dev.off()

# Plot Sensitivity
png("EnsSensitivity_secondary.png", width = 10, height = 5, units = "in", res = 300, pointsize = 4)
df$Name <- factor(df$Name, levels = df$Name)
ggplot(data = df, aes(x = Name, y = Sensitivity)) +
  geom_point(size = 3, color = "#0f7ab5") +
  geom_ribbon(aes(ymin = Sensitivity - SenSE, ymax = Sensitivity + SenSE, group = 1), fill = "#0f7ab5", alpha = 0.2) +
  geom_line(aes(group = 1), color = "#0f7ab5") +
  scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits = c(0.5, 1), oob = rescale_none) +
  custom_theme +
  xlab("")
dev.off()

# Plot Specificity
png("EnsSpecificity_secondary.png", width = 10, height = 5, units = "in", res = 300, pointsize = 4)
df$Name <- factor(df$Name, levels = df$Name)
ggplot(data = df, aes(x = Name, y = Specificity)) +
  geom_point(size = 3, color = "#0f7ab5") +
  geom_ribbon(aes(ymin = Specificity - SpeSE, ymax = Specificity + SpeSE, group = 1), fill = "#0f7ab5", alpha = 0.2) +
  geom_line(aes(group = 1), color = "#0f7ab5") +
  scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits = c(0.5, 1), oob = rescale_none) +
  custom_theme +
  xlab("")
dev.off()

# Plot F1
png("EnsF1_secondary.png", width = 10, height = 5, units = "in", res = 300, pointsize = 4)
df$Name <- factor(df$Name, levels = df$Name)
ggplot(data = df, aes(x = Name, y = F1)) +
  geom_point(size = 3, color = "#0f7ab5") +
  geom_ribbon(aes(ymin = F1 - F1SE, ymax = F1 + F1SE, group = 1), fill = "#0f7ab5", alpha = 0.2) +
  geom_line(aes(group = 1), color = "#0f7ab5") +
  scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits = c(0.5, 1), oob = rescale_none) +
  custom_theme +
  xlab("")
dev.off()