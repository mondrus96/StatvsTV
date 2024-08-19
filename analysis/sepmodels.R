# Load libraries
library(ggplot2)
library(plot3D)
library(reshape2)
library(scales)
library(caret)

### Plotting All Bar Graph Results ###
# Define the directories to go through
dirs = list.dirs("../classifiers")
models = "linsvm"

# Loop through directory and append results
all_preds = cols_all_preds = vector("list", length(models))
all_results = c()
names(all_preds) = names(cols_all_preds) = models
k = 1
for(i in 2:length(dirs)){
  # Current directory
  curdir = dirs[i]
  
  # List the files
  curfiles = list.files(curdir)
  curfiles = curfiles[grepl(".rda", curfiles)]
  
  # Read in results and relabel
  for(j in 1:length(curfiles)){
    # Load the file
    load(paste0(curdir, "/", curfiles[j]))
    
    # Create a name and also split into groups
    class = strsplit(curdir, "/")[[1]][3]
    nme = strsplit(curfiles[j], ".rda")[[1]][1]
    if(grepl("wDFC", nme)){
      grp = "wDFC"
    } else if(grepl("cpDFC", nme)){
      grp = "cpDFC"
    } else if(grepl("SFC", nme)){
      grp = "SFC"
    }
    
    # Find the number of unique features
    uniqfeat = sapply(strsplit(as.vector(CVfeatures), "_"), "[[", 2)
    uniqfeat = length(unique(uniqfeat))
  
    # Append this to the current CVeval dataframe
    CVeval = cbind(class, grp, nme, CVeval, uniqfeat)
    
    # Append this to the final dataframes
    all_results = rbind(all_results, CVeval)
    for(i in 1:length(models)){
      all_preds[[i]] = cbind(all_preds[[i]], CVpreds[,colnames(CVpreds) == models[i]])
      colnames(all_preds[[i]])[k] = paste0(class, ".", nme) 
    }
  k = k + 1
  }
}

# All the wDFC results
wDFC = all_results$Accuracy[all_results$grp == "wDFC" & all_results$class == "ADNI" & all_results$Model == "linsvm"]
wDFC = as.numeric(wDFC)
df = data.frame(wDFC)
base_size = 15
png("wDFC.png", width = 8, height = 5, units = "in", res = 120, pointsize = 10)
ggplot(df, aes(x = wDFC)) +
  geom_histogram(binwidth = 0.05, fill = "#0f7ab5", color = "black") +
  labs(x = "Accuracy", y = "Count", title = "Distribution of wDFC Accuracies (ADNI)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = base_size * 1.5),
    axis.title.x = element_text(size = base_size * 1.5),
    axis.title.y = element_text(size = base_size * 1.5),
    axis.text.x = element_text(size = base_size * 1.5),
    axis.text.y = element_text(size = base_size * 1.5)
  )
dev.off()

# Convert all_results columns to numeric
all_results[,5:ncol(all_results)] = apply(all_results[,5:ncol(all_results)], 2, as.numeric)

# Select only linsvm
all_results = all_results[all_results$Model == "linsvm",]

# Calculate SE for each of the quantities
all_results$AccSE = sqrt((all_results$Accuracy*(1-all_results$Accuracy))/68)
all_results$SenSE = sqrt((all_results$Sensitivity*(1-all_results$Sensitivity))/68)
all_results$SpeSE = sqrt((all_results$Specificity*(1-all_results$Specificity))/68)
all_results$F1SE = sqrt((all_results$F1*(1-all_results$F1))/68)

# Define a custom theme with larger text sizes
custom_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16),  # Increase legend title size
    plot.title = element_text(size = 18),  # Increase plot title size
    panel.background = element_blank()
  )

# Specify colours
cols <- c("#c6475b", "#67A968", "#0f7ab5")

# Plot Accuracy
df <- all_results[all_results$class == "ADNI" & !all_results$nme %in% c("CRMT_cpDFC1", "CRMT_cpDFC2", "NCPD_cpDFC1", "NCPD_cpDFC2"),]
df <- df[df$Accuracy > 0.6359212 | df$nme == "SFC",]

# Sort dataframe by F1 score
df <- df[order(df$F1, decreasing = TRUE), ]

# Plot Accuracy
png("Accuracy.png", width = 8, height = 5, units = "in", res = 120, pointsize = 10)
ggplot(data = df, aes(x = reorder(nme, -F1), y = Accuracy, fill = grp)) + 
  labs(fill = "Group") +
  geom_bar(stat = "identity") + 
  scale_y_continuous(breaks = seq(0.3, 1, 0.1), limits = c(0.3, 1), oob = rescale_none) +
  scale_fill_manual(values = cols) +
  custom_theme + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") +
  geom_errorbar(aes(ymin = Accuracy - AccSE, ymax = Accuracy + AccSE), width = 0.2, position = position_dodge(0.9))
dev.off()

# Plot Sensitivity
png("Sensitivity.png", width = 8, height = 5, units = "in", res = 120, pointsize = 10)
ggplot(data = df, aes(x = reorder(nme, -F1), y = Sensitivity, fill = grp)) + 
  labs(fill = "Group") +
  geom_bar(stat = "identity") + 
  scale_y_continuous(breaks = seq(0.3, 1, 0.1), limits = c(0.3, 1), oob = rescale_none) +
  scale_fill_manual(values = cols) +
  custom_theme + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") +
  geom_errorbar(aes(ymin = Sensitivity - SenSE, ymax = Sensitivity + SenSE), width = 0.2, position = position_dodge(0.9))
dev.off()

# Plot Specificity
png("Specificity.png", width = 8, height = 5, units = "in", res = 120, pointsize = 10)
ggplot(data = df, aes(x = reorder(nme, -F1), y = Specificity, fill = grp)) + 
  labs(fill = "Group") +
  geom_bar(stat = "identity") + 
  scale_y_continuous(breaks = seq(0.3, 1, 0.1), limits = c(0.3, 1), oob = rescale_none) +
  scale_fill_manual(values = cols) +
  custom_theme + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") +
  geom_errorbar(aes(ymin = Specificity - SpeSE, ymax = Specificity + SpeSE), width = 0.2, position = position_dodge(0.9))
dev.off()

# Plot F1 Score
png("F1.png", width = 8, height = 5, units = "in", res = 120, pointsize = 10)
ggplot(data = df, aes(x = reorder(nme, -F1), y = F1, fill = grp)) + 
  labs(fill = "Group") +
  geom_bar(stat = "identity") + 
  scale_y_continuous(breaks = seq(0.3, 1, 0.1), limits = c(0.3, 1), oob = rescale_none) +
  scale_fill_manual(values = cols) +
  custom_theme + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") +
  geom_errorbar(aes(ymin = F1 - F1SE, ymax = F1 + F1SE), width = 0.2, position = position_dodge(0.9))
dev.off()

# Plot Unique Features
png("Unique.png", width = 8, height = 5, units = "in", res = 120, pointsize = 10)
ggplot(data = df, aes(x = reorder(nme, -F1), y = uniqfeat, fill = grp)) + 
  labs(fill = "Group") +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cols) +
  custom_theme + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") + 
  ylab("Number of Unique Features")
dev.off()

### Plotting Heatmaps for wDFC ###
# Prepare data for plotting
df = all_results[all_results$grp == "wDFC",]
df = df[df$class == "ADNI",]
# Create vectors for window and step size, which will be the x and y
nmes = strsplit(gsub("wDFC", "", df$nme), "_")
nmes = do.call(rbind, nmes)
x = (gsub("win", "", nmes[,1]))
y = (gsub("step", "", nmes[,2]))

# Loop through different values and plot heatmaps
vals = c("Accuracy", "Sensitivity", "Specificity", "F1")
nms = c("Accuracy", "Sensitivity", "Specificity", "F1 Score")
for(i in 1:length(vals)){
  heatdf = as.data.frame(cbind(x, y, 100*df[,colnames(df) == vals[i]]))
  heatdf$V3 = as.numeric(heatdf$V3)
  png(paste0(vals[i], "_win.png"), width = 600, height = 300, res = 80)
  ggplot(heatdf, aes(x, y = reorder(y, as.numeric(y)), fill=V3)) + 
    geom_tile() + custom_theme + 
    xlab("Window Size") + ylab("Step Size") + labs(fill = nms[i]) +
    geom_text(aes(label = round(V3, 1))) + 
    scale_fill_gradientn(colors = c("red", "green"))
  dev.off()
}

### Statistical Tests ###
# Comparing 8ens to cpDFC2
mod1 = all_results$Accuracy[all_results$class == "ensADNImain" & all_results$nme == "8ens"]
mod2 = all_results$Accuracy[all_results$class == "ADNI" & all_results$nme == "FBS_cpDFC2"]
print(paste0("8ens > cpDFC2: ", prop.test(mod1*68, 68, mod2, "greater")$p.value)) # reject null

# Comparing wDFC60win_15step to cpDFC2
mod1 = all_results$Accuracy[all_results$class == "ADNI" & all_results$nme == "FBS_cpDFC2"]
mod2 = all_results$Accuracy[all_results$class == "ADNI" & all_results$nme == "wDFC60win_15step"]
print(paste0("cpDFC2 > wDFC60win_15_step: ", prop.test(mod1*68, 68, mod2, "greater")$p.value)) # fail to reject null

# Comparing cpDFC2 to wDFC condition
vals = c("Accuracy", "Sensitivity", "Specificity", "F1")
for(i in 1:length(vals)){
  cpDFC2 = all_results[all_results$class == "ADNI"
                       & all_results$nme == "FBS_cpDFC2",
                       colnames(all_results) == vals[i]]
  wDFC = all_results[all_results$class == "ADNI"
                     & all_results$grp == "wDFC",
                     colnames(all_results) == vals[i]]
  muP = mean(wDFC)
  sigP = ((muP*(1-muP))/length(wDFC))^0.5
  print(paste0(vals[i], ": ", pnorm(cpDFC2, muP, sigP, lower.tail = FALSE))) # reject all nulls
}
