# This script is for the analysis of features across groups
# Load libraries
library(ggplot2)

### 1 Feature Distances ###
# Load data
load("alldist.rda")

# Combine for plotting
df = append(list(cpDFCdist), allwDFCdist)
names(df)[1] = "cpDFC"
# Transform into dataframe
df = as.data.frame(do.call(cbind, df))
# Then stack
df = stack(df)

# Plot
png("Boxplots.png", width = 20, height = 10, units = "in", res = 1200, pointsize = 4)
ggplot(df, aes(x = ind, y = values)) + geom_boxplot() + xlab("DFC Condition") + 
  ylab("Mean Euclidean Distance Across States") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### 2 Feature Stability ###
# Load data
# Define the directories to go through
dirs = list.dirs("../classifiers")

# Loop through directory and append results
all_results = c()
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
    nme = paste0(strsplit(curdir, "_")[[1]][2],
                 gsub("_LOOCVresults.rda", "", gsub("0.5", "", curfiles[j])))
    if(grepl("win", nme)){
      grp = "wDFC"
    } else if(grepl("cpd", nme)){
      grp = "cpDFC"
    } else {
      grp = "SFC"
    }
  }
  
  # Get all features
  allfeat = sapply(strsplit(as.vector(CVfeatures), "_"), "[[", 2)
  
  # Save the number of unique features
  allfeat = length(unique(allfeat))
}

# Length/number of features
length(unique(as.vector(CVfeatures)))
