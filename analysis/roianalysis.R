# This script is for looking at which brain regions are important for classification
### Plotting All Bar Graph Results ###
library(fabisearch)
library(rgl)

# Loop through directory and append results
all_results = list()
j = 1
files = list.files("../classifiers/ADNI/", ".rda")
files = files[!grepl("NCPD", files) & !grepl("CRMT", files)]
# Read in results and relabel
for(i in 1:length(files)){
  # Load the file
  load(paste0("../classifiers/ADNI/", files[i]))
  
  # Get a frequency table of the features
  CVfeatures = table(c(CVfeatures))
  all_results[[j]] = CVfeatures[order(CVfeatures, decreasing = TRUE)]
  
  # Iterate through k
  j = j + 1
}
# Add names
names(all_results) = sapply(strsplit(files, ".rda"), "[[", 1)

### Plotting the nodes ###
# Select the nodes which appear in all folds for cpDFC2
FBS_cpDFC2 = all_results$FBS_cpDFC2
allnames = names(FBS_cpDFC2[1:10])
split_names = strsplit(names(FBS_cpDFC2), "\\._")
split_metric_number = function(part) {
  match = regmatches(part, regexec("([a-zA-Z]+)(\\d+)", part))
  return(match[[1]][-1])
}
processed_data = lapply(split_names, function(x) c(x[1], split_metric_number(x[2])))
df = data.frame(Index = sapply(processed_data, function(x) x[1]),
                 Metric = sapply(processed_data, function(x) x[2]),
                 Number = sapply(processed_data, function(x) x[3]),
                 Value = FBS_cpDFC2, stringsAsFactors = FALSE)
df = df[,c(1:3, 5)]
FBS_cpDFC2 = sapply(strsplit(names(FBS_cpDFC2), "_"), "[[", 2)
FBS_cpDFC2 = as.numeric(unlist(regmatches(FBS_cpDFC2, gregexpr("[[:digit:]]+", FBS_cpDFC2))))[1:10]
FBS_cpDFC2 = unique(FBS_cpDFC2)
# Import the AAL atlas coordinates
AALatlas = read.table("../data/AALatlas/AAL120.txt")
coordROIs = AALatlas[,c(5,2:4)]
# Plot the nodes with net.3dplot
A = matrix(0, 120, 120)
net.3dplot(A, ROIs = FBS_cpDFC2, coordROIs = coordROIs, labels = TRUE)
rglwidget()

### Summarizing the features ###
load("../data/ADNI/FBS.rda")
filesrcs = list.files(path = "../functions", pattern="*.R$")
sapply(paste0("../functions/", filesrcs), source)
load("../data/ADNI/ADNIlist.rda")
# Estimate networks and graphs
FCdata = est.FC(ADNIlist, CPDlist, 2)
FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
FCdata = cbind(subj = FCdata$subj, class = FCdata$class,
               graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
# Convert dataframe to be side by side
FCdata = concat.windows(FCdata)
FCdata = FCdata[,colnames(FCdata) %in% c("class", allnames)]
colMeans(FCdata[FCdata$class == 1,])
colMeans(FCdata[FCdata$class == 0,])