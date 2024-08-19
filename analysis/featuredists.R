# Script is for comparing features
# Define a function for extracting the distance between adjacent time windows
adjdist = function(dists){
  # Convert to matrix and get off-diagonals
  dists = as.matrix(dists)
  dists = mean(dists[row(dists) == (col(dists) - 1)])
  
  # Return output
  return(dists)
}

# Load all the functions
filesrcs = list.files(path = "../functions", pattern="*.R$")
sapply(paste0("../functions/", filesrcs), source)
# Load the data
load("../data/ADNI/ADNIlist.rda")

### 1 - Changepoint Detection with two CPs ###
# Estimate networks and graphs
FCdata = est.FC(ADNIlist, 2)
FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
FCdata = cbind(FCdata[,!grepl("_", colnames(FCdata))],
               graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
# Separate by class
cpDFCdist = split(FCdata, FCdata$subj)
# Find the average distance across adjacent time windows
cpDFCdist = lapply(cpDFCdist, dist)
cpDFCdist = sapply(cpDFCdist, adjdist)

### 2 - Window based DFC ###
# Loop through all possible combinations
wsize = seq(10, 70, 5)
wstep = c(1, 2, 3, 5, 8, 10, 15, 20)
grid = expand.grid(wsize, wstep)
allwDFCdist = vector(mode = "list", nrow(grid))
# Loop grid variable
for(i in 1:nrow(grid)){
  # Print current step
  print(paste0("window", grid[i,1], "step", grid[i,2]))
  
  # Estimate networks/graphs using windows
  FCdata = est.FC(ADNIlist, params = unlist(grid[i,]))
  FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
  FCdata = cbind(FCdata[,!grepl("_", colnames(FCdata))],
                 graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
  # Separate by subject
  wDFCdist = split(FCdata[,-(1:3)], FCdata$subj)
  # Find the average distance across adjacent time windows
  wDFCdist = lapply(wDFCdist, dist)
  wDFCdist = sapply(wDFCdist, adjdist)
  # Append to dataframe
  allwDFCdist[[i]] = wDFCdist
  names(allwDFCdist)[i] = paste0("win", grid[i,1], "step", grid[i,2])
}

# Save as table(s)
save(cpDFCdist, allwDFCdist, file = "alldist.rda")