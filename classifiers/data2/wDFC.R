### For All Experiments ###
# Inputs
load("workspace.RData")

### Window based DFC ###
# Define window and step sizes
wsize = seq(10, 70, 5)
wstep = c(1, 2, 3, 5, 8, 10, 15, 20)
# Loop through window sizes
for(i in wsize){
  # Loop across steps
  for(j in wstep){
    # Print current step
    print(paste0("window", i, "step", j))
    
    # Estimate networks/graphs using windows
    FCdata = est.FC(datalist, params = c(i, j))
    FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
    FCdata = cbind(FCdata[,!grepl("_", colnames(FCdata))],
                   graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
    
    # Convert dataframe to be side by side
    FCdata = concat.windows(FCdata)
    
    # Run classification
    outname = paste0("wDFC", i, "win_", j, "step")
    testing.classifiers(FCdata, outname, dimred, dimredparams)
  }
}
