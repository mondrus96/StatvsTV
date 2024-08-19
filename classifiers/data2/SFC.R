### For All Experiments ###
# Inputs
load("workspace.RData")

### SFC ###
# Estimate networks and graphs
FCdata = est.FC(datalist)
FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
FCdata = cbind(subj = FCdata$subj, class = FCdata$class,
               graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
# Run classification
outname = "SFC"
testing.classifiers(FCdata, outname, dimred, dimredparams)