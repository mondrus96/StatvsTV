### For All Experiments ###
# Inputs
load("../../data/data2/CRMT.rda")
load("workspace.RData")

### 1 Change Point ###
# Estimate networks and graphs
FCdata = est.FC(datalist, CPDlist, 1)
FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
FCdata = cbind(subj = FCdata$subj, class = FCdata$class,
               graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
# Convert dataframe to be side by side
FCdata = concat.windows(FCdata)
# Run classification
outname = "CRMT_cpDFC1"
testing.classifiers(FCdata, outname, dimred, dimredparams)


### 2 Change Points ###
# Estimate networks and graphs
FCdata = est.FC(datalist, CPDlist, 2)
FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
FCdata = cbind(subj = FCdata$subj, class = FCdata$class,
               graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
# Convert dataframe to be side by side
FCdata = concat.windows(FCdata)
# Run classification
outname = "CRMT_cpDFC2"
testing.classifiers(FCdata, outname, dimred, dimredparams)