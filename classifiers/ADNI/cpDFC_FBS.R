### For All Experiments ###
# Inputs
load("../../data/ADNI/FBS.rda")
load("workspace.RData")

### 1 Change Point (FaBiSearch) ###
# Estimate networks and graphs
FCdata = est.FC(ADNIlist, CPDlist, 1)
FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
FCdata = cbind(subj = FCdata$subj, class = FCdata$class,
               graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
# Convert dataframe to be side by side
FCdata = concat.windows(FCdata)
# Run classification
outname = "FBS_cpDFC1"
testing.classifiers(FCdata, outname, dimred, dimredparams)


### 2 Change Points (FaBiSearch) ###
# Estimate networks and graphs
FCdata = est.FC(ADNIlist, CPDlist, 2)
FCdata[,grepl("_", colnames(FCdata))] = 1*(abs(FCdata[,grepl("_", colnames(FCdata))]) > 0.5)
FCdata = cbind(subj = FCdata$subj, class = FCdata$class,
               graph.stats(FCdata[,grepl("_", colnames(FCdata))], "all"))
# Convert dataframe to be side by side
FCdata = concat.windows(FCdata)
# Run classification
outname = "FBS_cpDFC2"
testing.classifiers(FCdata, outname, dimred, dimredparams)