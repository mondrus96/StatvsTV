### For All Experiments ###
# Inputs
dimred = "SIS"
dimredparams = list(iter = FALSE, nsis = 10)
# Load all the functions
filesrcs = list.files(path = "../../functions", pattern="*.R$")
sapply(paste0("../../functions/", filesrcs), source)
# Load the data
load("../../data/ADNI/ADNIlist.rda")
save.image(file = "workspace.RData")