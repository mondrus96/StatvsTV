# This function will estimate functional connectivity (static or dynamic) and return a dataframe which contains
library(parallel)
library(doParallel)

# The necessary items for training/testing the models
est.FC = function(data, cpds = NULL, params = NULL){
  # data = the input data
  # cpds = change point detection data, if relevant
  # params = cutoff value used for change point detection, if relevant
  # cutoff = cutoff used to estimate networks for FC
  # window = whether to use window approach or not
  
  # Define tabdata, which will be returned by the function
  tabdata = c()
  
  # If params is a vector with two components then window = TRUE
  if(length(params) == 2){
    window = TRUE
  } else {
    window = FALSE
  }
  
  # Make cpds length of data if null
  if(is.null(cpds)){
    cpds = rep(NA, length(data))
  }
  
  # Use mclapply to parallelize correlations
  tabdata = mcmapply(est.subj.FC, data, cpds, SIMPLIFY = FALSE,
                     MoreArgs = list(params = params, window = window))
  
  # Add subject numbers to each component in the list
  tabdata = mapply(cbind, subj = 1:length(tabdata), 
                   class = sapply(data, "[[", 1),
                   tabdata, SIMPLIFY = FALSE)
  
  # Flatten the results
  tabdata = do.call(rbind, tabdata)
  
  # Remove the row names from the tabular data
  row.names(tabdata) = NULL
  tabdata = as.data.frame(tabdata)
  return(tabdata)
}
