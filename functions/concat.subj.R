# This function is used to concatenate edge weights or graph stats for each subject
# such that it will create one line for each subject
library(parallel)

concat.subj = function(data){
  # data = input classifier data, with subject column
  
  # Split the data into a list where each element is a subject
  data = split(data, data$subj)
  
  # Apply the concat.subj.each function over list elements
  data = mclapply(data, concat.subj.each)
  data = do.call("rbind", data)
}

# This is a helper function over which concat.subj parallelizes over
concat.subj.each = function(subjdata){
  # subjdata = one particular subjects data
  
  # Save subject and class information separately, drop from the other
  subjclass = subjdata[1, colnames(subjdata) %in% c("subj", "class")]
  indvars = subjdata[, !colnames(subjdata) %in% c("subj", "class")]
  
  # Save column names
  indvarsnames = colnames(indvars)
  indvarsnames = paste0(rep(indvarsnames, nrow(indvars)), ".", 
                        sort(rep(1:nrow(indvars), length(indvarsnames))))
  
  # Create a vector of the values
  indvars = as.vector(t(as.matrix(indvars)))
  indvars = matrix(indvars, nrow = 1)
  colnames(indvars) = indvarsnames
  indvars = as.data.frame(indvars)
  
  # Bind to the subject and class data
  output = cbind(subjclass, indvars)
  
  # Return the final row
  return(output)
}
