# This script will concatenate different windows into one line per subject
concat.windows = function(data, perm = FALSE){
  # data = input data to concatenate over
  # perm = whether to randomly permute over time
  
  # Set seed for replicable results
  set.seed(123)
  
  # Find the number of rows per subject
  rowvec = 1:(sum(data$subj == data$subj[1]))
  
  # Perm false or true
  if(perm == FALSE){
    # Concatenate columns together
    vars = split(data[,grepl("_", colnames(data))], rowvec)
    vars = do.call(cbind, vars)
    outdata = cbind(split(data[,!grepl("_", colnames(data))], rowvec)[[1]], vars)
  } else if(perm == TRUE){
    # Loop through all and randomly reshuffle
    outdata = c()
    for(i in 1:length(unique(data$subj))){
      # Select current subject and permute rows
      cursubj = data[data$subj == i,]
      cursubj = cursubj[order(runif(length(rowvec))),]
      
      # Concatenate columns together
      vars = split(cursubj[,grepl("_", colnames(cursubj))], rowvec)
      vars = do.call(cbind, vars)
      vars = cbind(cursubj[,!grepl("_", colnames(cursubj))][1,], vars)
      outdata = rbind(outdata, vars)
    }
  }
  
  # Return data
  return(outdata)
}
