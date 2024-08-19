# This function returns a list of the hold out samples for 10 fold CV when the number of samples is not divisible by 10
# Each element in the list is a vetor of indices to be used as the holdout sample

smpl.size = function(data.size){
  # Find the sample sizes that work with the data.size
  low = floor(data.size/10)
  high = ceiling(data.size/10)
  prop = as.integer(as.character((data.size/10 - low)*10))
  samplesize = c(0, rep.int(low, 10-prop), rep.int(high, prop))
  
  # Turn this into a vector of lists
  samplelist = list()
  for(i in 1:10){
    samplelist[[i]] = (sum(samplesize[1:i])+1):(sum(samplesize[1:i+1]))
  }
  
  return(samplelist)
}