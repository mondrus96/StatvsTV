# This file contains helper functions used to parallelize the est.FC main function
# Functions are organized to increasing levels of abstraction further down, ie.
# est.block.FC < est.scan.FC < est.subj.FC

# This is a helper function used to estimate individual blocks of FC
est.block.FC = function(blockdata){
  # blockdata = the input block of data
  
  # Estimate the functional connectivity using correlation
  edgemat = cor(blockdata)
  
  # Define the FC matrix
  FC = abs(edgemat[upper.tri(edgemat,diag=F)])
  
  # Define the flattened column 
  rowcol = expand.grid(rownames(edgemat), colnames(edgemat))
  labs = rowcol[as.vector(upper.tri(edgemat,diag=F)),]
  labs = paste0("_", labs[,1], ":", labs[,2])
  
  # Add the names to the FC vector
  names(FC) = labs
  
  # Return the FC estimate
  return(FC)
}

# This is a helper function used to estimate FC for each scan
est.scan.FC = function(scandata, cpts = NULL, type = "cpts"){
  # scandata = the subjects scan data
  # cpts = change points for the scan, if NULL, assumed stationary estimation
  # type = type of analysis to run - if cpts, then changepoints are used by default,
  #        otherwise, if "window" supplied, then window based method will be used
  
  if(type == "cpts"){
    # Separate based on type of FC estimation
    if(is.null(cpts)){
      # Estimate FC for the scan with no change points
      scanFC = est.block.FC(scandata)
      
      # Turn into a flat dataframe
      scanFC = t(data.frame(scanFC))
    } else if(!is.null(cpts)){
      # Add end points to the change point data
      cpts = c(0, cpts, nrow(scandata))
      
      # Separate out the scan into individual blocks based on change points
      scandatalist = list()
      for(i in 1:(length(cpts)-1)){
        # Define the beginning and end
        Tstart = cpts[i] + 1
        Tend = cpts[i+1]
        
        # Select the block of data and append to output
        scandatalist[[i]] = scandata[Tstart:Tend,]
      }
      # Apply the est.block.FC function
      scanFC = mapply(est.block.FC, scandatalist, SIMPLIFY = FALSE)
      scanFC = do.call(rbind, scanFC)
    }
  } else if(type == "window"){
    # Separate out the scan into individual blocks based on windows
    scandatalist = list()
    Tstart = seq(1, (nrow(scandata) - cpts[1]), cpts[2])
    for(i in 1:length(Tstart)){
      # Select the block of data and append to output
      scandatalist[[i]] = scandata[Tstart[i]:(Tstart[i] + cpts[1]),]
    }
    # Apply the est.block.FC function
    scanFC = mapply(est.block.FC, scandatalist, SIMPLIFY = FALSE)
    scanFC = do.call(rbind, scanFC)
  }
  # Flatten and return the output
  return(scanFC)
}

# This is a helper function used to estimate FC for each subject
est.subj.FC = function(subjdata, cpd, params = NULL, window = FALSE){
  # subjdata = the subjects data
  # cpd = changepoint detection data
  # params = cutoff value used for change point detection, if relevant
  # window = whether to use window approach or not
  
  # If params is not NULL, then estimate change points otherwise assign NULL
  if(window == FALSE){
    if(!is.null(params)){
      # Loop through and retrieve change points
      if(params < 1){
        # Use cutoff value
        cpts = cpd$change_points$T[cpd$change_points$stat_test < params]
      } else if(params >= 1){
        # Keep only the first param number of changepoints
        cpts = cpd$change_points$T[order(cpd$change_points$stat_test) <= params]
      }
    } else if(is.null(params)){
      # Define cpts as NULL
      cpts = NULL
    }
    # Estimate the subjects FC
    subjFC = est.scan.FC(subjdata$fmri, cpts)
  } else if(window == TRUE){
    # Pass along params into cpts
    cpts = params
    
    # Estimate the subjects FC
    subjFC = est.scan.FC(subjdata$fmri, cpts, type = "window")
  }
  
  # Return the output
  return(subjFC)
}
