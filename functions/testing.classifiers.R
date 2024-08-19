# This function will run the classifier method
library(parallel)

testing.classifiers = function(FCdata, outname, dimred, dimredparams){
  # FCdata = data to run
  # outname = output name
  # dimred = dimension reduction method
  # dimredparams = list of relevant dimension reduction parameters
  
  # Define data length
  datalen = max(unique(FCdata$subj))

  # Register the parallel backend
  numcores = detectCores()
  registerDoParallel(cl = numcores)
  print(paste("Number of Cores:", numcores))
  
  # Perform CV
  folds = unique(FCdata$subj)
  CVpreds = c()
  CVfitted = vector("list", length(folds))
  CVeval = c()
  CVfeatures = c()
  hyperall = list()
  for(i in folds){
    # Set the seed based on the index
    set.seed(i*12345)
    
    # Print progress
    print(paste0("Evaluting CV Fold: ", i))
    
    # Split into train and test sets, keep only the first scan
    trainsubj = (1:length(folds))[-i]
    testsubj = i
    
    ### Use subset selection ###
    # Print the dimension reduction technique being used
    print(paste0("Performing Variable Screening"))
    
    # Define the data to be subselected for
    data2subset = FCdata[FCdata$subj %in% trainsubj, (grepl("_", colnames(FCdata)) | 
                           colnames(FCdata) == "class")]
    
    if(dimred == "SIS"){
      # Use the SIS function to screen
      finalvars = SIS(as.matrix(data2subset[, grepl("_", colnames(data2subset))]),
                      data2subset$class, family = "binomial", iter = dimredparams$iter,
                      standardize = TRUE, nsis = dimredparams$nsis)
      finalvars = finalvars$ix0 
      keepcols = data2subset[, grepl("_", colnames(data2subset))]
      CVfeatures = rbind(CVfeatures, colnames(keepcols[, finalvars]))
      keepcols = c("subj", "class", colnames(keepcols[, finalvars]))
      
      # Apply this to the train and test subject data
      train = as.data.frame(FCdata[trainsubj, colnames(FCdata) %in% keepcols, drop = FALSE])
      test = as.data.frame(FCdata[testsubj, colnames(FCdata) %in% keepcols, drop = FALSE])
    } else if(dimred == "none"){
      # Define the train and test dataframes as they are without variable selection
      train = as.data.frame(FCdata[trainsubj, ])
      test = as.data.frame(FCdata[testsubj, ])
    }
    # Run the models
    splitFCdata = list(train, test)
    names(splitFCdata) = c("train", "test")
    outCV = train.test.CV(splitFCdata)
    
    # Run the models and save the predicted and actual classes
    CVpreds = rbind(CVpreds, outCV$preds)
    
    # Also save the fitted values at each iteration
    CVfitted[[i]] = outCV$fitted
    
    # Run evaluation
    CVeval = eval.CV(CVpreds)
    print(CVeval)
    
    # Save results
    save(CVeval, CVfeatures, CVpreds, CVfitted, file = paste0(outname, ".rda"))
  }
}
