# Load libraries
library(caret)
library(doParallel)
library(parallel)
library(SIS)

### This function is used to run the main classifier protocol ###
run.classifiers = function(datalist, outname, regvals, nrun = NULL, alpha = NULL, CVtype = NULL, 
                           cutoff = NULL, models = NULL, dimred = NULL, dimredparams = NULL, 
                           netstats = FALSE){
  # datalist = list of data to use as an input into the classifier models
  # outname = name of the output
  # regvals = regularization parameters for classifiers that require it
  # nrun = nrun used by NMF, if non NULL NMF is used, correlation is default
  # alpha = cutoff value for change point detection, if NULL then assumed static networks
  # CVtype = type of cross validation to use
  # cutoff = cutoff used to estimate networks
  # models = models to try, by default is set to NULL and tries all of them
  # dimred = which dimension reduction to use
  # dimredparams = dimension reduction parameters/inputs to use for the selected method
  # netstats = whether to estimate network statistics (TRUE) or just use plain network edges (FALSE)
  
  # Set seed
  set.seed(123)
  
  # Register the parallel backend
  numcores = detectCores()
  registerDoParallel(cl = numcores)
  print(paste("Number of Cores:", numcores))
  
  # For 10 fold CV with a non-divisible number of subjects, split into two different hold out sizes
  testsmpls = smpl.size(length(datalist))
  
  # Print the functional connectivity estimation procedure
  if(is.null(alpha)){
    print("Estimating Static Networks")
  } else if(!is.null(alpha)){
    print("Estimating CPD Networks")
  }
  
  # Estimate functional connectivity for all subjects
  FCdata = est.FC(datalist, alpha, nrun, cutoff)
  
  # Keep only the first scan and remove unneeded columns
  FCdata = FCdata[FCdata$scan == 1, !colnames(FCdata) %in% c("scan", "rank", "Tstart", "Tend")]
  
  # If relevant, estimation graph summary statistics
  if(netstats == TRUE){
    # Select only the columns containing edge information
    graphstats = graph.stats(FCdata[, grepl("_", colnames(FCdata))])
    
    # Swap with the raw edge weights from FCdata
    FCdata = cbind(FCdata[, !grepl("_", colnames(FCdata))], graphstats)
  }
  
  # If alpha > 1, concatenate variables from each subject into one line
  if(!is.null(alpha) && alpha >= 1){
    FCdata = concat.subj(FCdata)
  }
  
  # Perform CV
  if (CVtype == "LOO"){
    folds = unique(FCdata$subj)
  } else if (CVtype == "10fold"){
    folds = 1:10
  }
  CVpreds = c()
  CVeval = c()
  hyperall = list()
  for(i in folds){
    # Set the seed based on the index
    set.seed(i*12345)
    
    # Define the randomly reshuffled indices
    randind = sample(1:length(datalist))
    
    # Print progress
    print(paste0("Evaluting CV Fold: ", i))
    
    # Split into train and test sets, keep only the first scan
    if(CVtype == "LOO"){
      trainsubj = (1:95)[-i]
      testsubj = i
    } else if (CVtype == "10fold"){
      trainsubj = randind[!1:length(datalist) %in% testsmpls[[i]]]
      testsubj = randind[testsmpls[[i]]]
    }
    
    ### Use subset selection ###
    # Print the dimension reduction technique being used
    print(paste0("Performing Variable Screening"))
    
    # Define the data to be subselected for
    data2subset = FCdata[FCdata$subj %in% trainsubj, grepl("_", colnames(FCdata)) | 
                              colnames(FCdata) == "class"]
    
    # Apply dimension reduciton technique
    if (is.null(dimred)){
      train = FCdata[FCdata$subj %in% trainsubj,]
      test = FCdata[FCdata$subj %in% testsubj,] 
    } else if (dimred == "SIS"){
      # Use the SIS function to screen
      finalvars = SIS(as.matrix(data2subset[,grepl("_", colnames(data2subset))]), 
                      data2subset$class, family = "binomial", varISIS = "cons",
                      standardize = TRUE, iter = dimredparams$iter)
      
      if(dimredparams$iter == TRUE){
        finalvars = finalvars$ix
      } else if (dimredparams$iter == FALSE){
        finalvars = finalvars$sis.ix0
      }
      finalvars = colnames(data2subset[,grepl("_", colnames(data2subset))])[finalvars]
      print(paste0("Number of Variables: ", length(finalvars)))
      
      # If the number of variables selected is 0, then go to next CV fold
      if(length(finalvars) == 0){
        next
      }
      
      # Apply this to the train and test subject data
      train = FCdata[FCdata$subj %in% trainsubj, !grepl("_", colnames(FCdata)) | 
                              colnames(FCdata) %in% c(finalvars)]
      test = FCdata[FCdata$subj %in% testsubj, !grepl("_", colnames(FCdata)) | 
                            colnames(FCdata) %in% c(finalvars)]
    } else if (dimred == "PCA"){
      # Fit PCA
      pcaout = prcomp(as.data.frame(data2subset[,grepl("_", colnames(data2subset))]), 
                                    center=TRUE, scale=TRUE)
      
      # Use PCA on train and testdata
      pcadata = predict(pcaout, as.data.frame(FCdata[,grepl("_", colnames(FCdata))]))
      colnames(pcadata) = paste0("_", colnames(pcadata))
      
      # Pick the appropriate cutoff for variance explained
      varexpl = cumsum(pcaout$sdev^2/sum(pcaout$sdev^2)) <= dimredparams$varexp
      pcadata = cbind(FCdata[,!grepl("_", colnames(FCdata))], pcadata[,varexpl])
      
      # Append to other information
      train = pcadata[pcadata$subj %in% trainsubj,]
      test = pcadata[pcadata$subj %in% testsubj,]
    }
    
    ### Tune hyperparameters ###
    hyperparams = tune.hyperparams(train, regvals, models)
    
    # Run the models and save the predicted and actual classes
    CVpreds[[i]] = train.test.CV(train, test, hyperparams, models)
    
    if (CVtype == "10fold"){
      # Evaluate the models for that fold
      CVeval = rbind(CVeval, eval.CV(CVpreds[[i]]))
      hyperall[[i]] = hyperparams
      
      # Save results
      FCtype = ifelse(is.null(alpha), "static", alpha)
      save(CVeval, hyperall, file = paste0(outname, "_", FCtype, "_", "10foldCVresults.rda"))
    } else if (CVtype == "LOO"){
      # Evaluate the models across all LOOCV
      modnames = names(CVpreds[[1]])
      CVeval = list()
      for(j in 1:length(modnames)){
        curreval = do.call(rbind, sapply(CVpreds, `[`, j))
        rownames(curreval) = NULL
        CVeval[[modnames[j]]] = curreval
      }
      
      # Run evaluation
      CVeval = eval.CV(CVeval)
      
      # Save results
      FCtype = ifelse(is.null(alpha), "static", alpha)
      save(CVeval, hyperall, file = paste0(outname, "_", FCtype, "_", "LOOCVresults.rda"))
    }
  }
 }