# Load libraries
library(SIS)

# This function is for building a meta-learning which takes in multiple outputs
mean.combine.preds = function(dir, modacc, preds){
  # dir = directory name which contains results to pull from
  # modacc = model accuracies
  # preds = which predictions to combine over
  
  # Find which predictors to use
  prednames = c()
  for(i in 2:ncol(modacc)){
    prednames = cbind(prednames, modacc[order(modacc[,i], decreasing = TRUE),1][preds])
    colnames(prednames)[i-1] = colnames(modacc)[i]
  }
  
  # Loop through prediction names and get results
  filesrcs = list.files(path = paste0("../", dir), pattern="*.rda$")
  allpreds = vector("list", ncol(prednames))
  names(allpreds) = colnames(prednames)
  for(i in 1:ncol(prednames)){
    for(j in 1:nrow(prednames)){
      # Load relevant file
      load(paste0("../", dir, "/", filesrcs[grepl(prednames[j,i], filesrcs)]))
      
      # Append to predictions
      allpreds[[i]] = cbind(allpreds[[i]], CVpreds[,colnames(CVpreds) == names(allpreds)[i]])  
    }
    allpreds[[i]] = 1*(rowMeans(data.frame(allpreds[[i]])) > 0.5)
  }
  allpreds = do.call(cbind, allpreds)
  allpreds = cbind(CVpreds$class, allpreds)
  colnames(allpreds)[1] = "class"
  
  # Run evaluation
  CVeval = eval.CV(allpreds)
  print(CVeval)
  
  # Save results
  save(CVeval, file = paste0(length(preds), "ens.rda"))
}

# This function will use SIS to build the best meta classifier model
SIS.combine.preds = function(all_fitted, all_preds, nvar){
  # all_fitted = list of all fiited values, each element is a different FC condition
  # all_preds = list of all predicted values, each element is a different FC condition
  # nvars = number of variables to recruit from SIS
  
  # Save cpDFC2
  cpDFC2 = all_fitted[names(all_fitted) == "cpDFC2"]
  nobs = length(cpDFC2[[1]])
  
  # Loop through all the models
  models = names(all_fitted[[1]][[1]])
  models = models[models != "class"]
  CVpreds = CVeval = vector("list", 2*length(models))
  for(i in 1:length(models)){
    # Loop through all CV folds
    for(j in 1:nobs){
      # Select the current models x values
      x_train = all_fitted[!names(all_fitted) == "cpDFC2"]
      x_train = do.call(cbind, lapply(lapply(x_train, "[[", j), "[[", models[i]))
      
      # Define the y_train vector
      y_train = do.call(cbind, lapply(lapply(cpDFC2, "[[", j), "[[", 1))
      
      # Run SIS
      sisout = SIS(x_train, y_train, family = "binomial", nsis = nvar-1, iter = TRUE, tune = "cv", type.measure = "auc", penalty = "lasso")
      finalvars = c(which(names(all_fitted) == "cpDFC2"), sisout$ix0)
      x_train = do.call(cbind, lapply(lapply(all_fitted[finalvars], "[[", j), "[[", models[i]))
      
      # Define the y_test and x_test matrices
      x_test = do.call(cbind, lapply(all_preds[finalvars], "[[", models[i]))[j,]
      y_test = all_preds[finalvars][[1]][j,1]
      
      # Append all the relevant x values
      train = as.data.frame(cbind(y_train, x_train))
      test = as.data.frame(rbind(c(y_test, x_test)))
      colnames(train)[1] = colnames(test)[1] = "class"
      colnames(train)[2:ncol(train)] = colnames(test)[2:ncol(test)] = paste0("_", colnames(train)[2:ncol(train)])
      
      # Run the models
      outCV = train.test.CV(train, test, models[i])
      meanpred = rbind(c(test[,1], rowMeans(test[,2:ncol(test)])))
      colnames(meanpred) = c("class", paste0("mean", models[i]))
      
      # Run the models and save the predicted and actual classes
      CVpreds[[(i-1)*2+1]] = rbind(CVpreds[[(i-1)*2+1]], outCV$preds)
      CVpreds[[(i-1)*2+2]] = as.data.frame(rbind(CVpreds[[(i-1)*2+2]], meanpred))
      
      # Run evaluation
      CVeval[[(i-1)*2+1]] = eval.CV(CVpreds[[(i-1)*2+1]], models[i])
      CVeval[[(i-1)*2+2]] = eval.CV(CVpreds[[(i-1)*2+2]], paste0("mean", models[i]))
    }
    names(CVpreds)[(i-1)*2+1] = names(CVeval)[(i-1)*2+1] = models[i]
    names(CVpreds)[(i-1)*2+2] = names(CVeval)[(i-1)*2+2] = paste0("mean", models[i])
    print(CVeval)
  }
  # Save the results
  save(CVeval, CVpreds, file = paste0("sis", nvar, "ens.rda"))
}

# This function will use in sample fit to build the best meta classifier model
fit.combine.preds = function(all_fitted, all_preds, nvar){
  # all_fitted = list of all fiited values, each element is a different FC condition
  # all_preds = list of all predicted values, each element is a different FC condition
  # nvars = number of variables to recruit from SIS
  
  # Save cpDFC2
  cpDFC2 = all_fitted[names(all_fitted) == "cpDFC2"]
  nobs = length(cpDFC2[[1]])
  
  # Loop through all the models
  models = names(all_fitted[[1]][[1]])
  models = models[models != "class"]
  CVpreds = CVeval = vector("list", 2*length(models))
  for(i in 1:length(models)){
    # Loop through all CV folds
    for(j in 1:nobs){
      # Select the current models x values
      x_train = all_fitted[!names(all_fitted) == "cpDFC2"]
      x_train = do.call(cbind, lapply(lapply(x_train, "[[", j), "[[", models[i]))
      
      # Define the y_train vector
      y_train = do.call(cbind, lapply(lapply(cpDFC2, "[[", j), "[[", 1))
      colnames(y_train) = "class"
      
      # Select based on F1
      results = cbind(y_train, x_train)
      results = eval.CV(results, colnames(results)[colnames(results) != "class"])
      results = results$Model[order(results$F1, decreasing = TRUE)][1:nvar]
      finalvars = names(all_preds) %in% results
      x_train = x_train[,colnames(x_train) %in% results]
      
      # Define the y_test and x_test matrices
      x_test = do.call(cbind, lapply(all_preds[finalvars], "[[", models[i]))[j,]
      y_test = all_preds[finalvars][[1]][j,1]
      
      # Append all the relevant x values
      train = as.data.frame(cbind(y_train, x_train))
      test = as.data.frame(rbind(c(y_test, x_test)))
      colnames(train)[1] = colnames(test)[1] = "class"
      colnames(train)[2:ncol(train)] = colnames(test)[2:ncol(test)] = paste0("_", colnames(train)[2:ncol(train)])
      
      # Run the models
      outCV = train.test.CV(train, test, models[i])
      meanpred = rbind(c(test[,1], rowMeans(test[,2:ncol(test)])))
      colnames(meanpred) = c("class", paste0("mean", models[i]))
      
      # Run the models and save the predicted and actual classes
      CVpreds[[(i-1)*2+1]] = rbind(CVpreds[[(i-1)*2+1]], outCV$preds)
      CVpreds[[(i-1)*2+2]] = as.data.frame(rbind(CVpreds[[(i-1)*2+2]], meanpred))
      
      # Run evaluation
      CVeval[[(i-1)*2+1]] = eval.CV(CVpreds[[(i-1)*2+1]], models[i])
      CVeval[[(i-1)*2+2]] = eval.CV(CVpreds[[(i-1)*2+2]], paste0("mean", models[i]))
    }
    names(CVpreds)[(i-1)*2+1] = names(CVeval)[(i-1)*2+1] = models[i]
    names(CVpreds)[(i-1)*2+2] = names(CVeval)[(i-1)*2+2] = paste0("mean", models[i])
    print(CVeval)
  }
  # Save the results
  save(CVeval, CVpreds, file = paste0("fit", nvar, "ens.rda"))
}

# This function will use LOOCV and fit to determine which models to combine together
loocv.combine.preds = function(FCdata, nens, dimredparams){
  # FCdata = list of all FCdata
  # nens = vector of number of ensembles to try 
  # dimredparams = dimension reduction parameters
  
  # Create a list which will hold the results of the ensembles
  ensdata = vector("list", length(nens))
  names(ensdata) = nens
  
  # Pre-processing
  for(i in 1:length(FCdata)){
    FCdata[[i]] = FCdata[[i]][,!(grepl("scan", colnames(FCdata[[i]])) | grepl("rank", colnames(FCdata[[i]])) | 
                         grepl("start", colnames(FCdata[[i]])) | grepl("Tend", colnames(FCdata[[i]])))]  
  }
  
  # Define data length
  datalen = max(unique(FCdata[[1]]$subj))
  
  # Register the parallel backend
  numcores = detectCores()
  registerDoParallel(cl = numcores)
  print(paste("Number of Cores:", numcores))
  
  # Perform CV
  models = c("rbfsvm", "linsvm", "logit", "tree")
  folds = unique(FCdata[[1]]$subj)
  modresults = ensresults = ensevals = vector("list", length(models))
  names(modresults) = models
  class = c()
  for(i in folds){
    # Set the seed based on the index
    set.seed(i*12345)
    
    # Print progress
    print(paste0("###CV Fold: ", i))
    
    # Split into train and test sets, keep only the first scan
    trainsubj = (1:length(folds))[-i]
    testsubj = i
    
    # Loop through all training subjects
    CVpreds = CVeval = vector("list", length(FCdata))
    names(CVpreds) = names(CVeval) = names(FCdata)
    k = 1
    for(j in trainsubj){
      # Print progress
      print(paste0("inner CV Fold: ", k))
      
      # Split into inner train and test sets
      inner_trainsubj = trainsubj[trainsubj != j]
      inner_testsubj = j
      
      # Apply variable selection
      splitFCdata = mclapply(FCdata, var.sel, inner_trainsubj, inner_testsubj, dimredparams)
      
      # Apply model training
      outCV = mclapply(splitFCdata, train.test.CV)
      
      # Save the predicted and actual classes
      currpreds = lapply(outCV, "[[", "preds")
      CVpreds = Map(rbind, CVpreds, currpreds)
    
      # Iterate
      k = k + 1
    }
    # Evaluate results
    CVeval = mclapply(CVpreds[!names(CVpreds) %in% c("SFC", "cpDFC1", "cpDFC2")], eval.CV)
    CVeval = Map(cbind, names(CVeval), CVeval)
    CVeval = lapply(CVeval, setNames, c("Group", colnames(CVeval[[1]])[2:ncol(CVeval[[1]])]))
    CVeval = do.call(rbind, CVeval)
    rownames(CVeval) = NULL
    
    # Loop through and use these results to select new model
    print("Selecting and training best models")
    models = unique(CVeval$Model)
    for(j in 1:length(models)){
      # Select the current model, and organize by Accuracy scores
      currmodel = CVeval[CVeval$Model == models[j],]
      currmodel = currmodel[order(currmodel$Accuracy, decreasing = TRUE),]
      
      # Find the maximum ensemble size
      maxnens = max(nens) - 1
      
      # Subselect FCdata
      splitFCdata = FCdata[names(FCdata) %in% c("cpDFC2", currmodel$Group[1:maxnens])]
      
      # Apply variable selection
      splitFCdata = mclapply(splitFCdata, var.sel, trainsubj, testsubj, dimredparams)
      
      # Apply model training
      outCV = mclapply(splitFCdata, train.test.CV, models[j])
      preds = sapply(lapply(outCV, "[[", 1), "[[", 2)
      
      # Append to all results
      modresults[[j]] = append(modresults[[j]], list(preds))
    }
    # Save the class of the current subject
    class = c(class, splitFCdata[[1]]$test$class)
  }
  # Loop through ensembles
  print("Evaluating all ensembles")
  currpreds = class
  k = 2
  for(i in nens){
    # Loop through all models
    for(j in 1:length(modresults)){
      currpreds = cbind(currpreds, rowMeans(do.call(rbind, modresults[[j]])[,1:i]))
      colnames(currpreds)[k] = paste0(i, "ens_", names(modresults)[j])
      k = k + 1
    }
  }
  colnames(currpreds)[1] = "class"
  
  # Evaluate results
  CVeval = eval.CV(currpreds)
  
  # Save all results
  save(CVeval, file = "ensresults.rda")
}

# This function is used for variable selection
var.sel = function(FCdata_iter, inner_trainsubj, inner_testsubj, dimredparams){
  # FCdata_iter = ith selection of FCdata
  # inner_trainsubj = training subjects for the inner CV fold
  # inner_testsubj = test subjects for the inner CV fold
  # dimredparams = parameters for SIS
  
  # Define the data to be subselected for
  data2subset = FCdata_iter[FCdata_iter$subj %in% inner_trainsubj, (grepl("_", colnames(FCdata_iter)) | 
                                                      colnames(FCdata_iter) == "class")]
  
  # Use the SIS function to screen
  finalvars = SIS(as.matrix(data2subset[, grepl("_", colnames(data2subset))]),
                  data2subset$class, family = "binomial", iter = dimredparams$iter,
                  standardize = TRUE, nsis = dimredparams$nsis)
  finalvars = finalvars$ix0 
  keepcols = data2subset[, grepl("_", colnames(data2subset))]
  keepcols = c("subj", "class", colnames(keepcols[, finalvars]))
  
  # Apply this to the train and test subject data
  train = as.data.frame(FCdata_iter[inner_trainsubj, colnames(FCdata_iter) %in% keepcols, drop = FALSE])
  test = as.data.frame(FCdata_iter[inner_testsubj, colnames(FCdata_iter) %in% keepcols, drop = FALSE])
  
  # Return as a list
  splitFCdata_iter = list(train, test)
  names(splitFCdata_iter) = c("train", "test")
  return(splitFCdata_iter)
}
