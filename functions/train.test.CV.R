# Load libraries
library(e1071)
library(rpart)

# This function takes in data and performs cross validation across the data
train.test.CV = function(splitFCdata, models = NULL){
  # splitFCdata = contains train and test data
  # models = models to run
  
  # Define train and test sets
  train = splitFCdata$train
  test = splitFCdata$test
  
  # Check on models
  if(is.null(models)){
    models = c("rbfsvm", "linsvm", "logit", "tree")
  }
  
  # Remove the subject number
  train = train[,colnames(train) != "subj"]
  test = test[,colnames(test) != "subj"]
  
  # Define the x and y for training and test splits
  xtrain = as.matrix(train[,grepl("_", colnames(train))])
  ytrain = as.factor(train$class)
  xtest = as.matrix(test[,grepl("_", colnames(train))])
  ytest = as.factor(test$class)
  
  # Fit models to the training set, test on test set
  fitted = preds = c()
  for(i in 1:length(models)){
    currmod = models[i]
    if(currmod == "rbfsvm"){
      fitted_model = svm(x = xtrain, y = ytrain, type = "nu-classification", probability = TRUE)
      pred = predict(fitted_model, newdata = xtest, probability = TRUE)
      pred = attr(pred, "probabilities")[,1]
      fit = predict(fitted_model, newdata = xtrain, probability = TRUE)
      fit = attr(fit, "probabilities")[,1]
    } else if(currmod == "linsvm"){
      fitted_model = svm(x = xtrain, y = ytrain, type = "nu-classification", 
                         kernel = "linear", probability = TRUE)
      pred = predict(fitted_model, newdata = xtest, probability = TRUE)
      pred = attr(pred, "probabilities")[,1]
      fit = predict(fitted_model, newdata = xtrain, probability = TRUE)
      fit = attr(fit, "probabilities")[,1]
    } else if(currmod == "logit"){
      fitted_model = glm(class ~ ., family = binomial, data = train)
      pred = predict(fitted_model, test, type = "response")
      fit = predict(fitted_model, train, type = "response")
    } else if(currmod == "tree"){
      fitted_model = rpart(class ~ ., method = "class", data = train)
      pred = predict(fitted_model, test, type = "prob")[,2]
      fit = predict(fitted_model, train, type = "prob")[,2]
    }
    
    # Append to preds and fitted
    preds = data.frame(cbind(preds, t(pred)))
    fitted = data.frame(cbind(fitted, fit))
    colnames(preds)[i] = colnames(fitted)[i] = models[i]
  }
  # Append true class to each
  preds = cbind(test$class, preds)
  fitted = cbind(train$class, fitted)
  colnames(preds)[1] = colnames(fitted)[1] = "class"
  
  # Append to all results
  outCV = list(preds, fitted)
  names(outCV) = c("preds", "fitted")
  
  # Return all results
  return(outCV)
}
