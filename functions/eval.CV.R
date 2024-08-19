# Load libraries
library(caret)

# This function evaluates the results of the cross validation fold
eval.CV = function(results, models = NULL){
  # results = results of prediction
  # models = models to run
  
  # Initialize the output
  CVevals = c()
  if(is.null(models)){
    models = colnames(results)
    models = models[models != "class"]
  }
  
  # Loop through
  for(i in models){
    # Select relevant results for model
    model_results = results[,colnames(results) %in% c("class", i)]
    
    # Evaluate with confusion matrix
    curr.metrics = confusionMatrix(factor((1*(model_results[,2] > 0.5)), levels = c(0,1)),
                                   factor(model_results[,1], levels = c(0,1)), mode = "everything", positive="1")
    
    # Put into a dataframe format
    curr.row = as.data.frame(t(c(i, curr.metrics$overall, curr.metrics$byClass)))
    CVevals = rbind(CVevals, curr.row)
  }
  # Rename the first column
  colnames(CVevals)[1] = "Model"
  
  # Return the final output 
  return(CVevals)
}