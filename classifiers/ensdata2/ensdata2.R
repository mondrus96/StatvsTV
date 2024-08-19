### For All Experiments ###
# Define the directory
dir = "data2"

# Load all the functions
filesrcs = list.files(path = "../../functions", pattern="*.R$")
sapply(paste0("../../functions/", filesrcs), source)

# Order the predictions by accuracy
allres = list.files(paste0("../", dir))
allres = allres[grepl(".rda", allres)]
allres = allres[!grepl("CRMT", allres)] # remove CRMT
allres = allres[!grepl("NCPD", allres)] # remove NCPD
modacc = c()
for(j in 1:length(allres)){
  load(paste0("../", dir, "/", allres[j]))
  crit = CVeval$F1
  names(crit) = CVeval$Model
  nme = strsplit(allres[j], ".rda")[[1]][1]
  modacc = rbind(modacc, rbind(c(nme, crit)))
}
# Change crit to numeric and order by it
modacc = as.data.frame(modacc)
modacc[,2:ncol(modacc)] = apply(modacc[,2:ncol(modacc)], 2, as.numeric)
colnames(modacc)[1] = "Name"
modacc[order(modacc$linsvm, decreasing=TRUE),][1:10,]

# Loop through predictions
for(j in 2:10){
  mean.combine.preds(dir, modacc, 1:j)
}

# Find which predictors to use
preds = 1:10
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
}
allpreds = cbind(CVpreds$class, allpreds$linsvm)

# Convert matrix to data frame
df <- as.data.frame(allpreds)
colnames(df) <- c("Class", paste0(1:(ncol(df) - 1)))

# Check statistical test
numtests <- (ncol(df) - 1)
testpreds <- allpreds[,2:ncol(allpreds)]
p_values <- matrix(NA, nrow = numtests, ncol = numtests)
rownames(p_values) <- colnames(testpreds)
colnames(p_values) <- colnames(testpreds)

# Iterate through all combinations of columns
k = 0
for (i in 1:(numtests - 1)){
  for (j in (i + 1):numtests){
    test_result <- cor.test(testpreds[, i], testpreds[, j])
    p_values[i, j] <- test_result$p.value
    p_values[j, i] <- test_result$p.value  # Symmetric matrix
    k = k + 1
  }
}
tests <- (p_values < 0.05/k)
sum(tests[upper.tri(tests)])
tests

# Print the matrix of p-values
print(p_values)

# Convert Class labels
df$Class <- factor(df$Class, levels = c(0, 1), labels = c("CN", "eMCI"))

# Exclude the first column
df_subset <- df[,colnames(df) != "Class"]

# Compute the correlation matrix
cor_matrix <- cor(df_subset)

# Plot the correlation matrix
png("predcorrs.png", width = 7, height = 5, units = "in", res = 300, pointsize = 10)
ggcorrplot(cor_matrix, 
           method = "circle", 
           type = "lower", 
           lab = TRUE) +
  scale_fill_gradient2(low = "#0C6291", 
                       high = "#A63446", 
                       mid = "white", 
                       midpoint = 0, 
                       limit = c(-1, 1), 
                       name = "Correlation") +
  scale_x_continuous(breaks = 1:ncol(df_subset), labels = 1:ncol(df_subset)) +
  scale_y_continuous(breaks = 1:ncol(df_subset), labels = 1:ncol(df_subset))
dev.off()
