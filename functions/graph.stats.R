# Load library
library(igraph)
library(parallel)
library(DirectedClustering)

# This function will be used to generate graph summary statistics from the estimated networks
graph.stats = function(data, graphtype = "all"){
  # data = input data in tabular format, where each row is a stationary graph, and columns are edges
  # graphtype = type of graph stats to calculate, by default is local clustering coefficient
  
  # Set the seed based on the length of data
  set.seed(54321)
  
  # Change data to list format
  varnames = colnames(data)
  data = as.list(as.data.frame(t(data)))
  names(data) = NULL
  for(i in 1:length(data)){
    names(data[[i]]) = varnames
  }
  
  # Apply function
  output = mclapply(data, graph.stats.each, graphtype)
  output = do.call(rbind, output)
  
  # Clean up and make it a dataframe
  rownames(output) = NULL
  output = as.data.frame(output)
  
  # Return all the stats
  return(output)
}

# This is a helper function to be applied over each row to calculate graphical statistics
graph.stats.each = function(rowvals, graphtype){
  # rowvals = a vector of edge weights for the overall graph
  # graphtype = type of graph stats to calculate, by default is local clustering coefficient
  
  # Convert data into adjacency matrix
  nodes = unlist(strsplit(unlist(strsplit(names(rowvals), ":")), "_"))
  nodes = unique(nodes[nodes != ""])
  adjmat = matrix(0, length(nodes), length(nodes))
  adjmat[upper.tri(adjmat)] = rowvals
  adjmat = t(adjmat)
  adjmat[upper.tri(adjmat)] = rowvals
  colnames(adjmat) = rownames(adjmat) = nodes
  diag(adjmat) = 0
  
  # Convert the adjacency matrix into a graph
  g = graph_from_adjacency_matrix(adjmat, mode = "upper")
  
  # Compute all the relevant metrics required
  clcoef = transitivity(g, type = "local", isolates = "zero")
  names(clcoef) = paste0("_clustcoef", 1:length(clcoef))
  asst = assortativity_degree(g, directed = FALSE)
  names(asst) = "_assort"
  localeff = local_efficiency(g, directed = FALSE)
  names(localeff) = paste0("_localeff", 1:length(localeff))
  deg = degree(g)
  names(deg) = paste0("_deg", 1:length(deg))
  short = all_shortest_paths(g, V(g))$nrgeo
  names(short) = paste0("_shorpath", 1:length(short))
  betw = betweenness(g)
  names(betw) = paste0("_betw", 1:length(betw))
  
  # Put them together into a final vector
  graphstats = c(clcoef, asst, localeff, deg, short, betw)
  
  # Return the output
  return(graphstats)
}
