#============================================================================
# Parallel Functions
# These functions are for the cluster

#============================================================================
# This function calculates spectral clustering for a given data set and
# number of clusters

spectral.clustering = function(ny,nK){

## Inputs: 
# ny     = data set to be analyzed
# nK     = prespecified integer for the number of clusters expected
# 	     N.B. to be safe, one can set K slightly larger

sub.y      = ny
sigma      = cor(sub.y)
eig.adj    = eigen(sigma)
u.mat      = eig.adj$vectors[,c(1:nK)]
fit.kmeans = kmeans(u.mat, centers = nK)

temp1       = matrix(0, nrow = length(fit.kmeans$cluster), ncol = dim(fit.kmeans$centers)[2])
for(jj in 1:length(fit.kmeans$cluster)){
	set.seed(jj)
	ind = fit.kmeans$cluster[jj]
	temp1[jj,] = fit.kmeans$centers[ind,]
	}
return(temp1)
		
}

#===========================================================================
# This functions is used in the spectral.clustering function above

rep.kmeans = function(cluster, center){

## Inputs: 
# cluster = clusters
# center  = centre of cluster

temp2 = matrix(0, nrow = length(cluster), ncol = dim(center)[2])
for(jj in 1:length(cluster)){
	ind = cluster[jj]
	temp2[jj,] = center[ind,]
	}
return(temp2)
}
#============================================================================
