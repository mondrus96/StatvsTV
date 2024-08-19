#===========================================================================
# Libraries
#===========================================================================
library(MASS)
library(matlab)
library(parallel)
library(igraph)
library(glasso)

#===========================================================================
# Functions
#===========================================================================
# The main function that calls all other functions

spectralfunction = function(x,y,type,try.K,T,lower,upper,cp.index,maxn,quan,n.boot,block,quantile1,permtype,actualcps=NULL,bootcps=NULL,signcps=NULL,elapsed.time=NULL){

## Inputs: 
# x               = column matrix of time                
# y               = data set to be analyzed (timexROIs)
# type            = 1 the algorithm is exhausted first, then inference
#                 = 2 inference is performed after each cp is found       
# try.K           = prespecified integer for the number of clusters expected
# T               = number of time points
# lower           = first time point for the algorithm
# upper           = last time point for the algorithm
# cp.index        = list for output
# maxn            = minimum distance between change points
# quan            = quantile threshold to locate outliers
# n.boot          = number of bootstrap\permutation replicates
# block           = size of average block length for stationary bootstrap
# quantile1       = quantile for conf bounds adjusting for multiple comparisons
# permtype        = 1 for permutations
#                 = 2 stationary bootstrap
# actualcps       = file for actual cps
# bootcps         = file for bootstrap distribution at cps
# signcps         = file for significant cps
# elapsed.time    = file for elapsed time

# Record the starting time
t1 = proc.time()[3]

# Estimate change points
if (type == 1){
	cps = cp.spectral(x,y,try.K,lower,upper,cp.index,maxn,quan)
} else {
	cps = cp.spectral.all(x,y,try.K,lower,upper,cp.index,maxn,quan,n.boot,block,permtype)
}

# If there are any change points, continue to inference procedure
if (length(cps)>0){
	if (type == 1){
		cps = stat.boot.spect(cps,x,y,try.K,T,n.boot,block,permtype)
		for (ijc in 1:dim(cps$actual)[1]){
			print(paste("CP at =",cps$actual[ijc,1]))
			print(paste("result.diff =",cps$actual[ijc,2]))
		}
	}
	# save the change points and the bootstrap samples
	#write(t(cps$actual),file = actualcps,ncolumns = 2,sep = "\t",append=TRUE)
	#write(t(cps$boot),file = bootcps,ncolumns = 2,sep = "\t",append=TRUE)

	# Find signifcant change points and plot the graph between each pair of significant change points
	results  = sign.cps(y,T,cps,quantile1,n.boot)
	
	# Save significant cps
	if (length(results$sign.cps) > 0){
		#write(t(results$sign.cps),file = signcps,ncolumns = 3,sep = "\t",append=TRUE)
		} else {
			print("No Significant CPs")
	}
} else {
	print("No Splits")

}
# Record the starting time
#elapsed.t = proc.time()[3] - t1
#cat('Elapsed time: ', elapsed.t, ' seconds', '\n')
#write(elapsed.t,file = elapsed.time,ncolumns = 1,sep = "\t",append=TRUE)
return(list(cps = cps, results = results))
}

#===========================================================================
# This function detects change points based the angle between the
# eigenvectors by exhausting the algorithm first.

cp.spectral = function(x,y,try.K,lower,upper,cp.index,maxn,quan){
  
  ## Inputs: 
  # x               = column matrix of time                
  # y               = data set to be analyzed (timexROIs)      
  # try.K           = prespecified integer for the number of clusters expected
  # lower           = first time point for the algorithm
  # upper           = last time point for the algorithm
  # cp.index        = list for output
  # maxn            = minimum distance between change points
  
  low.s       = lower+maxn
  upp.s       = upper-maxn
  
  if (low.s <= upp.s) {
    
    indices = low.s:upp.s
    # Functions for below
    inputglassol = function(i,px,py,lower){
      crapl = which(px<=i & px>lower)
      return(py[crapl,])
    }
    inputglassor = function(i,qx,qy,upper){
      crapr = which(qx<=upper & qx>i)
      return(qy[crapr,])
    }
    
    if (length(indices) > 1){
      
      # Left indices 
      # For parallel programming
      clusterExport(cl, varlist=c("x","y","try.K"))
      clusterEvalQ(cl, source("../../functions/ParallelFunctions.R"))
      left.ind = parLapply(cl,indices,inputglassol,x,y,lower)
      temp_l   = parLapply(cl,left.ind,spectral.clustering,try.K)
      
      # Right indices
      right.ind = parLapply(cl,indices,inputglassor,x,y,upper)
      temp_r    = parLapply(cl,right.ind,spectral.clustering,try.K)
      
      rr = length(temp_l)
      mul.list = list()
      for (ll in 1:rr){
        mul.list[[ll]] = t(temp_l[[ll]])%*%temp_r[[ll]]
      }
      temp.svd = lapply(mul.list,function(xx){junk1=svd(xx);return(sum(junk1$d))})
      temp.svd = unlist(temp.svd)
      loc      = which.min(temp.svd)
      
      # Add new changepoint
      loc2 = loc
      new.cp = indices[loc]
      cc  = c(new.cp,temp.svd[loc2])
      cp.index$actual = rbind(cp.index$actual,cc)
      
      cp.index  = Recall(x,y,try.K,lower,new.cp,cp.index,maxn,quan)
      cp.index  = Recall(x,y,try.K,new.cp,upper,cp.index,maxn,quan)
      
      #if (length(temp.svd)>2){
      #    loc2     = select(temp.svd, quan)
      #} else {loc2 = loc}
      
      #if(min(temp.svd[loc2])>= 2){
        
      #  cp.index = cp.index
        
      #} else {
      #  new.cp = indices[loc2]
        #print(paste("CP at =", new.cp))
        #print(paste("result.diff =",temp.svd[loc2]))
      #  cc  = c(new.cp,temp.svd[loc2])
      #  cp.index$actual = rbind(cp.index$actual,cc)
        
      #  cp.index  = Recall(x,y,try.K,lower,new.cp,cp.index,maxn,quan)
      #  cp.index  = Recall(x,y,try.K,new.cp,upper,cp.index,maxn,quan)
        
      #} # if(min(temp.svd[loc2])>= 2) loop
      
    } else {cp.index=cp.index} # if (length(indices)>1) loop
    
  } # if low.s<upp.s loop
  
  else {cp.index=cp.index}	
  return(cp.index)	
  
}

#===========================================================================
# This function detects change points based the angle between the
# eigenvectors. It carries out the inference on each cp as they are found

cp.spectral.all = function(x,y,try.K,lower,upper,cp.index,maxn,quan,n.boot,block,permtype){

## Inputs: 
# x               = column matrix of time                
# y               = data set to be analyzed (timexROIs)      
# try.K           = prespecified integer for the number of clusters expected
# lower           = first time point for the algorithm
# upper           = last time point for the algorithm
# cp.index        = list for output
# maxn            = minimum distance between change points
# quan            = quantile threshold to locate outliers
# n.boot          = number of bootstrap\permutation replicates
# block           = size of average block length for stationary bootstrap
# permtype        = 1 for permutations
#                 = 2 stationary bootstrap

low.s       = lower+maxn
upp.s       = upper-maxn

if (low.s <= upp.s) {

	indices = low.s:upp.s
	# Functions for below
	inputglassol = function(i,px,py,lower){
		crapl = which(px<=i & px>lower)
		return(py[crapl,])
	}
	inputglassor = function(i,qx,qy,upper){
		crapr = which(qx<=upper & qx>i)
		return(qy[crapr,])
	}
	if (length(indices)>2){

		# Left indices 
		# For parallel programming
		clusterExport(cl, varlist=c("x","y","try.K"))
		clusterEvalQ(cl, source("../../functions/ParallelFunctions.R"))
		left.ind = parLapply(cl,indices,inputglassol,x,y,lower)
		temp_l   = parLapply(cl,left.ind,spectral.clustering,try.K)
	
		# Right indices
		right.ind = parLapply(cl,indices,inputglassor,x,y,upper)
		temp_r    = parLapply(cl,right.ind,spectral.clustering,try.K)
	
		rr = length(temp_l)
		mul.list = list()
		for (ll in 1:rr){
			mul.list[[ll]] = t(temp_l[[ll]])%*%temp_r[[ll]]
		}
		temp.svd = lapply(mul.list,function(xx){junk1=svd(xx);return(sum(junk1$d))})
		temp.svd = unlist(temp.svd)
		loc      = which.min(temp.svd)
		loc2     = select(temp.svd, quan)
	
		if(min(temp.svd[loc2])>= 2){cp.index = cp.index}
	
		else{
			new.cp = indices[loc2]
			print(paste("CP at =", new.cp))
			print(paste("result.diff =",temp.svd[loc2]))
			cc  = c(new.cp,temp.svd[loc2])
			cp.index$actual = rbind(cp.index$actual,cc)
			
			#####################################################
		      # Now do the permutation/bootstrapping for this cp

			ind.interval   = which(x<=upper & x>lower)

			bootsamp       = list()
			for (bj in 1:n.boot){
				set.seed(bj*2)
				if (permtype == 1){
					junk =  stat_boot_cp(y[ind.interval,],block)
				} else {
					junk = permute_cp(y[ind.interval,])
				}
				row.names(junk) = ind.interval
				bootsamp[[bj]] = junk
				}
			# Left and right subsets
			left.ind      = which(x<=new.cp & x>lower)
			l.ind         = length(left.ind)
			right.ind     = which(x<=upper & x>new.cp)
			r.ind         = length(right.ind)
			
			# Left samples
			lbootsampdata1 = parLapply(cl,bootsamp,function(yy,l.ind){return(yy[1:l.ind,])},l.ind)
			# Left estimate
			clusterEvalQ(cl, source("../../functions/ParallelFunctions.R"))
			hot1 = parLapply(cl,lbootsampdata1,spectral.clustering,try.K)
			
			remove(lbootsampdata1)
			# Right samples
			rbootsampdata1 = parLapply(cl,bootsamp,function(yy,l.ind,r.ind){return(yy[((l.ind+1):(l.ind+r.ind)),])},l.ind,r.ind)
			remove(bootsamp)
			# Right estimate
			clusterEvalQ(cl, source("../../functions/ParallelFunctions.R"))
			hot2 = parLapply(cl,rbootsampdata1,spectral.clustering,try.K)
			remove(rbootsampdata1)
			
			rr = length(hot2)
			mul.list = list()
			for (ll in 1:rr){
				mul.list[[ll]] = t(hot1[[ll]])%*%hot2[[ll]]
			}
			temp.svd      = lapply(mul.list,function(xx){junk1=svd(xx);return(sum(junk1$d))})
			temp.svd      = unlist(temp.svd)
			ddd           = cbind(rep(new.cp,n.boot),temp.svd) 
			cp.index$boot = rbind(cp.index$boot,ddd)

			cp.index  = Recall(x,y,try.K,lower,new.cp,cp.index,maxn,quan,n.boot,block,permtype)
			cp.index  = Recall(x,y,try.K,new.cp,upper,cp.index,maxn,quan,n.boot,block,permtype)
			
			} # else loop
			
		} # if (length(indices)>1) loop
	} # if low.s<upp.s loop

else{cp.index=cp.index}	
return(cp.index)	
}

#===========================================================================
# This function calculates the stationary or permutation distribution at 
# each candidate change point

stat.boot.spect = function(ncps,nx,ny,ntry.K,nT,nnboot=1000,nblock,npermtype){

## Inputs: 
# ncps     = candidate change points
# nx       = column matrix of time                
# ny       = data set to be analyzed (timexROIs)     
# ntry.K   = prespecified integer for the number of clusters expected
# nT       = number of time points
# nnboot   = number of bootstrap\permutation replicates
# block    = size of average block length for stationary bootstrap
# permtype = 1 for permutations
#          = 2 stationary bootstrap

results = list()
cptemp  = ncps$actual[,1]

order.cp = c()
ccc = order(cptemp)
for (i in 1:length(cptemp)){	
	order.cp = c(order.cp,cptemp[ccc[i]])
	}

# Order them
order.cp = c(0,order.cp)
order.cp  = c(order.cp,nT)

attributes(order.cp) = NULL

# Parallel programming
indices1   = 1:(length(order.cp)-1)
inputglasso1 = function(i,xx,yy,order.cp){
	crap1 = which(xx<=order.cp[i+1] & xx>order.cp[i])
	return(yy[crap1,])
}
# Reestimating the angle between candidate change points
clusterEvalQ(cl, source("../../functions/ParallelFunctions.R"))
data1 = parLapply(cl,indices1,inputglasso1,nx,ny,order.cp)

clusterEvalQ(cl, source("../../functions/ParallelFunctions.R"))
hot1 = parLapply(cl,data1,spectral.clustering,ntry.K)

rr = length(hot1)
mul.list = list()
for (ll in 1:(rr-1)){
	mul.list[[ll]] = t(hot1[[ll]])%*%hot1[[(ll+1)]]
}
temp.svd = lapply(mul.list,function(xx){junk1=svd(xx);return(sum(junk1$d))})
temp.svd = unlist(temp.svd)

results$actual = cbind(order.cp[2:(length(order.cp)-1)],temp.svd) 

# Now perform inference using stationary bootstrap
for (aa in 1:(length(order.cp)-2)){
	
	lower1         = order.cp[aa]
	upper1         = order.cp[aa+2]
	ind.interval   = which(x<=upper1 & x>lower1)

	bootsamp       = list()
	for (bj in 1:nnboot){
		set.seed(bj*2)
		if (npermtype == 1){
			junk =  stat_boot_cp(ny[ind.interval,],nblock)
		} else {
			junk = permute_cp(ny[ind.interval,])
		}
		row.names(junk) = ind.interval
		bootsamp[[bj]] = junk
	}
	# Left and right subsets
	left.ind      = which(nx<=order.cp[aa+1] & nx>order.cp[aa])
	l.ind         = length(left.ind)
	right.ind     = which(nx<=order.cp[aa+2] & nx>order.cp[aa+1])
	r.ind         = length(right.ind)

	# Left samples
	lbootsampdata1 = parLapply(cl,bootsamp,function(yy,l.ind){return(yy[1:l.ind,])},l.ind)
	# Left estimate
	clusterEvalQ(cl, source("../../functions/ParallelFunctions.R"))
	hot1 = parLapply(cl,lbootsampdata1,spectral.clustering,ntry.K)

	# Right samples
	rbootsampdata1 = parLapply(cl,bootsamp,function(yy,l.ind,r.ind){return(yy[((l.ind+1):(l.ind+r.ind)),])},l.ind,r.ind)
	# Right estimate
	clusterEvalQ(cl, source("../../functions/ParallelFunctions.R"))
	hot2 = parLapply(cl,rbootsampdata1,spectral.clustering,ntry.K)
	
	rr = length(hot2)
	mul.list = list()
	for (ll in 1:rr){
		mul.list[[ll]] = t(hot1[[ll]])%*%hot2[[ll]]
	}
	temp.svd     = lapply(mul.list,function(xx){junk1=svd(xx);return(sum(junk1$d))})
	temp.svd     = unlist(temp.svd)
	ddd          = cbind(rep(order.cp[aa+1],nnboot),temp.svd) 
	results$boot = rbind(results$boot,ddd)

	} # aa loop
return(results)
}


#=============================================================================
# This function calculates the evaluation value for each partition point

partition.value = function(Y,try.K,t.l.min,t.l.max,t.r.min,t.r.max){

## Inputs: 
# Y        = data set to be analyzed (T by p matrix)
# try.K    = prespecified integer for the number of clusters expected
# t.l.min  = minimum time for segment 1
# t.l.max  = maximum time for segment 1
# t.r.min  = minimum time for segment 2
# t.r.max  = maximum time for segment 2

V = spectral.clustering(Y = Y, try.K = try.K, t.min = t.l.min, t.max = t.l.max)
W = spectral.clustering(Y = Y, try.K = try.K, t.min = t.r.min, t.max = t.r.max)
	
temp = t(V) %*% W
temp.svd = svd(temp)
return(sum(temp.svd$d))
}

#============================================================================
# This functions selects the minimum value based 

select = function(ix,iquan){
	
## Inputs:
# ix    = data matrix
# iquan = quantile threshold
 
len = length(ix)
diff.seq = abs(ix[-1] - ix[-len])
thresh = quantile(diff.seq, prob = c(iquan))
#print(thresh)
outlier = c()
counter = 1
for(i in 2:(len - 1)){
	if(abs(ix[i-1] - ix[i]) >= thresh & abs(ix[i+1] - ix[i]) >= thresh){
		outlier[counter] = i
		counter = counter + 1
		}
	}
max.value = max(ix) + 10
temp = ix
temp[outlier] = max.value
return(which.min(temp))
	
}

#===========================================================================
# This function permutes the data set
# The first element of the permutation will be a 1*n matrix where n is the 
# number of ROIs

permute_cp = function(yayadata){

## Inputs:
# yayadata = data set to be analyzed

n1   = dim(yayadata)[1]
k1   = dim(yayadata)[2] # Number of regions

b.y = c()
index45 = sample(seq(1,n1,1))
for (i in 1:n1){
	b.y    = rbind(b.y,yayadata[index45[i],])
	}
return(as.matrix(b.y))
}

#============================================================================
# This function stationary bootstraps the whole data set 
# The first block element of the resampling will be a block*n matrix where n 
# is the number of ROIs

stat_boot_cp = function(data,block){

## Inputs:
# data   = data set to be analyzed 
# block  = size of average block length for stationary bootstrap
#          (given in % terms)

n1   = dim(data)[1]
k1   = dim(data)[2] # Number of subjects

block      = round(n1*block)
ytemp      = data
b.x        = NULL
b.y        = NULL
r.num      = runif(n1, min=0, max=1)
index      = sample(seq(1,n1,1))[1]
temp_index = c(seq(1,n1,1),seq(1,n1,1))
b.x        = c(b.x,index)

for (i in 1:(n1-1)){
	if (r.num[i] < 1 - (1/block)){
			
		b.x   = c(b.x,temp_index[index + 1])
		index = index + 1
			
		}
	else
	if (r.num[i] >= 1/block){
		index  = sample(seq(1,n1,1))[1]
		b.x    = c(b.x,index)
            }
      }
b.x   = b.x[1:n1]
b.y   = rbind(b.y,ytemp[b.x,])
return(b.y)
}

#===========================================================================
# Finds significant change points using the stationary bootstrap or
# permutation distribution

sign.cps = function(ydata,yT,ycp.index,yquantile1,ynnboot=1000){

## Inputs: 
# ydata        = data set to be analyzed (timexROIs)
# yT           = number of time points
# ycp.index    = actual cps and boot cps results from above functions
# yquantile1   = quantile for conf bounds adjusting for multiple comparisons  
# ynnboot      = number of bootstrap\permutation replicates

newresults = list()
output     = as.matrix(ycp.index$actual)

junk = rep(0,yT)
n    = dim(output)[1]
for (i in 1:n){
	junk[c(output[i,1])] = output[i,2]
	}

stat     = as.matrix(ycp.index$boot)
cp.times = output[,1]
pvalues  = c()
for (i in c(cp.times)){
	index1 = which(stat[,1] == i)
	index2 = which(output[,1] == i)
	pvalues = c(pvalues,length(which(stat[index1,2]<= output[index2,2]))/ynnboot)
	}

# Multiple comparisons
pp = p.adjust(pvalues, method = "BH", n = length(pvalues))

junk.output = c()
for (i in 1:dim(output)[1]){
	if (pp[i] < yquantile1){
		seb1 = c(output[i,],pp[i])
		junk.output = rbind(junk.output,seb1)
		}					
	}

newresults$sign.cps = junk.output
return(newresults)
}


#============================================================================
#============================================================================
