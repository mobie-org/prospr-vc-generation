################################################################################################
#Splits a matrix into two
#
#Parameters:
# - mat -> matrix to cluster
#Returns the indexes of the matrix (rownames) divided
################################################################################################
library(vegan)
splitMat <- function(mat){
  # set custom distance and clustering functions
  hclustfunc <- function(x) hclust(x, method="complete")
  #distfunc <- function(x) dist(x,method="euclidean")
  distfunc <-function(x) vegdist(x, method="jaccard")
  #distfunc <-function(x) bcdist(x)
  
  if(nrow(mat)<20000){ #the matrix size can be computed with vegdist
    distance_matrix <- distfunc(mat)
  }else{ #use another function to compute the matrix
    #calculate the blocks subdivision
    for(nB in c(5:100)){
      if (nrow(mat) %% nB == 0){break}
    }
    if (nB==100){ #no possible subdivision
      distance_matrix <- distfunc(mat)
    }else{
      distance_matrix <- 1-bigcor(t(mat),nblocks = nB)
      colnames(distance_matrix) <- rownames(mat)
      rownames(distance_matrix) <- rownames(mat)
      distance_matrix <- as.dist(distance_matrix)
    }
  }
  
  # obtain the clusters
  nclusters <- 2
  fit <- hclustfunc(distance_matrix)
  clusters <- cutree(fit, nclusters) 
  
  # return two vectors in a list
  list2return <- list(names(clusters[clusters==1]),names(clusters[clusters==2]))
  names(list2return) <- c("A","B")
  
  list2return
}

################################################################################################
#Calculates the correlation of a list of supervoxels
#
#Parameters:
# - SVvec -> vector containing Supervoxels IDs
# - mat -> expression matrix
#Returns the correlation value of the supervoxels
################################################################################################
MyCorrFunc <- function(SVvec,mat){
  submat <- mat[SVvec,]
  cormat <- cor(t(initMat[SVvec,]),use="complete.obs",method="kendall")
  diag(cormat) <- NA
  mean(cormat,na.rm=T)  
}


################################################################################################
#Calculates the spatial coherence of a list of supervoxels
#
#Parameters:
# - SVvec -> vector containing Supervoxels IDs
# - coord -> coordinates of supervoxels
#Returns a spatial correlation value of the supervoxels, taking on account the bilaterality
################################################################################################
SpatCoheFunc <- function(SVvec,coord){
  #"fold" the x axis to make it bilateral, and calculate the distances between the coordinates
  #the coords are ordered y,x,z
  #calculate y axis mean point. THIS ASSUMES THE REFERENCE IS SIMETRICAL
  yaxMP <- mean(c(max(coord[,2]),min(coord[,2])))
  #recalibrate y coordinates to overlap both halves
  coord[,2] <- abs(coord[,2]-yaxMP) + yaxMP
  
  #calculate a distance matrix between the coordinates and get the mean
  subcoord <- coord[SVvec,,drop=F]
  distmat <- dist(subcoord, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  mean(distmat,na.rm=T)
}



################################################################################################
#Splits a big matrix to calculate correlations
#slightly modified to not use the ff package
#From: http://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/
################################################################################################

bigcor <- function(x, nblocks = 10, verbose = TRUE, ...)
{
  #library(ff, quietly = TRUE)
  NCOL <- ncol(x)
  
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  #corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  corMAT <- matrix(NA,NCOL,NCOL)
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2], ...)
    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }
  
  gc()
  return(corMAT)
}
