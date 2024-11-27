################################################################################
## Function for the simulation for one iteration; the function calculate      ##
## the following variogram estimator for the fourth direciton (S-N, E-W,      ##
## SW-NE, SE-NW) for one data set: Matheron, Genton, raw/reweighted MCD.org,  ##
## raw/reweighted MCD.diff, modified raw/reweighted MCD.org and modified raw/ ##
## reweighted MCD.diff; the modified versions are only calculated for the     ##
## S-N and E-W direction and only if the argument mod is TRUE (no default)    ##
################################################################################

## Input: data:    vector   - nx X ny dimensional; containing the data 
##        grid:    matrix   - (nx X ny) rows; 2 columns; first column contains the 
##                            x-coordinate of the observations and the second column
##                            contains the y-coordinate of the observations
##        hmax:    vector   - 4 dimensional; containing for each direction the number
##                            of lags to be estimated (order: S-N, E-W, SE-NW, SW-NE)
##        cp:      numeric  - correction factor; default: 1
##        mod:     logical  - TRUE if the modified versions of the MCD variogram 
##                            estimators should be calculated (see section 4.3);
##                            in this case all estimators are only calculated for 
##                            the S-N and E-W direction; default: FALSE
##        mx:       numeric - dependency in the E-W direction (along the x-axis);
##                            for the modified version mx dependency in the E-W
##                            direction is assumed (see section 4.3); default: NULL
##        my:       numeric - dependency in the N-S direction (along the y-axis);
##                            for the modified version my dependency in the N-S
##                            direction is assumed (see section 4.3); default: NULL
##        missing:  logical - TRUE if the data contains missing values; the 
##                            columns with missing values will be delete;
##                            default: FALSE
##        det:      logical - TRUE if the deterministic MCD should be used
##                            instead of FASTMCD (see also explanations in 
##                            covMcd of the package robustbase);
##                            default: FALSE
## Output: res: list - with four entries (if mod = TRUE only two entries); for
##                     each direction one entry (order: S-N, E-W, SE-NW, SW-NE); 
##                     each entry is a matrix with hmax rows and six columns 
##                     (if mod = TRUE ten columns);each column contains the result 
##                     of one estimator in the following order: MCD.diff, MCD.diff.re, 
##                     MCD.org, MCD.org.re, Matheron, Genton, (MCD.diff.mod, 
##                     MCD.diff.mod.re, MCD.org.mod, MCD.org.mod.re); 


sim.it <- function(data, 
                   grid, 
                   hmax, 
                   cp = 1, 
                   mod = FALSE, 
                   mx = NULL, 
                   my = NULL, 
                   missing = FALSE, 
                   det = FALSE){
  
  
  nx <- length(unique(grid[,1])) # number of x-coordinates, i. e. number of coordinates in 
                                 # east-west direction
  ny <- length(unique(grid[,2])) # number of y-coordinates, i. e. number of coordinates in 
                                 # north-south direction
  
  # Save the data as matrix with ny rows and nx columns
  # the element in the i-th row and j-th column is than the value of the process
  # with s = (j, i), i.e. with x-coordinate equal j and y-coordinate equal i
  data.mat <- matrix(NA, nrow = ny, ncol = nx)
  for(i in 1:nrow(grid)){
    data.mat[grid[i,2], grid[i,1]] <- data[i]
  } # 1. column in grid is the x-coordinate, i.e. along the x-axis
    # 2. column in grid is the y-coordinate, i.e. along the y-axis
  
  # build the vectors for the MCD variogram estimator based on differences
  # (see section 3.2); Output: S-N, E-W, SW-NE, SE-NW
  vec.diff <- build.diff(data.mat, hmax)
  
  
  # build the vectors for the MCD variogram estimator based on the original data
  # (see section 3.1); Output: S-N, E-W, SW-NE, SE-NW
  vec.org <- build.org(data.mat, hmax)
  
  # if mod = TRUE; build the vetors for the modified MCD variogram estimators 
  # (see section 4.3);  Output: S-N, E-W,
  if(mod){
    mod.vec <- build.mod(data.mat, hmax, mx, my)
  }
  
  # build the differences needed for the Matheron and Genton variogram estimator;
  # Output: S-N, E-W, SW-NE, SE-NW
  mat.gen <- build.mat(data, grid, hmax)
  
  # calculate the variogram estimations based on the different estimators for 
  # all four directions
  diff <- sapply(1:4, function(x)  MCD.diff(vec.diff[[x]], reweighting = FALSE, missing = missing, det = det))
  diff.re <- sapply(1:4, function(x)  MCD.diff(vec.diff[[x]], reweighting = TRUE, missing = missing, det = det))
  org <- sapply(1:4, function(x) MCD.org(vec.org[[x]], reweighting = FALSE, missing = missing, det = det))
  org.re <- sapply(1:4, function(x) MCD.org(vec.org[[x]], reweighting = TRUE, missing = missing, det = det))
  Math <- Matheron(mat.gen) 
  Gen <- Genton(mat.gen)
  
  # if mod = TRUE; calculate also the variogram estimations based on the modified
  # MCD variogram estimators; these will be only calculated for the directions:
  # S-N and E-W!
  if(mod){  
    diff.mod <- list()
    diff.mod.re <- list()
    org.mod <- list()
    org.mod.re <- list()
    
    # MCD.diff.mod
    for(r in 1:2){
      mod.diff.null <- lapply(mod.vec[[r]][[2]], function(l) if(!is.null(l)){MCD.diff(l, reweighting = FALSE, missing = missing, det = det)})
      mod.diff <- c()
      for(l in 1:length(mod.diff.null)){
        if(!is.null(mod.diff.null[[l]])) mod.diff <- cbind(mod.diff, mod.diff.null[[l]])
      }
      diff.mod[[r]] <- rowMeans(mod.diff)
    }
    
    # MCD.diff.mod.re
    for(r in 1:2){
      mod.diff.null <- lapply(mod.vec[[r]][[2]], function(l) if(!is.null(l)){MCD.diff(l, reweighting = TRUE, missing = missing, det = det)})
      mod.diff <- c()
      for(l in 1:length(mod.diff.null)){
        if(!is.null(mod.diff.null[[l]])) mod.diff <- cbind(mod.diff, mod.diff.null[[l]])
      }
      diff.mod.re[[r]] <- rowMeans(mod.diff)
    }
    
    # Mod.org.mod
    for(r in 1:2){
      mod.org.null <- lapply(mod.vec[[r]][[1]], function(l) if(!is.null(l)){MCD.org(l, reweighting = FALSE, missing = missing, det = det)})
      mod.org <- c()
      for(l in 1:length(mod.org.null)){
        if(!is.null(mod.org.null[[l]])) mod.org <- cbind(mod.org, mod.org.null[[l]])
      }
      org.mod[[r]] <- rowMeans(mod.org)
    }
    
    # Mod.org.mod
    for(r in 1:2){
      mod.org.null <- lapply(mod.vec[[r]][[1]], function(l) if(!is.null(l)){MCD.org(l, reweighting = TRUE, missing = missing, det = det)})
      mod.org <- c()
      for(l in 1:length(mod.org.null)){
        if(!is.null(mod.org.null[[l]])) mod.org <- cbind(mod.org, mod.org.null[[l]])
      }
      org.mod.re[[r]] <- rowMeans(mod.org)
    }
  }

  # save the results
  resS.N <- cbind("MCD.diff" = diff[[1]], "MCD.diff.re" = diff.re[[1]], "MCD.org" = org[[1]],
                  "MCD.org.re" = org.re[[1]], "Matheron" = Math[[1]], "Genton" = Gen[[1]])
  resE.W <- cbind("MCD.diff" = diff[[2]], "MCD.diff.re" = diff.re[[2]], "MCD.org" = org[[2]],
                  "MCD.org.re" = org.re[[2]], "Matheron" = Math[[2]], "Genton" = Gen[[2]])
  resSW.NE <- cbind("MCD.diff" = diff[[3]], "MCD.diff.re" = diff.re[[3]], "MCD.org" = org[[3]],
                  "MCD.org.re" = org.re[[3]], "Matheron" = Math[[3]], "Genton" = Gen[[3]])
  resSE.NW <- cbind("MCD.diff" = diff[[4]], "MCD.diff.re" = diff.re[[4]], "MCD.org" = org[[4]],
                  "MCD.org.re" = org.re[[4]], "Matheron" = Math[[4]], "Genton" = Gen[[4]])
  
  if(mod){
    resS.N <- cbind(resS.N, "MCD.diff.mod" = diff.mod[[1]], "MCD.diff.mod.re" = diff.mod.re[[1]], "MCD.org.mod" = org.mod[[1]],
                    "MCD.org.mod.re" = org.mod.re[[1]])
    resE.W <- cbind(resE.W, "MCD.diff.mod" = diff.mod[[2]], "MCD.diff.mod.re" = diff.mod.re[[2]], "MCD.org.mod" = org.mod[[2]],
                    "MCD.org.mod.re" = org.mod.re[[2]])
  }
  
  
  res <- list(S.N = resS.N, E.W = resE.W, SW.NE = resSW.NE, SE.NW = resSE.NW)
  return(res)
}



