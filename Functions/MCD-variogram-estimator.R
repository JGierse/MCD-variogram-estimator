################################################################################
## Function for the calculation of the MCD variogram estimator                ##
## based on differences (section 3.2)                                         ##
################################################################################

## Input: mat:          matrix - (hmax x nv) dimensional; including the differences
##                               for one direction
##        reweighting: logical - whether the reweighted MCD should be calculated;
##                               default: FALSE
##        cp:          numeric - correction factor; default: 1
##        missing:     logical - TRUE if the data contains missing values; the 
##                               columns with missing values will be delete;
##                               default: FALSE
##        det:         logical - TRUE if the deterministic MCD should be used
##                               instead of FASTMCD (see also explanations in 
##                               covMcd of the package robustbase);
##                               default: FALSE
## Output: estimate:   vector  - hmax dimensional; containing the variogram 
##                               estiamtion

MCD.diff <- function(mat, 
                     reweighting = FALSE, 
                     cp = 1, 
                     missing = FALSE, 
                     det = FALSE){
  
  ## covMcd of the package robustbase is used here
  library(robustbase)
  

  hmax <- nrow(mat)
  estimate <- rep(0, hmax)
  
  ## if missing = TRUE; delete all vectors (columns) with missing errors
  if(missing){
    mat.miss <- c()
    for(i in 1:ncol(mat)){
      if(!any(is.na(mat[,i]))){
        mat.miss <- cbind(mat.miss, mat[,i])
      }
    }
    mat <- mat.miss
  }
  mat <- t(mat)
  
  ## estimation of the variogram
  if(det == TRUE){ # use Deterministic MCD Algorithm of Hubert et al. (2012)
                   # (see documentation of covMcd in robustbase)
    
    if(reweighting == FALSE){ # MCD.diff.re
      covest <- covMcd(mat, raw.only = TRUE, use.correction = FALSE, nsamp = "deterministic")
      estimate <- cp*(diag(covest$cov))
    } else{ # MCD.diff
      covest <- covMcd(mat, raw.only = FALSE, use.correction = FALSE, nsamp = "deterministic")
      estimate <- cp*(diag(covest$cov))
    }
  } 
  if(det == FALSE){  # use FastMCD Algorithm of Rousseeuw and van Driessen (1999)
                     # (see documentation of covMcd in robustbase)
    
    if(reweighting == FALSE){ # MCD.diff.re
      covest <- covMcd(mat, raw.only = TRUE, use.correction = FALSE)
      estimate <- cp*(diag(covest$cov))
    } else{ # MCD.diff
      covest <- covMcd(mat, raw.only = FALSE, use.correction = FALSE)
      estimate <- cp*(diag(covest$cov))
    }
  }
  
  return(estimate)
}

################################################################################
## Function for the calculation of the MCD variogram estimator                ##
## based on the original data (section 3.1)                                   ##
################################################################################

## Input: mat:          matrix - ((hmax + 1) x nv) dimensional; including the 
##                               vectors for one direction
##        reweighting: logical - whether the reweighted MCD should be calculated;
##                               default: FALSE
##        cp:          numeric - correction factor; default: 1
##        missing:     logical - TRUE if the data contains missing values; the 
##                               columns with missing values will be delete;
##                               default: FALSE
##        det:         logical - TRUE if the deterministic MCD should be used
##                               instead of FASTMCD (see also explanations in 
##                               covMcd of the package robustbase);
##                               default: FALSE
## Output: estimate:   vector  - hmax dimensional; containing the variogram 
##                               estiamtion

MCD.org <- function(mat, 
                    reweighting = FALSE, 
                    cp = 1, 
                    missing = FALSE, 
                    det = FALSE){
  
  ## covMcd of the package robustbase is used here
  library(robustbase)
  
  hmax <- nrow(mat) - 1
  estimate <- rep(0, hmax)
  
  ## if missing = TRUE; delete all vectors (columns) with missing errors
  if(missing){
    mat.miss <- c()
    for(i in 1:ncol(mat)){
      if(!any(is.na(mat[,i]))){
        mat.miss <- cbind(mat.miss, mat[,i])
      }
    }
    mat <- mat.miss
  }
  mat <- t(mat)
  
  ## estimation of the Covarainzmatrix of the vectors
  if(det == TRUE){ # use Deterministic MCD Algorithm of Hubert et al. (2012)
    # (see documentation of covMcd in robustbase)
    
    if(reweighting == FALSE){ # raw
      est1 <- covMcd(mat, raw.only = TRUE, use.correction = FALSE, nsamp = "deterministic")$cov
    } else { # reweighted
      est1 <- covMcd(mat, raw.only = FALSE, use.correction = FALSE,  nsamp = "deterministic")$cov
    }
  }
  if(det == FALSE){ # use FastMCD Algorithm of Rousseeuw and van Driessen (1999)
    # (see documentation of covMcd in robustbase)
    
    if(reweighting == FALSE){ # raw 
      est1 <- covMcd(mat, raw.only = TRUE, use.correction = FALSE)$cov
    } else { # reweighted
      est1 <- covMcd(mat, raw.only = FALSE, use.correction = FALSE)$cov
    }
  }
  
  # Calculate the variogram estimation based on the formula (1) (section 3.1)
  X <- est1
  gamma.0 <- mean(diag(X)) # average the different variance estimations (C(0))
                           # contained on the main diagonal
  Xr <- Xl <- X 
  for(j in 1:(hmax - 1)){ # average the different estimations of c(1), ... , C(hmax-1)
                          # see explanations in section 3.1
    Xr <- Xr[-1,-(hmax + 2 -j)] 
    Xl <- Xl[-(hmax + 2 -j),-1]
    cov <- mean(c(diag(Xr), diag(Xl)))
    estimate[j] <- 2 * gamma.0 - 2 * cov # claculate Variogram for lags 1, ..., hmax
  }
  estimate[hmax] <- (2 * gamma.0 - 2 * mean(Xr[1,2], Xr[2,1])) # calculate variogram for lag hmax
  
  return(cp * estimate)
}

################################################################################
## Function for calculation of the Matheron variogram estimator               ##
## (Matheron, 1962)                                                           ##
################################################################################

## Input: datas:     list - 4 entries (for each direction one entry; 
##                          order: S-N, E-W, SE-NW, SW-NE)); 
##                          each entry is an array with 8 columns containing
##                          all difference for all directions; the columns are:
##                          x.1 (x-coordinate of the first data point (s1)), 
##                          y.1 (y-coordinate of the first data point (s1)),
##                          x.2 (x-coordinate of the second data point (s2)),
##                          y.2 (y-coordinate of the second data point (s2)),
##                          data.1 (observed value of the first data point (z(s1))),
##                          data.2 (observed value of the second data point (z(s2))),
##                          lag (distance between s1 and s2),
##                          diff (difference between z(s1) and z(s2))
## Output: estimate: list - 4 entries (for each direction one entry; 
##                          order: S-N, E-W, SE-NW, SW-NE)); 
##                          each entry a hmax dimensional vector containing the
##                          Matheron variogram estimation in the corresponding 
##                          direction


Matheron <- function(datas){
  
  estimate <- list()
  for(i in 1:4){ # for each direction (i = 1: S-N; i = 2: E-W; i = 3: SE-NW; i = 4: SW-NE)
    data.split <- split(datas[[i]], datas[[i]][,7]) # split for the different lags in the i-th direction
    
    # calculate the variogram direction for all lags in direction i
    estimate[[i]] <- sapply(1:length(data.split), function(x) mean(data.split[[x]]$diff^2, na.rm = TRUE))
  }
  
  return(estimate)
}

################################################################################
## Function for calculation of the Genton variogram estimator                 ##
## (Genton, 1998a)                                                            ##
################################################################################

## Input: datas:     list - 4 entries (for each direction one entry; 
##                          order: S-N, E-W, SE-NW, SW-NE)); 
##                          each entry is an array with 8 columns containing
##                          all difference for all directions; the columns are:
##                          x.1 (x-coordinate of the first data point (s1)), 
##                          y.1 (y-coordinate of the first data point (s1)),
##                          x.2 (x-coordinate of the second data point (s2)),
##                          y.2 (y-coordinate of the second data point (s2)),
##                          data.1 (observed value of the first data point (z(s1))),
##                          data.2 (observed value of the second data point (z(s2))),
##                          lag (distance between s1 and s2),
##                          diff (difference between z(s1) and z(s2))
## Output: estimate: list - 4 entries (for each direction one entry; 
##                          order: S-N, E-W, SE-NW, SW-NE)); 
##                          each entry a hmax dimensional vector containing the
##                          Genton variogram estimation in the corresponding 
##                          direction


Genton <- function(datas){
  
  ## Qn of the package robustbase is used here
  library(robustbase)
  
  estimate <- list()
  for(i in 1:4){ # for each direction (i = 1: S-N; i = 2: E-W; i = 3: SE-NW; i = 4: SW-NE)
    data.split <- split(datas[[i]], datas[[i]][,7]) # split for the different lags in the i-th direction
    
    # calculate the variogram direction for all lags in direction i (using function Qn of package robustbase)
    estimate[[i]] <- sapply(1:length(data.split), function(x) Qn(na.omit(data.split[[x]]$diff))^2)
  }
  
  return(estimate)
}
