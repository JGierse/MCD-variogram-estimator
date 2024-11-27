################################################################################
## Function for building the vectors (W_1, ... W_(n_W)) needed for the        ##
## MCD variogram estimator based on differences for all directions            ##
## (see section 3.2)                                                          ##
################################################################################

## Input: data.mat: matrix - ny x nx dimensional; containing the data; for example
##                           the first row contains all data points z(s) with 
##                           s = (i,1) for i = 1, ..., nx and the first contains 
##                           all data points z(s) with s = (1, j) 
##                           for j = 1, ... , ny
##        hmax:     vector - 4 dimensional; containing for each direction the number
##                           of lags to be estimated (order: S-N, E-W, SE-NW, SW-NE)
## Output: res: list - 4 entries (for each direction one entry; order: S-N, E-W, 
##                     SE-NW, SW-NE)); each entry a ((hmax + 1) x nv) dimensional
##                     matrix including the vectors for the corresponding direction


build.diff <- function(data.mat, 
                       hmax){
  
  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in 
                       # east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in 
                       # north-south direction
  
  
  ## S-N: build the vectors for the estimation in north-south direction
  # build a matrix for storing the vectors: in each column one vector
  vecS.N <- matrix(NA, nrow = hmax[1], ncol = nx*(ny - hmax[1])) 
  ind.vec <- 1 # counts the vectors already build
  
  for(i in 1:nx){ # for each column, i.e. for each x-coordinate, 
                  # build the following vectors (see also explanations in 
                  # section 3.2)
    for(j in ny:(hmax[1] + 1)){
      vecS.N[,ind.vec] <- data.mat[j,i] - data.mat[(j-1):(j-hmax[1]),i]
      ind.vec <- ind.vec + 1
    }
  }
  
  ## E-W: build the vectors for the estimation in east-west direction
  # build a matrix for storing the vectors: in each column one vector
  vecE.W <- matrix(NA, nrow = hmax[2], ncol = ny*(nx - hmax[2]))
  ind.vec <- 1 # counts the vectors already build
  
  for(j in 1:ny){ # for each row, i.e. for each y-coordinate, 
                  # build the following vectors (see also explanations in 
                  # section 3.2)
    for(i in 1:(nx - hmax[2])){
      vecE.W[,ind.vec] <- data.mat[j,i] - data.mat[j,(i+1):(i+hmax[2])]
      ind.vec <- ind.vec + 1
    }
  }
  
  ## SW-NE: build the vectors for the estimation in southwest-northeast direction
  # build a matrix for storing the vectors: in each column one vector
  vecSW.NE <- matrix(NA, nrow = hmax[3], ncol = (nx - hmax[3])*(ny-hmax[3]))
  ind.vec <- 1 # counts the vectors already build
  
  for(i in 1:(nx - hmax[3])){ # i: x-coordinate (east-west)
    for(j in ny:(hmax[3]+1)){ # j: y-coordinate (north-south)
                              # (see also explanations in section 3.2)
      vecSW.NE[,ind.vec] <- data.mat[j, i ] - diag(data.mat[(j-1):(j-hmax[3]), (i+1):(i+hmax[3])])
      ind.vec <- ind.vec + 1
    }
  }
  
  ## SE-NW: build the vectors for the estimation in southeast-northwest direction
  # build a matrix for storing the vectors: in each column one vector
  vecSE.NW <- matrix(NA, nrow = hmax[4], ncol = (nx - hmax[4])*(ny-hmax[4]))
  ind.vec <- 1 # counts the vectors already build
  
  for(i in 1:(nx - hmax[4])){    # i: x-coordinate (east-west)
    for(j in (ny -  hmax[4]):1){ # j: y-coordinate (north-south)
                                 # (see also explanations in section 3.2)
      vecSE.NW[,ind.vec] <- data.mat[j, i ] - diag(data.mat[(j+1):(j+hmax[4]), (i+1):(i+hmax[4])])
      ind.vec <- ind.vec + 1
    }
  }
  
  res <- list(vecS.N = vecS.N, vecE.W = vecE.W, vecSW.NE = vecSW.NE, vecSE.NW = vecSE.NW)
  return(res)
}


################################################################################
## Function for building the vectors (V_1, ... V_(n_V)) needed for the        ##
## MCD variogram estimator based on the original data for all directions      ##
## (see section 3.1)                                                          ##
################################################################################

## Input: data.mat: matrix - ny x nx dimensional; containing the data; for example
##                           the first row contains all data points z(s) with 
##                           s = (i,1) for i = 1, ..., nx and the first contains 
##                           all data points z(s) with s = (1, j) 
##                           for j = 1, ... , ny
##        hmax:     vector - 4 dimensional; containing for each direction the number
##                           of lags to be estimated (order: S-N, E-W, SE-NW, SW-NE)
## Output: res: list - 4 entries (for each direction one entry; order: S-N, E-W, 
##                     SE-NW, SW-NE)); each entry a ((hmax + 1) x nv) dimensional
##                     matrix including the vectors for the corresponding direction

build.org <- function(data.mat, 
                      hmax){
  
  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in 
                       # east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in 
                       # north-south direction
  
  
  ## S-N: build the vectors for the estimation in north-south direction
  # build a matrix for storing the vectors: in each column one vector
  vecS.N <- matrix(NA, nrow = hmax[1] + 1, ncol = nx*(ny - hmax[1]))
  ind.vec <- 1  # counts the vectors already build
  
  for(i in 1:nx){ # for each column, i.e. for each x-coordinate, 
                  # build the following vectors (see also explanations in 
                  # section 3.1)
    for(j in ny:(hmax[1] + 1)){
      vecS.N[,ind.vec] <- data.mat[j:(j-hmax[1]), i]
      ind.vec <- ind.vec + 1
    }
  }
  
  ## E-W: build the vectors for the estimation in east-west direction
  # build a matrix for storing the vectors: in each column one vector
  vecE.W <- matrix(NA, nrow = hmax[2] + 1, ncol = ny*(nx - hmax[2]))
  ind.vec <- 1 # counts the vectors already build
  
  for(j in 1:ny){  # for each row, i.e. for each y-coordinate, 
                   # build the following vectors (see also explanations in 
                   # section 3.1)
    for(i in 1:(nx - hmax[2])){
      vecE.W[,ind.vec] <- data.mat[j,i:(i+hmax[2])]
      ind.vec <- ind.vec + 1
    }
  }

  ## SW-NE: build the vectors for the estimation in southwest-northeast direction
  # build a matrix for storing the vectors: in each column one vector
  vecSW.NE <- matrix(NA, nrow = hmax[3] + 1, ncol = (nx - hmax[3])*(ny-hmax[3]))
  ind.vec <- 1  # counts the vectors already build
  
  for(i in 1:(nx - hmax[3])){ # i: x-coordinate (east-west)
    for(j in ny:(hmax[3]+1)){ # j: y-coordinate (north-south)
                              # (see also explanations in section 3.1)
      vecSW.NE[,ind.vec] <- diag(data.mat[j:(j-hmax[3]), i:(i+hmax[3])])
      ind.vec <- ind.vec + 1
    }
  }
  
  ## SE-NW: build the vectors for the estimation in southeast-northwest direction
  # build a matrix for storing the vectors: in each column one vector
  vecSE.NW <- matrix(NA, nrow = hmax[4] + 1, ncol = (nx - hmax[4])*(ny-hmax[4]))
  ind.vec <- 1   # counts the vectors already build
  
  for(i in 1:(nx - hmax[4])){    # i: x-coordinate (east-west)
    for(j in (ny -  hmax[4]):1){ # j: y-coordinate (north-south)
                                 # (see also explanations in section 3.1)
      vecSE.NW[,ind.vec] <- diag(data.mat[j:(j+hmax[4]), i:(i+hmax[4])])
      ind.vec <- ind.vec + 1
    }
  }
  
  res <- list(vecS.N = vecS.N, vecE.W = vecE.W, vecSW.NE = vecSW.NE, vecSE.NW = vecSE.NW)
  return(res)
}


################################################################################
## Function for building all differences of the process for all directions    ##
## and all lags; these differences are needet for the estimation of the       ## 
## Matheron variogram estimator and for the Genton variogram estimator        ##
################################################################################

## Input: data: vector - nx X ny dimensional; containing the data 
##        grid: matrix - (nx X ny) rows; 2 columns; first column contains the 
##                       x-coordinate of the observations and the second column
##                       contains the y-coordinate of the observations
##        hmax: vector - 4 dimensional; containing for each direction the number
##                       of lags to be estimated (order: S-N, E-W, SE-NW, SW-NE)
## Output: res: list - 4 entries (for each direction one entry; order: S-N, E-W, 
##                     SE-NW, SW-NE)); each entry 


build.mat <- function(data, 
                      grid, 
                      hmax){
  
  # build two equal datasets; for each data point a row with the following columns
  # first: an ID; second: x-coordinate; third: y-coordinate; fourth: data value
  grid.1 <- cbind(1:nrow(grid), grid, data)
  colnames(grid.1) <- c("ID.1", "x.1", "y.1", "data.1")  
  grid.2 <- cbind(1:nrow(grid), grid, data)
  colnames(grid.2) <- c("ID.2", "x.2", "y.2", "data.2")
  
  # combine the two equal data set to calculate the distance of each point of the
  # grid to all other points of the grid
  grid.ind <- expand.grid(1:nrow(grid), 1:nrow(grid))
  grid.comp <- as.data.frame(cbind(grid.1[grid.ind[,1],], grid.2[grid.ind[,2],]))
  
  # delete double rows: e.g. ID.1 == 1, ID.2 == 2 has the same distance between 
  # the data and the same difference of the data points as ID.1 == 2 and ID.2 == 1
  grid.comp <- subset(grid.comp, grid.comp$ID.1 < grid.comp$ID.2)
  
  
  # Calculate the distance between the two data points in the different directions
  grid.comp$lagS.N <- ifelse(grid.comp$x.1 == grid.comp$x.2, abs(grid.comp$y.1 - grid.comp$y.2), NA)
  grid.comp$lagE.W <- ifelse(grid.comp$y.1 == grid.comp$y.2, abs(grid.comp$x.1 - grid.comp$x.2), NA)
  grid.comp$lagSE.NW <- ifelse((grid.comp$x.1 - grid.comp$x.2) == -(grid.comp$y.1 - grid.comp$y.2), 
                               sqrt((grid.comp$x.1-grid.comp$x.2)^2 + (grid.comp$y.1 - grid.comp$y.2)^2), NA)
  grid.comp$lagSW.NE <- ifelse((grid.comp$x.1 - grid.comp$x.2) == (grid.comp$y.1 - grid.comp$y.2), 
                               sqrt((grid.comp$x.1-grid.comp$x.2)^2 + (grid.comp$y.1 - grid.comp$y.2)^2), NA)
  
  # delete all unnecessary rows: rows with distance greater than the 
  # maxium lags which should be considered in the corresponding direction
  hmaxSW.NE <- unique(grid.comp$lagSW.NE)[hmax[3] + 1] # maximum distance in the SW.NE direction
  hmaxSE.NW <- unique(grid.comp$lagSE.NW)[hmax[4] + 1] # maximum distance in the SE.NW direction
  grid.comp <- subset(grid.comp, (grid.comp$lagS.N <= hmax[1])|(grid.comp$lagE.W <= hmax[2])|(grid.comp$lagSW.NE <= hmaxSW.NE)|(grid.comp$lagSE.NW <= hmaxSE.NW))
  
  # build the differences of the two data points
  grid.comp$diff <- grid.comp$data.1 - grid.comp$data.2
  
  # split into four different data sets (for each direction one dataset)
  # delete unnecessary columns
  S.N <- subset(grid.comp, !is.na(grid.comp$lagS.N))  
  S.N <- subset(S.N, select = -c(lagE.W, lagSW.NE, lagSE.NW, ID.1, ID.2))
  E.W <- subset(grid.comp, !is.na(grid.comp$lagE.W))  
  E.W <- subset(E.W, select = -c(lagS.N, lagSW.NE, lagSE.NW, ID.1, ID.2))
  SE.NW <- subset(grid.comp, !is.na(grid.comp$lagSE.NW))  
  SE.NW <- subset(SE.NW, select = -c(lagE.W, lagSW.NE, lagS.N, ID.1, ID.2))
  SW.NE <- subset(grid.comp, !is.na(grid.comp$lagSW.NE))  
  SW.NE <- subset(SW.NE, select = -c(lagE.W, lagS.N, lagSE.NW, ID.1, ID.2))
  
  
  res <- list(S.N = S.N, E.W = E.W, SE.NW = SE.NW, SW.NE = SW.NE)
  return(res)
}



################################################################################
## Function for building the vectors for the modified MCD variogram           ##
## estimators based on differences and based on the original data             ##
## as described in section 4.3 (for the S-N direction and the E-W direction)  ##
################################################################################

## Input: data.mat: matrix  - ny x nx dimensional; containing the data; for example
##                            the first row contains all data points z(s) with 
##                            s = (i,1) for i = 1, ..., nx and the first contains 
##                            all data points z(s) with s = (1, j) 
##                            for j = 1, ... , ny
##        hmax:     vector  - 2 dimensional; containing for each direction the number
##                            of lags to be estimated (order: S-N, E-W)
##        mx:       numeric - dependency in the E-W direction (along the x-axis);
##                            for the modified version mx dependency in the E-W
##                            direction is assumed (see section 4.3)
##        my:       numeric - dependency in the N-S direction (along the y-axis);
##                            for the modified version my dependency in the N-S
##                            direction is assumed (see section 4.3)
## Output: res: list - 2 entries (for each direction one entry; order: S-N, E-W)); 
##                     each entry is a list with 2 elements; !!!

build.mod <- function(data.mat, 
                      hmax, 
                      mx, 
                      my){
  
  nx <- ncol(data.mat) # number of x-coordinates, i. e. number of coordinates in 
                       # east-west direction
  ny <- nrow(data.mat) # number of y-coordinates, i. e. number of coordinates in 
                       # north-south direction
  

  ## S-N: build the vectors for the estimation in north-south direction
  vec.listS.N <- list() # list for the vectors based on the original data, i.e. 
                        # vectors for MCD.org.mod (see section 4.3)
  diff.listS.N <- list() # list for the vectors based on differences of the data,
                         # i.e. vectors for MCD.diff.mod (see section 4.3)
  
  # as described in section 4.3, it exists different possibilities to choose 
  # independent vectors for the modified MCD estimators. We build here all possible
  # sets of independent vectors
  ind <- 1 # counts the vectors already build
  for(i.S in 1:(mx+1)){ 
    vecs <- c()
    diffs <- c()
    cols <- seq(i.S, nx, mx + 1)
    
    for(j.S in 1:min((hmax[1] + my + 1),  ny - hmax[1])){
      rows <- seq(j.S, ny - hmax[1], hmax[1] + my + 1)
      for(i.V in cols){
        for(j.V in rows){
          vecs <- cbind(vecs, data.mat[j.V:(j.V + hmax[1]), i.V])
          diffs <- cbind(diffs, data.mat[j.V, i.V] - data.mat[(j.V+1):(j.V+hmax[1]), i.V])
        }
      }
     vec.listS.N[[ind]] <- vecs
     diff.listS.N[[ind]] <- diffs
     ind <- ind + 1
    }
  }
  # delete all sets of independent vectors that have too few vectors for the 
  # estimation
  vec.listS.N <- lapply(vec.listS.N, function(l){if(ncol(l) > 2* hmax[1]) l}) 
  diff.listS.N <- lapply(diff.listS.N, function(l){if(ncol(l) > 2* hmax[1] ) l})
  
  
  
  ## E-W: build the vectors for the estimation in east-west direction
  vec.listE.W <- list() # list for the vectors based on the original data, i.e. 
                        # vectors for MCD.org.mod (see section 4.3)
  diff.listE.W <- list() # list for the vectors based on differences of the data,
                         # i.e. vectors for MCD.diff.mod (see section 4.3)
  
  # as described in section 4.3, it exists different possibilities to choose 
  # independent vectors for the modified MCD estimators. We build here all possible
  # sets of independent vectors
  ind <- 1 # counts the vectors already build
  for(j.S in 1:(my+1)){
    vecs <- c()
    diffs <- c()
    rows <- seq(j.S, ny, my + 1)
    
    for(i.S in 1:min((hmax[2] + mx + 1),  nx - hmax[2])){
      cols <- seq(i.S, nx - hmax[2], hmax[2] + mx + 1)
      for(j.V in rows){
        for(i.V in cols){
          vecs <- cbind(vecs, data.mat[j.V, i.V:(i.V+hmax[2])])
          diffs <- cbind(diffs, data.mat[j.V, i.V] - data.mat[j.V, (i.V+1):(i.V+hmax[2])])
        }
      }
      vec.listE.W[[ind]] <- vecs
      diff.listE.W[[ind]] <- diffs
      ind <- ind + 1
    }
  }
  # delete all sets of independent vectors that have too few vectors for the 
  # estimation
  vec.listE.W <- lapply(vec.listE.W, function(l){if(ncol(l) > 2* hmax[1]) l})
  diff.listE.W <- lapply(diff.listE.W, function(l){if(ncol(l) > 2* hmax[1]) l})
  
  
  res <- list(S.N = list(vec.listS.N = vec.listS.N, diff.listS.N = diff.listS.N),
              E.W = list(vec.listE.W = vec.listE.W, diff.listE.W = diff.listE.W)) 
  
  return(res)
}

