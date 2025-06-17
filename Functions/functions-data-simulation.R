################################################################################
## Function for the simulation of the data, i.e. function to simulate the     ##
## Gaussian random fields                                                     ##
################################################################################

## Input: gridsize:        vector           - 2-dim vector contain the gridsize in x-direction (nx) and
##                                            the gridsize in y-direction (ny); default: nx = ny = 15
##        dist.outlier:    function         - function to draw random numbers from the distribution of the 
##                                            outliers; default: rnorm
##        variogram:       RMmodelgenerator - semivariogram model, notation of package RandomFields; 
##        param.variogram: vector           - including the sill (var) and the range (scale) of
##                                            the variogram; default: var = 1, scale = 5
##        aniso.param:     vector           - parameters of the geometric anisotropy if the
##                                            process is anisotropic, notation of RandomFields
##                                            if the process should be isotrop it has to be null;
##                                            default: NULL -> isotropy
##        out.type:        charackter       - type of the outliers, possibilities: "block", "isolated";
##                                            default: NULL
##        amount:          numeric          - amount of outlier; default: NULL 
##        param.outlier:   vector           - parameters of the outlier distribution and the
##                                            distribution; default: NULL
##        mixture:         logical          - whether the outlier should be genereted as mixture 
##                                            model or as additive outliers; default: TRUE -> mixture
##        nugget:          numeric          - if not NULL, value of the nugget effekt of the semivariogram;
##                                            default: NULL
##        n.it:            numeric          - number of iterations; default: 1000
## Output: res: list - first element a matrix: in each column is one data set 
##                   - second element: a matrix with the coordinates


sim.data <- function(gridsize = c("nx" = 15, "ny" = 15),
                     dist.outlier = rnorm, 
                     variogram = "spherical", 
                     param.variogram = c("var" = 1, "scale" = 5),
                     aniso.param = NULL, 
                     out.type = NULL, 
                     amount = NULL,
                     param.outlier = NULL, 
                     mixture = TRUE, 
                     nugget = NULL,
                     n.it = 1000){
  
  # needed for simulation
  library(geoR)
  
  # simulate data without outliers
  if(is.null(aniso.param)){ # isotropy
    if(is.null(nugget)){ # without nugget effect
      data <- grf(grid = expand.grid(1:gridsize["nx"], 1:gridsize["ny"]),
                  cov.model = variogram, cov.pars = param.variogram, 
                  message = FALSE, nsim = n.it)
    } else{ # nugget effect
      data <- grf(grid = expand.grid(1:gridsize["nx"], 1:gridsize["ny"]),
                  cov.model = variogram, cov.pars = param.variogram, nugget = nugget,
                  message = FALSE, nsim = n.it)
    }
  } else{ # anisotropy
    if(is.null(nugget)){ # without nugget effect
      data <- grf(grid = expand.grid(1:gridsize["nx"], 1:gridsize["ny"]),
                  cov.model = variogram, cov.pars = param.variogram, aniso.pars = aniso.param,
                  message = FALSE, nsim = n.it)
    } else{ # nugget effect
      data <- grf(grid = expand.grid(1:gridsize["nx"], 1:gridsize["ny"]),
                  cov.model = variogram, cov.pars = param.variogram, nugget = nugget, aniso.pars = aniso.param,
                  message = FALSE, nsim = n.it)
    }
  }
  
  grid.coord <- expand.grid(1:gridsize["nx"], 1:gridsize["ny"]) # save the grid of the data
  data <- data$data # save the simulated data (for each iteration one column)
  
  # if desired include outliers
  if(!is.null(out.type)){
    for(n in 1:n.it){ # for each iteration 
      
      # generate the coordinates of the outliers
      coords.outlier <- gen.coords(type = out.type, amount = amount, gridsize = gridsize)
      row.out <- apply(coords.outlier, 1, function(i) which(as.numeric(i[1]) == grid.coord[,1] & as.numeric(i[2]) == grid.coord[,2]))
      
      # simulate the values of the outliers 
      if(length(param.outlier) == 2){
        cont <- dist.outlier(nrow(coords.outlier), param.outlier[1], param.outlier[2])
      } else{
        cont <- dist.outlier(nrow(coords.outlier), param.outlier[1])
      }
      
      # add the outliers to the simulated data
      if(mixture){ # replaced the original data by the outliers (mixture model)
        data[row.out, n] <- cont
      } else{ # add the outliers to the original data (additive model)
        data[row.out, n] <- data[row.out, n] + cont
      }
      
    }
  }
  
  return(res <- list(data = data, grid = grid.coord))
}


################################################################################
## Function to determine the coordinates of the outlier                       ##
################################################################################

## input: 
##        type:     charackter - type of the outliers, possibilities: "block", "isolated";
##                               default: "block"
##        amount:   numeric    - amount of outlier 
##        gridsize: vector     - 2-dim vector contain the gridsize in x-direction (nx) and
##                               the gridsize in y-direction (ny)
## output: coords: 


gen.coords <- function(type = "block", 
                       amount, 
                       gridsize){
  
  # determine the grid
  grid.coord <- expand.grid(1:gridsize["nx"], 1:gridsize["ny"])
  # determine the gridsize
  n <- nrow(grid.coord) 
  
  ## isolated outliers: draw the coordinates at random (see also section 4.5)
  if(type == "isolated"){
    row.ind <- sample(1:n, ceiling(n*amount))
    coords <- grid.coord[row.ind, ]
  }
  
  ## block outliers: build a nearly quadratic block of coordinates (see also section 4.4)
  if(type == "block"){
    
    # amount of outliers 
    n.outlier <- ceiling(n * amount) 
    
    # determine the block size of the largest possible square block
    dim.block <- c(floor(sqrt(n.outlier)), floor(sqrt(n.outlier)))
    
    # change block dimension, if the grid in direction is smaller than the 
    # largest possible square block
    if(gridsize["nx"] < dim.block[1]){
      dim.block[2] <- dim.block[2] + (dim.block[1] - gridsize["nx"])
      dim.block[1] <- gridsize["nx"]
    }
    if(gridsize["ny"] < dim.block[2]){
      dim.block[1] <- dim.block[1] + (dim.block[2] - gridsize["ny"])
      dim.block[2] <- gridsize["ny"]
    }
    
    # determine how much smaller the largest possible square block is compared to 
    # the required number of outliers 
    rest <- n.outlier - dim.block[1] * dim.block[2]
    
    # Determine a random grid point (s0) around which the block is to be constructed
    # (see also section 4.4)
    start.point <- grid.coord[sample(1:n, 1),]
    
    # 1. Determine the x-coordinates of the largest possible square block  
    # with s0 as centred as possible
    dims <- sample(c(floor((dim.block[1]-1)/2), ceiling((dim.block[1]-1)/2)), 2)
    x.coords <- as.numeric((start.point[1] - dims[1])):as.numeric((start.point[1] + dims[2]))
    
    x.true.u <- sum(x.coords < 1)
    x.true.o <- sum(x.coords > gridsize["nx"])
    
    if(x.true.u > 0){
      x.coords <- x.coords + x.true.u
    }
    if(x.true.o > 0){
      x.coords <- x.coords - x.true.o
    }
    
    
    # 2. Determine the y-coordinates of the largest possible square block  
    # with s0 as centred as possible
    dims <- sample(c(floor((dim.block[2]-1)/2), ceiling((dim.block[2]-1)/2)), 2)
    y.coords <- as.numeric((start.point[2] - dims[1])):as.numeric((start.point[2] + dims[2]))
    
    y.true.u <- sum(y.coords < 1)
    y.true.o <- sum(y.coords > gridsize["ny"])
    
    if(y.true.u > 0){
      y.coords <- y.coords + y.true.u
    }
    if(y.true.o > 0){
      y.coords <- y.coords - y.true.o
    }
    
    # Assemble the coordinates of the largest possible square block around s0 
    # from the two components (1. & 2.) 
    coords <- expand.grid(x.coords, y.coords)
    
    
    # Add the remaining grid points so that the outlier block contains exactly 
    # n.outlier points; the block should be as square as possible
    if(rest > 0){
      dims.rest <- sample(c(floor(rest/2), ceiling(rest/2)), 2)
      
      if(min(x.coords) == 1 & max(x.coords) == gridsize["nx"]){
        dims.rest[2] <- rest
        dims.rest[1] <- 0
      }
      
      if(min(y.coords) == 1 & max(y.coords) == gridsize["ny"]){
        dims.rest[1] <- rest
        dims.rest[2] <- 0
      }
      
      
      if(!any(dims.rest == 0)){
        # x-direction
        if(min(x.coords) == 1){
          x.rest <- max(x.coords) + 1
        }
        if(max(x.coords) == gridsize["nx"]){
          x.rest <- min(x.coords) - 1
        }
        if(min(x.coords) != 1 & max(x.coords) != gridsize["nx"]){
          x.rest <- sample(c((max(x.coords) + 1), (min(x.coords) - 1)), 1)
        }
        if(dims.rest[1] == dim.block[2]){
          y.rest <- y.coords 
        } else {
          r <- dim.block[2] - dims.rest[1]
          help <- sample(c(0,1), 1)
          
          if(help == 1){
            y.rest <- y.coords[-(1:r)]
          } else{
            y.rest <- y.coords[dim.block[2]:1]
            y.rest <- y.rest[-(1:r)]
          }
        }
        coords <- rbind(coords, cbind(Var1  = x.rest, Var2 = y.rest))
        
        # y-direction
        if(min(y.coords) == 1){
          y.rest <- max(y.coords) + 1
        }
        if(max(y.coords) == gridsize["ny"]){
          y.rest <- min(y.coords) - 1
        }
        if(min(y.coords) != 1 & max(y.coords) != gridsize["ny"]){
          y.rest <- sample(c((max(y.coords) + 1), (min(y.coords) - 1)), 1)
        }
        if(dims.rest[2] == dim.block[1]){
          x.rest <- x.coords 
        } else {
          r <- dim.block[1] - dims.rest[2]
          help <- sample(c(0,1), 1)
          
          if(help == 1){
            x.rest <- x.coords[-(1:r)]
          } else{
            x.rest <- x.coords[dim.block[1]:1]
            x.rest <- x.rest[-(1:r)]
          }
        }
        coords <- rbind(coords, cbind(Var1  = x.rest, Var2 = y.rest))
      } else{
        dims.rest2 <- c(ceiling(rest/2), floor(rest/2))
        
        if(dims.rest[1] == 0){
          if(min(y.coords) == 1){
            y.rest <- c(max(y.coords) + 1, max(y.coords) + 2)
          }
          if(max(y.coords) == gridsize["ny"]){
            y.rest <- c(min(y.coords) - 1, min(y.coords) - 2)
          }
          if(min(y.coords) != 1 & max(y.coords) != gridsize["ny"]){
            y.rest <- c((max(y.coords) + 1), (min(y.coords) - 1))
          }
          
          if(dims.rest2[1] == dim.block[1]){
            x.rest <- x.coords 
          } else {
            r1 <- dim.block[1] - dims.rest2[1]
            r2 <- dim.block[1] - dims.rest2[2]
            help <- sample(c(0,1), 1)
            
            if(help == 1){
              x.rest <- list(x.coords[-(1:r1)], x.coords[-(1:r2)])
            } else{
              x.rev <- x.coords[dim.block[1]:1]
              x.rest <- list(x.rev[-(1:r1)], x.rev[-(1:r2)])
            }
          }
          
          coords <- rbind(coords, cbind(Var1  = unlist(x.rest[1]), Var2 = y.rest[1]), cbind(Var1 = unlist(x.rest[2]), Var2 = y.rest[2]))
          
        }
        
        if(dims.rest[2] == 0){
          if(min(x.coords) == 1){
            x.rest <- c(max(x.coords) + 1, max(x.coords) + 2)
          }
          if(max(x.coords) == gridsize["ny"]){
            x.rest <- c(min(x.coords) - 1, min(x.coords) - 2)
          }
          if(min(x.coords) != 1 & max(x.coords) != gridsize["ny"]){
            x.rest <- c((max(x.coords) + 1), (min(x.coords) - 1))
          }
          
          if(dims.rest2[1] == dim.block[2]){
            y.rest <- y.coords 
          } else {
            r1 <- dim.block[2] - dims.rest2[1]
            r2 <- dim.block[2] - dims.rest2[2]
            help <- sample(c(0,1), 1)
            
            if(help == 1){
              y.rest <- list(y.coords[-(1:r1)], y.coords[-(1:r2)])
            } else{
              y.rev <- y.coords[dim.block[1]:1]
              y.rest <- list(y.rev[-(1:r1)], y.rev[-(1:r2)])
            }
          }
          
          coords <- rbind(coords, cbind(Var1  = x.rest[1], Var2 = unlist(y.rest[1])), cbind(Var1 = x.rest[2], Var2 = unlist(y.rest[2])))
          
        }
      }  
    }
  }
  return(coords)
}
