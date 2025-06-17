################################################################################
## Simulation of the correction factors (see section 4.1)                     ##
################################################################################

##############################
## Load functions and packages
##############################

library(geoR) 
library(future) # for parallelisation 
library(future.apply) # for parallelisation  
library(readxl) # to read excel files 
source("Functions/MCD-variogram-estimator.R") 
source("Functions/function-for-simulation.R")
source("Functions/Auxiliary-functions.R")
source("Parameters.R")

######################################################################
## Parameter combinations for the simulation of the correction factors
######################################################################

combs <- read_excel("Parametercombinations.xlsx")
combs.corr <- subset(combs, combs$type == 1) 

#############
## Simulation
#############

# needed to run the different scenarios on the HPC Cluster
ind <- as.integer(Sys.getenv("PBS_ARRAYID")) 
print(ind)

# load the data for the corresponding scenario
load(file = paste0("Data/Correctionfactors/data_grid",
                   combs.corr$gridsize[ind], "_vario",
                   combs.corr$variogram[ind], "_aniso",
                   combs.corr$aniso[[ind]], "_",
                   outlier.types[combs.corr$type[ind]], ".RData"))


plan(multicore) # need for the HPC cluster
res <- future_lapply(1:1000, function(i){
  sim.it(data$data[,i], expand.grid(1:grids[combs.corr$gridsize[ind],"nx"], 1:grids[combs.corr$gridsize[ind],"ny"]), hmaxs, cp = 1)
}, future.seed = TRUE)
save(res, file = paste0("Correctionfactors/res_corr_grid", 
                        combs.corr$gridsize[ind], "_vario",
                        combs.corr$variogram[ind], "_aniso", 
                        combs.corr$aniso[[ind]], "_", 
                        outlier.types[combs.corr$type[ind]], "_hmax",
                        combs.corr$hmax[[ind]], ".RData"))
rm(res)


#############
## Evaluation
#############

# to save the results
res.corr <- array( dim = c(6, 4, 5, 3),
                   dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                "Matheron", "Genton"),
                                   direction = c("S-N", "E-W", "SW-NE",  "SE-NW"),
                                   grid = c("15,15", "25,25", "50,50", "60, 60", "75,75"),
                                   variogram = c("sph", "exp", "gauss")))

for(ind in 1:nrow(combs.corr)){ # for each szenario
  
  # load the output of the simulation
  load(file = paste0("Correctionfactors/res_corr_grid",
                     combs.corr$gridsize[ind], "_vario",
                     combs.corr$variogram[ind], "_aniso",
                     combs.corr$aniso[[ind]], "_",
                     outlier.types[combs.corr$type[ind]], "_hmax",
                     combs.corr$hmax[[ind]], ".RData"))

  
  ## Calculate the true variogram values in the four different directions
  
  # Geometric anisotropy: scaling and rotation (see section 4)
  Ts <- diag(c(1, 1/2)) # scaling
  R <- matrix(c(cos(3*pi/8), -sin(3*pi/8), sin(3*pi/8), cos(3*pi/8)), nrow = 2) # rotation

  # lag vectors in S-N direction
  lag.vec.SN <- cbind(rep(0, hmaxs[1]), 1:hmaxs[1])
  lags.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))

  # lag vectors in E-W direction
  lag.vec.EW <- cbind(1:hmaxs[2], rep(0, hmaxs[2]))
  lags.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))

  # lag vectors in SW-NE direction
  lag.vec.SWNE <- cbind(1:hmaxs[3], 1:hmaxs[3])
  lags.SWNE <- apply(lag.vec.SWNE, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))

  # lag vectors in SE-NW direction
  lag.vec.SENW <- cbind(1:hmaxs[4], -(1:hmaxs[4]))
  lags.SENW <- apply(lag.vec.SENW, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))


  # combs.corr$variogram = 1: spherical 
  sp.var <- function(h, a, c){ 
    if(h < a){res <- c*(((3*h)/(2*a)) - (1/2) * (h/a)^3)}
    if(h >= a){res <- c}
    if(h == 0){res <- 0}
    return(res)
  }
  true.sph <- list()
  true.sph[[1]] <- sapply(lags.SN, function(l) 2*(nugget + sp.var(l, a = 5, c = 0.4)))
  true.sph[[2]] <- sapply(lags.EW, function(l) 2*(nugget + sp.var(l, a = 5, c = 0.4)))
  true.sph[[3]] <- sapply(lags.SWNE, function(l) 2*(nugget + sp.var(l, a = 5, c = 0.4)))
  true.sph[[4]] <- sapply(lags.SENW, function(l) 2*(nugget + sp.var(l, a = 5, c = 0.4)))
  
  # combs.corr$variogram = 2: exponential
  exp.var <- function(h, a, c){
    res <- c*(1-exp(-abs(h)/a))
    return(res)
  }
  true.exp <- list()
  true.exp[[1]] <- sapply(lags.SN, function(l) 2*(nugget + exp.var(l, a = 3, c = 0.4)))
  true.exp[[2]] <- sapply(lags.EW, function(l) 2*(nugget + exp.var(l, a = 3, c = 0.4)))
  true.exp[[3]] <- sapply(lags.SWNE, function(l) 2*(nugget + exp.var(l, a = 3, c = 0.4)))
  true.exp[[4]] <- sapply(lags.SENW, function(l) 2*(nugget + exp.var(l, a = 3, c = 0.4)))
  
  # combs.corr$variogram = 3: gauss
  gauss.var <- function(h, a, c){
    res <- c*(1-exp(-abs(h)^2/a^2))
    return(res)
  }
  true.gauss <- list()
  true.gauss[[1]] <- sapply(lags.SN, function(l) 2*(nugget + gauss.var(l, a = 3, c = 0.4)))
  true.gauss[[2]] <- sapply(lags.EW, function(l) 2*(nugget + gauss.var(l, a = 3, c = 0.4)))
  true.gauss[[3]] <- sapply(lags.SWNE, function(l) 2*(nugget + gauss.var(l, a = 3, c = 0.4)))
  true.gauss[[4]] <- sapply(lags.SENW, function(l) 2*(nugget + gauss.var(l, a = 3, c = 0.4)))

  # save the different true variograms
  true.var <- list(true.sph, true.exp, true.gauss)
  

  ## Calculate correction factors (see section 4.1)
  # exclude the last lag, see explanations in section 4.1 
  # first: calculate estimation/true
  corr <- lapply(1:1000, function(l){ 
                 lapply(1:4, function(i){ 
                   apply(res[[l]][[i]][-nrow(res[[l]][[i]]),], 2, function(x){ 
                     x/true.var[[combs.corr$variogram[ind]]][[i]][-nrow(res[[l]][[i]])]
                     })
                   })
           })
  
  # second: average over the different lags
  corr <- lapply(1:1000, function(l) lapply(1:4, function(i) colMeans(corr[[l]][[i]])))
  
  for(i in 1:4){ # for each direction
    
    # third: average over all iterations
    val <- sapply(1:1000, function(l) corr[[l]][[i]])
    corr.val <- rowMeans(val)
    
    # calculate the final correction factor
    res.corr[,i, combs.corr$gridsize[ind], combs.corr$variogram[ind]] <- 1/corr.val
  }
}

# save the results
save(res.corr, file = "Correctionfactors/correctionfactors.RData")
