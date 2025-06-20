################################################################################
## Simulation of the GRF with an outlier block (see section 4.4)              ##
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
combs.block <- subset(combs, combs$type == 4) 

#############
## Simulation
#############

# needed to run the different scenarios on the HPC Cluster
ind <- as.integer(Sys.getenv("PBS_ARRAYID")) 
print(ind)

# load the data for the corresponding scenario
load(file = paste0("Data/Block/data_grid", 
                   combs.block$gridsize[ind], "_vario",
                   combs.block$variogram[ind], "_aniso", 
                   combs.block$aniso[[ind]], "_", 
                   outlier.types[combs.block$type[ind]], "_",
                   outlier.amount[combs.block$amount[ind]], "_param",
                   combs.block$param.outlier[[ind]], ".RData"))

plan(multicore) # need for the HPC cluster
res <- future_lapply(1:1000, function(i){
  sim.it(data$data[,i], expand.grid(1:grids[combs.block$gridsize[ind],"nx"], 1:grids[combs.block$gridsize[ind],"ny"]), hmaxs, cp = 1)
}, future.seed = TRUE)
save(res, file = paste0("Block/res_block_grid", 
                        combs.block$gridsize[ind], "_vario",
                        combs.block$variogram[ind], "_aniso", 
                        combs.block$aniso[[ind]], "_", 
                        outlier.types[combs.block$type[ind]], "_",
                        outlier.amount[combs.block$amount[ind]], "_param",
                        combs.block$param.outlier[[ind]], "_hmax",
                        combs.block$hmax[[ind]], ".RData"))
rm(res)


#############
## Evaluation
#############

# load the simulated correctionfactors
load(file = "Correctionfactors/correctionfactors.RData")

## Calculate bias, estimated mean and sqrtMSE
kons.bias <- array(dim = c(6, 4, 7,  4, 4, 2),
                   dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                 "Matheron", "Genton"),
                                   direction = c("S-N", "E-W", "SW-NW", "SE-NW"),
                                   lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                   amount = c(0.05, 0.10, 0.15, 0.25),
                                   distribution = 1:4,
                                   grid = c("15,15", "60,60")))


kons.mean  <- array(dim = c(6, 4, 7,  4, 4, 2),
                    dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                  "Matheron", "Genton"),
                                    direction = c("S-N", "E-W", "SW-NW", "SE-NW"),
                                    lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                    amount = c(0.05, 0.10, 0.15, 0.25),
                                    distribution = 1:4,
                                    grid = c("15,15", "60,60")))

kons.sqrtMSE <- array(dim = c(6, 4, 7,  4, 4, 2),
                      dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                    "Matheron", "Genton"),
                                      direction = c("S-N", "E-W", "SW-NW", "SE-NW"),
                                      lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                      amount = c(0.05, 0.10, 0.15, 0.25),
                                      distribution = 1:4,
                                      grid = c("15,15", "60,60")))

for(ind in 1:nrow(combs.block)){ # for each scenario
  load(file=  paste0("Block/res_block_grid",
                     combs.block$gridsize[ind], "_vario",
                     combs.block$variogram[ind], "_aniso",
                     combs.block$aniso[[ind]], "_",
                     outlier.types[combs.block$type[ind]], "_",
                     outlier.amount[combs.block$amount[ind]], "_param",
                     combs.block$param.outlier[[ind]], "_hmax",
                     combs.block$hmax[[ind]], ".RData"))

  ## Calculate the true variogram values in the four different directions
  # Geometric anisotropy: scaling and rotation (see section 4)
  Ts <- diag(c(1, 1/2))
  R <- matrix(c(cos(3*pi/8), -sin(3*pi/8), sin(3*pi/8), cos(3*pi/8)), nrow = 2)

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


  fac <- res.corr[,, combs.block$gridsize[ind], combs.block$variogram[ind]]
  res.c <- lapply(res, function(i){
    z.res <- list()
    for(r in 1:4){
      z.res[[r]] <- sapply(1:6, function(s) i[[r]][,s] * fac[s,r])
    }
    return(z.res)
  })

  # calculate the squared estimation errors ((estimation - true)^2)
  sq.error <- lapply(1:1000, function(l){ # needed for sqrtMSE
    lapply(1:4, function(r){
      apply(res.c[[l]][[r]], 2, function(x) (x - true.sph[[r]])^2)
    })})

  for(r in 1:4){ # for each direction
    total <- 0
    MSE <- 0
    for(l in 1:1000){
      total <- total + res.c[[l]][[r]] # calculate the sum of the estimates
      MSE <- MSE + sq.error[[l]][[r]]  # calculate the MSE
    }

    total <- total/1000
    MSE <- sqrt(MSE/1000)
    Bias <- total - true.sph[[r]] # calculate the bias

    # in the SW-NE and SE-NW direction the number of lags is smaller than in
    # the other two direction -> therefore fill up with NAs
    if(nrow(total) <= 7){
      d <- 7 - nrow(total)
      total <- rbind(total, matrix(rep(NA, 6*d), ncol = 6))
      MSE <- rbind(MSE, matrix(rep(NA, 6*d), ncol = 6))
      Bias <- rbind(Bias, matrix(rep(NA, 6*d), ncol = 6))
    }
    
    g <- ifelse(combs.block$gridsize == 1, 1, 2)

    t(total) -> kons.mean[, r, ,  combs.block$amount[ind]-1, combs.block$param.outlier[ind]-1, g]
    t(Bias) -> kons.bias[, r, ,  combs.block$amount[ind]-1, combs.block$param.outlier[ind]-1, g]
    t(MSE) -> kons.sqrtMSE[, r, ,  combs.block$amount[ind]-1, combs.block$param.outlier[ind]-1, g]
  }
}

save(kons.mean, file = "Block/estimation_mean.RData")
save(kons.bias, file = "Block/estimation_bias.RData")
save(kons.sqrtMSE, file = "Block/estimation_sqrtMSE.RData")

