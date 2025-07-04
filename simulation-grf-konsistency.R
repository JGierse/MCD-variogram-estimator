################################################################################
## Simulation of the GRF without outliers (see section 4.2 & 4.3)             ##
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
combs.cons <- subset(combs, combs$type == 2) 

#############
## Simulation
#############

# needed to run the different scenarios on the HPC Cluster
ind <- as.integer(Sys.getenv("PBS_ARRAYID")) 
print(ind)

# load the data for the corresponding scenario
load(file = paste0("Data/Consistency/data_grid",
                   combs.cons$gridsize[ind], "_vario",
                   combs.cons$variogram[ind], "_aniso",
                   combs.cons$aniso[[ind]], "_",
                   outlier.types[combs.cons$type[ind]], ".RData"))


plan(multicore)  # need for the HPC cluster
res <- future_lapply(1:1000, function(i){
  sim.it(data$data[,i], expand.grid(1:grids[combs.cons$gridsize[ind],"nx"], 1:grids[combs.cons$gridsize[ind],"ny"]), hmaxs, cp = 1)
}, future.seed = TRUE)
save(res, file = paste0("Consistency/res_cons_grid", 
                        combs.cons$gridsize[ind], "_vario",
                        combs.cons$variogram[ind], "_aniso", 
                        combs.cons$aniso[[ind]], "_", 
                        outlier.types[combs.cons$type[ind]], "_hmax",
                        combs.cons$hmax[[ind]], ".RData"))


#############
## Evaluation
#############

# load the simulated correctionfactors
load(file = "Correctionfactors/correctionfactors.RData")


# to save the estimation errors; needed for boxplot in 4.3 and calculation of
# bias and sqrtMSE
kons.est.error <- array(dim = c(6, 4, 7, 1000,  5),
                        dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                     "Matheron", "Genton"),
                                        direction = c("S-N", "E-W", "SW-NE", "SE-NW"),
                                        lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                        iteration = 1:1000,
                                        grid = c("15,15", "25,25", "50,50", "60,60", "75,75")))

for(ind in 1:nrow(combs.cons)){ # for each scenario
  load(file= paste0("Consistency/res_cons_grid",
                    combs.cons$gridsize[ind], "_vario",
                    combs.cons$variogram[ind], "_aniso",
                    combs.cons$aniso[[ind]], "_",
                    outlier.types[combs.cons$type[ind]], "_hmax",
                    combs.cons$hmax[[ind]], ".RData"))


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
  
  # multiplicate with the correction factors
  fac <- res.corr[,,combs.cons$gridsize[ind], combs.cons$variogram[ind]]
  res.c <- lapply(res, function(i){
    z.res <- list()
    for(r in 1:4){
      z.res[[r]] <- sapply(1:6, function(s) i[[r]][,s] * fac[s,r])
    }
    return(z.res)
  })


  # calculate the estimation errors (estimation - true)
  est.error <- lapply(1:1000, function(l){
    lapply(1:4, function(r){
      apply(res.c[[l]][[r]], 2, function(x) x - true.sph[[r]])
    })})

  # in the SW-NE and SE-NW direction the number of lags is smaller than in
  # the other two direction -> therefore fill up with NAs
  for(l in 1:1000){
    for(r in 1:4){
      if(nrow(est.error[[l]][[r]]) <= 7){
        d <- 7 - nrow(est.error[[l]][[r]])
        est.error[[l]][[r]] <- rbind(est.error[[l]][[r]], matrix(rep(NA, 6*d), ncol = 6))
      }


      t(est.error[[l]][[r]]) -> kons.est.error[, r, , l, combs.cons$gridsize[ind]]
    }
  }
}
save(kons.est.error, file = "Consistency/estimation_errors.RData")


## Calculate bias, estimated mean and sqrtMSE
kons.bias <- array( dim = c(6, 4, 7, 5),
                    dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                 "Matheron", "Genton"),
                                    direction = c("S-N", "E-W", "SW-NE", "SE-NW"),
                                    lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                    grid = c("15,15", "25,25", "50,50", "60,60", "75,75")))

kons.mean <- array( dim = c(6, 4, 7, 5),
                    dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                 "Matheron", "Genton"),
                                    direction = c("S-N", "E-W", "SW-NE", "SE-NW"),
                                    lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                    grid = c("15,15", "25,25", "50,50", "60,60", "75,75")))

kons.sqrtMSE <- array( dim = c(6, 4, 7, 5),
                       dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                    "Matheron", "Genton"),
                                       direction = c("S-N", "E-W", "SW-NE", "SE-NW"),
                                       lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                       grid = c("15,15", "25,25", "50,50", "60,60", "75,75")))


for(ind in 1:nrow(combs.cons)){ # for each scenario
  load(file= paste0("Consistency/res_cons_grid",
                    combs.cons$gridsize[ind], "_vario",
                    combs.cons$variogram[ind], "_aniso",
                    combs.cons$aniso[[ind]], "_",
                    outlier.types[combs.cons$type[ind]], "_hmax",
                    combs.cons$hmax[[ind]], ".RData"))

  ## Calculate the true variogram values in the four different directions
  # Geometric anisotropy: scaling and rotation (see section 4)
  Ts <- diag(c(1, 1/2)) # scaling
  R <- matrix(c(cos(3*pi/8), -sin(3*pi/8), sin(3*pi/8), cos(3*pi/8)), nrow = 2) #rotation

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

  # multiplicate with the correction factors
  fac <- res.corr[,, combs.cons$gridsize[ind], combs.cons$variogram[ind]]
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
    if(nrow(est.error[[l]][[r]]) <= 7){
      d <- 7 - nrow(total)

      total <- rbind(total, matrix(rep(NA, 6*d), ncol = 6))
      MSE <- rbind(MSE, matrix(rep(NA, 6*d), ncol = 6))
      Bias <- rbind(Bias, matrix(rep(NA, 6*d), ncol = 6))
    }

    t(total) -> kons.mean[, r, , combs.cons$gridsize[ind]]
    t(Bias) -> kons.bias[, r, , combs.cons$gridsize[ind]]
    t(MSE) -> kons.sqrtMSE[, r, , combs.cons$gridsize[ind]]
  }
}
save(kons.mean, file = "Consistency/estimation_mean.RData")
save(kons.bias, file = "Consistency/estimation_bias.RData")
save(kons.sqrtMSE, file = "Consistency/estimation_sqrtMSE.RData")
