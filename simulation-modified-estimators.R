################################################################################
## Simulation with the modified MCD estimators described in section 4.3       ##
################################################################################

##############################
## Load functions and packages
##############################

library(RandomFields) 
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
combs.mod <- subset(combs, combs$type != 3) # not for the case of isolated outliers  
combs.mod <- subset(combs.mod, combs.mod$variogram == 1) # only for the spherical variogram
# only for outlier distribution 1, 2 and 4 (scenarios of the paper)
combs.mod <- subset(combs.mod, combs.mod$param.outlier == 1 | combs.mod$param.outlier == 2| combs.mod$param.outlier == 4) 
# only for 0%, 5% and 15% outliers (scenarios of the paper)
combs.mod <- subset(combs.mod, combs.mod$amount == 1 | combs.mod$amount == 2 | combs.mod$amount == 4)

#############
## Simulation
#############

# needed to run the different scenarios on the HPC Cluster
ind <- as.integer(Sys.getenv("PBS_ARRAYID")) 
print(ind)

# simulation of correctionfactors of the modified estimators
if(combs.mod$type[ind] ==1){ 
# load the data for the corresponding scenario
load(file = paste0("Data/Correctionfactors/data_grid", 
                   combs.mod$gridsize[ind], "_vario",
                   combs.mod$variogram[ind], "_aniso", 
                   combs.mod$aniso[[ind]], "_", 
                   outlier.types[combs.mod$type[ind]], ".RData"))


plan(multicore) # need for the HPC cluster
res <- future_lapply(1:1000, function(i){
  sim.it(data$data[[i]], data$grid, hmaxs[[combs.mod$hmax[[ind]]]], cp = 1, mod = TRUE, mx = 8, my = 5)
}, future.seed = TRUE)

save(res, file = paste0("Modified/res_corr_mod_grid", 
                        combs.mod$gridsize[ind], "_vario",
                        combs.mod$variogram[ind], "_aniso", 
                        combs.mod$aniso[[ind]], "_", 
                        outlier.types[combs.mod$type[ind]], "_hmax",
                        combs.mod$hmax[[ind]], ".RData"))
rm(res)
} 

# simulations of the modified estimators for the grf without outliers
if(combs.mod$type[ind] ==2){ 
  # load the data for the corresponding scenario
  load(file = paste0("Data/Consistency/data_grid", 
                     combs.mod$gridsize[ind], "_vario",
                     combs.mod$variogram[ind], "_aniso", 
                     combs.mod$aniso[[ind]], "_", 
                     outlier.types[combs.mod$type[ind]], ".RData"))
  
  
  plan(multicore)# need for the HPC cluster
  res <- future_lapply(1:1000, function(i){
    sim.it(data$data[[i]], data$grid, hmaxs[[combs.mod$hmax[[ind]]]], cp = 1, mod = TRUE, mx = 8, my = 5)
  }, future.seed = TRUE)
  
  save(res, file = paste0("Modified/res_cons_mod_grid", 
                          combs.mod$gridsize[ind], "_vario",
                          combs.mod$variogram[ind], "_aniso", 
                          combs.mod$aniso[[ind]], "_", 
                          outlier.types[combs.mod$type[ind]], "_hmax",
                          combs.mod$hmax[[ind]], ".RData"))
  rm(res)
}

# simulations of the modified estimators for the grf with an outlier block
if(combs.mod$type[ind] ==4){ 
  # load the data for the corresponding scenario
  load(file = paste0("Data/Block/data_grid", 
                     combs.mod$gridsize[ind], "_vario",
                     combs.mod$variogram[ind], "_aniso", 
                     combs.mod$aniso[[ind]], "_", 
                     outlier.types[combs.mod$type[ind]], "_",
                     outlier.amount[combs.mod$amount[ind]], "_param",
                     combs.mod$param.outlier[[ind]], ".RData"))
  plan(multicore)
  res <- future_lapply(1:1000, function(i){
    sim.it(data$data[[i]], data$grid, hmaxs[[combs.mod$hmax[[ind]]]], cp = 1, mod = TRUE, mx = 8, my = 5)
  }, future.seed = TRUE)
  save(res, file = paste0("Modified/res_block_mod_grid", 
                          combs.mod$gridsize[ind], "_vario",
                          combs.mod$variogram[ind], "_aniso", 
                          combs.mod$aniso[[ind]], "_", 
                          outlier.types[combs.mod$type[ind]], "_",
                          outlier.amount[combs.mod$amount[ind]], "_param",
                          combs.mod$param.outlier[[ind]], "_hmax",
                          combs.mod$hmax[[ind]], ".RData"))
  rm(res)
  
  
}


################################
## Evaluation: Correctionfactors
################################

# to save the results
res.corr <- array( dim = c(10, 2, 4),
                   dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                "Matheron", "Genton", "MCD.diff.mod",
                                                "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re"),
                                   direction = c("S-N", "E-W"),
                                   grid = c("15,15", "25,25", "50,50", "75,75")))

for(ind in 1:4){ # for each scenarion, i.e. for each grid
  
  # load the output of the simulation
  load(file = paste0("Modified/res_corr_mod_grid",
                     combs.mod$gridsize[ind], "_vario",
                     combs.mod$variogram[ind], "_aniso",
                     combs.mod$aniso[[ind]], "_",
                     outlier.types[combs.mod$type[ind]], "_hmax",
                     combs.mod$hmax[[ind]], ".RData"))

  ## Calculate the true variogram values in the E-W and S-N direction
  # Geometric anisotropy: scaling and rotation (see section 4)
  Ts <- diag(c(1, 1/2))
  R <- matrix(c(cos(3*pi/8), -sin(3*pi/8), sin(3*pi/8), cos(3*pi/8)), nrow = 2)
  
  # lag vectors in S-N direction
  lag.vec.SN <- cbind(rep(0, hmaxs[1]), 1:hmaxs[1])
  lags.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
  
  # lag vectors in E-W direction
  lag.vec.EW <- cbind(1:hmaxs[2], rep(0, hmaxs[2]))
  lags.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
  
  # combs.corr$variogram = 1: spherical 
  sp.var <- function(h, a, c){ 
    if(h < a){res <- c*(((3*h)/(2*a)) - (1/2) * (h/a)^3)}
    if(h >= a){res <- c}
    if(h == 0){res <- 0}
    return(res)
  }
  true.sph <- list()
  true.sph[[1]] <- sapply(lags.SN, function(l) 2*sp.var(l, a = 5, c = 1))
  true.sph[[2]] <- sapply(lags.EW, function(l) 2*sp.var(l, a = 5, c = 1))

  
  ## Calculate correction factors (see section 4.1)
  # exclude the last lag, see explanations in section 4.1 
  # first: calculate estimation/true
  corr <- lapply(1:1000, function(l){ 
                 lapply(1:2, function(i){
                   apply(res[[l]][[i]][-nrow(res[[l]][[i]]),], 2, function(x){
                     x/true.sph[[i]][-nrow(res[[l]][[i]])]
                     })
                   })
          })
  
  # second: average over the different lags
  corr <- lapply(1:1000, function(l) lapply(1:2, function(i) colMeans(corr[[l]][[i]])))

  for(i in 1:2){ # for each direction (E-W and S-N)
   
    # third: average over all iterations
    val <- sapply(1:1000, function(l) corr[[l]][[i]])
    corr.val <- rowMeans(val)
    
    # calculate the final correction factor
    res.corr[,i, combs.mod$gridsize[ind]] <- 1/corr.val
  }
}

save(res.corr, file = "Modified/correctionfactors.RData")


###################################
## Evaluation: GRF without outliers
###################################

# load the simulated correctionfactors
load(file = "Modified/correctionfactors.RData")


# to save the estimation errors; needed for boxplot in 4.3 and calculation of 
# bias and sqrtMSE 
kons.est.error <- array(dim = c(10, 2, 7, 1000,  4),
                        dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                      "Matheron", "Genton", "MCD.diff.mod",
                                                      "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re"),
                                        direction = c("S-N", "E-W"),
                                        lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                        iteration = 1:1000,
                                        grid = c("15,15", "25,25", "50,50", "75,75")))

for(ind in 5:8){ # for each scenario, i.e for each grid
  load(file= paste0("Modified/res_cons_mod_grid",
                    combs.mod$gridsize[ind], "_vario",
                    combs.mod$variogram[ind], "_aniso",
                    combs.mod$aniso[[ind]], "_",
                    outlier.types[combs.mod$type[ind]], "_hmax",
                    combs.mod$hmax[[ind]], ".RData"))


  ## Calculate the true variogram values in the E-W and S-N direction
  # Geometric anisotropy: scaling and rotation (see section 4)
  Ts <- diag(c(1, 1/2)) # scaling
  R <- matrix(c(cos(3*pi/8), -sin(3*pi/8), sin(3*pi/8), cos(3*pi/8)), nrow = 2) # rotation
  
  # lag vectors in S-N direction
  lag.vec.SN <- cbind(rep(0, hmaxs[1]), 1:hmaxs[1])
  lags.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
  
  # lag vectors in E-W direction
  lag.vec.EW <- cbind(1:hmaxs[2], rep(0, hmaxs[2]))
  lags.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))

  # combs.corr$variogram = 1: spherical 
  sp.var <- function(h, a, c){ 
    if(h < a){res <- c*(((3*h)/(2*a)) - (1/2) * (h/a)^3)}
    if(h >= a){res <- c}
    if(h == 0){res <- 0}
    return(res)
  }
  true.sph <- list()
  true.sph[[1]] <- sapply(lags.SN, function(l) 2*sp.var(l, a = 5, c = 1))
  true.sph[[2]] <- sapply(lags.EW, function(l) 2*sp.var(l, a = 5, c = 1))

  # multiplicate with the correctionfactors
  fac <- res.corr[,,combs.mod$gridsize[ind]]
  res.c <- lapply(res, function(i){
    z.res <- list()
    for(r in 1:2){
      z.res[[r]] <- sapply(1:10, function(s) i[[r]][,s] * fac[s,r])
    }
    return(z.res)
  })

  # calculate the estimation errors (estimation - true)
  est.error <- lapply(1:1000, function(l){
    lapply(1:2, function(r){
      apply(res.c[[l]][[r]], 2, function(x) x - true.sph[[r]])
    })})

  for(l in 1:1000){
    for(r in 1:2){

      t(est.error[[l]][[r]]) -> kons.est.error[, r, , l, combs.mod$gridsize[ind]]
    }
  }
}
save(kons.est.error, file = "Modified/estimation_errors.RData")


## Calculate bias, estimated mean and sqrtMSE
kons.bias <- array( dim = c(10, 2, 7, 4),
                    dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                 "Matheron", "Genton", "MCD.diff.mod",
                                                 "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re"),
                                    direction = c("S-N", "E-W"),
                                    lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                    grid = c("15,15", "25,25", "50,50", "75,75")))

kons.mean <- array( dim = c(10, 2, 7, 4),
                    dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                 "Matheron", "Genton", "MCD.diff.mod",
                                                 "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re"),
                                    direction = c("S-N", "E-W"),
                                    lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                    grid = c("15,15", "25,25", "50,50", "75,75")))

kons.sqrtMSE <- array( dim = c(10, 2, 7, 4),
                       list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                         "Matheron", "Genton", "MCD.diff.mod",
                                         "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re"),
                            direction = c("S-N", "E-W"),
                            lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                            grid = c("15,15", "25,25", "50,50", "75,75")))


for(ind in 5:8){  # for each scenarion, i.e for each grid
  load(file= paste0("Modified/res_cons_mod_grid",
                    combs.mod$gridsize[ind], "_vario",
                    combs.mod$variogram[ind], "_aniso",
                    combs.mod$aniso[[ind]], "_",
                    outlier.types[combs.mod$type[ind]], "_hmax",
                    combs.mod$hmax[[ind]], ".RData"))


  ## Calculate the true variogram values in the E-W and S-N direction
  # Geometric anisotropy: scaling und rotation (see section 4)
  Ts <- diag(c(1, 1/2)) # scaling
  R <- matrix(c(cos(3*pi/8), -sin(3*pi/8), sin(3*pi/8), cos(3*pi/8)), nrow = 2) #rotation
  
  # lag vectors in S-N direction
  lag.vec.SN <- cbind(rep(0, hmaxs[1]), 1:hmaxs[1])
  lags.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
  
  # lag vectors in E-W direction
  lag.vec.EW <- cbind(1:hmaxs[2], rep(0, hmaxs[2]))
  lags.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
  # combs.corr$variogram = 1: spherical 
  sp.var <- function(h, a, c){ 
    if(h < a){res <- c*(((3*h)/(2*a)) - (1/2) * (h/a)^3)}
    if(h >= a){res <- c}
    if(h == 0){res <- 0}
    return(res)
  }
  true.sph <- list()
  true.sph[[1]] <- sapply(lags.SN, function(l) 2*sp.var(l, a = 5, c = 1))
  true.sph[[2]] <- sapply(lags.EW, function(l) 2*sp.var(l, a = 5, c = 1))
  
  # multiplicate with the correction factors
  fac <- res.corr[,, combs.mod$gridsize[ind]]
  res.c <- lapply(res, function(i){
    z.res <- list()
    for(r in 1:2){
      z.res[[r]] <- sapply(1:10, function(s) i[[r]][,s] * fac[s,r])
    }
    return(z.res)
  })

  # calculate the squared estimation errors ((estimation - true)^2)
  sq.error <- lapply(1:1000, function(l){ # needed for sqrtMSE
    lapply(1:2, function(r){
      apply(res.c[[l]][[r]], 2, function(x) (x - true.sph[[r]])^2)
    })})

  for(r in 1:2){ # for S-N and E-W
    total <- 0
    MSE <- 0
    for(l in 1:1000){

      total <- total + res.c[[l]][[r]] # calculate the sum of the estimates 
      MSE <- MSE + sq.error[[l]][[r]]  # calculate the MSE
    }
    total <- total/1000
    MSE <- sqrt(MSE/1000)
    Bias <- total - true.sph[[r]] # calculate the bias

    t(total) -> kons.mean[, r, , combs.mod$gridsize[ind]]
    t(Bias) -> kons.bias[, r, , combs.mod$gridsize[ind]]
    t(MSE) -> kons.sqrtMSE[, r, , combs.mod$gridsize[ind]]
  }
}
save(kons.mean, file = "Modified/estimation_mean.RData")
save(kons.bias, file = "Modified/estimation_bias.RData")
save(kons.sqrtMSE, file = "Modified/estimation_sqrtMSE.RData")

#####################################
## Evaluation: GRF with blockoutliers
#####################################

# load the simulated correctionfactors
load(file = "Modified/correctionfactors.RData")


## Calculate bias, estimated mean and sqrtMSE
kons.bias <- array(dim = c(10, 2, 7,  2, 2),
                   dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                "Matheron", "Genton", "MCD.diff.mod",
                                                "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re"),
                                   direction = c("S-N", "E-W"),
                                   lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                   amount = c(0.05, 0.15),
                                   distribution = 1:2))

kons.mean  <- array(dim = c(10, 2, 7,  2, 2),
                    dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                 "Matheron", "Genton", "MCD.diff.mod",
                                                 "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re"),
                                    direction = c("S-N", "E-W"),
                                    lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                    amount = c(0.05,  0.15),
                                    distribution = 1:2))

kons.sqrtMSE <- array(dim = c(10, 2, 7,  2, 2),
                      dimnames = list(estimator = c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re",
                                                   "Matheron", "Genton", "MCD.diff.mod",
                                                   "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re"),
                                      direction = c("S-N", "E-W"),
                                      lag = c("l1", "l2", "l3", "l4", "l5", "l6", "l7"),
                                      amount = c(0.05, 0.152),
                                      distribution = 1:2))


for(ind in 9:12){ # for each scenario
  load(file=  paste0("Modified/res_block_mod_grid",
                     combs.mod$gridsize[ind], "_vario",
                     combs.mod$variogram[ind], "_aniso",
                     combs.mod$aniso[[ind]], "_",
                     outlier.types[combs.mod$type[ind]], "_",
                     outlier.amount[combs.mod$amount[ind]], "_param",
                     combs.mod$param.outlier[[ind]], "_hmax",
                     combs.mod$hmax[[ind]], ".RData"))


  ## Calculate the true variogram values in the E-W and S-N direction
  # Geometric anisotropy: scaling and rotation (see section 4)
  Ts <- diag(c(1, 1/2))
  R <- matrix(c(cos(3*pi/8), -sin(3*pi/8), sin(3*pi/8), cos(3*pi/8)), nrow = 2)
  
  # lag vectors in S-N direction
  lag.vec.SN <- cbind(rep(0, hmaxs[1]), 1:hmaxs[1])
  lags.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
  
  # lag vectors in E-W direction
  lag.vec.EW <- cbind(1:hmaxs[2], rep(0, hmaxs[2]))
  lags.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))

  # combs.corr$variogram = 1: spherical 
  sp.var <- function(h, a, c){ 
    if(h < a){res <- c*(((3*h)/(2*a)) - (1/2) * (h/a)^3)}
    if(h >= a){res <- c}
    if(h == 0){res <- 0}
    return(res)
  }
  true.sph <- list()
  true.sph[[1]] <- sapply(lags.SN, function(l) 2*sp.var(l, a = 5, c = 1))
  true.sph[[2]] <- sapply(lags.EW, function(l) 2*sp.var(l, a = 5, c = 1))
  

  # multiplicate with the correction factors
  fac <- res.corr[,, combs.mod$gridsize[ind]]
  res.c <- lapply(res, function(i){
    z.res <- list()
    for(r in 1:2){
      z.res[[r]] <- sapply(1:10, function(s) i[[r]][,s] * fac[s,r])
    }
    return(z.res)
  })


  # calculate the squared estimation errors ((estimation - true)^2)
  sq.error <- lapply(1:1000, function(l){ # needed for sqrtMSE
    lapply(1:2, function(r){
      apply(res.c[[l]][[r]], 2, function(x) (x - true.sph[[r]])^2)
    })})

  for(r in 1:2){
    summe <- 0
    MSE <- 0
    for(l in 1:1000){
      summe <- summe + res.c[[l]][[r]] # calculate the sum of the estimates 
      MSE <- MSE + sq.error[[l]][[r]]  # calculate the MSE
    }

    summe <- summe/1000
    MSE <- sqrt(MSE/1000)
    Bias <- summe - true.sph[[r]] # calculate the bias

    param <- ifelse(combs.mod$param.outlier[ind] ==2, 1, 2)
    amo <- ifelse(combs.mod$amount[ind] ==2, 1, 2)
    t(summe) -> kons.mean[, r, ,  amo, param]
    t(Bias) -> kons.bias[, r, ,  amo, param]
    t(MSE) -> kons.sqrtMSE[, r, ,  amo, param]
  }
}
save(kons.mean, file = "Modified/estimation_mean_block.RData")
save(kons.bias, file = "Modified/estimation_bias_block.RData")
save(kons.sqrtMSE, file = "Modified/estimation_sqrtMSE_block.RData")
