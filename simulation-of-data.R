################################################################################
## Simulation of the data, i.e. simulation of the Gaussian random fields      ##
################################################################################

##############################
## Load functions and packages
##############################

library(RandomFields) 
library(readxl) # to read excel files  
source("/Functions/functions-data-simulation.R") 
source("/Parameters.R")


########################################################
## Parameter combinations for the scenarios of the paper
########################################################

combs <- read_excel("Parametercombinations.xlsx")
  

###################################################################
## Gaussian random fields: for simulation of the correction factors
###################################################################

combs.corr <- subset(combs, combs$type == 1)
seeds <- rep(0, nrow(combs.corr)) # to save the seeds
for(ind in 1:nrow(combs.corr)){ 
  s <- sample(1:10000, 1) # draw a random seed
  set.seed(s)
  seeds[ind] <- s
  data <- sim.data(gridsize = grids[combs.corr$gridsize[ind],],
                   dist = dists,
                   dist.outlier = dists.outlier,
                   variogram = variograms[combs.corr$variogram[[ind]]][[1]],
                   param.variogram = params.variogram[combs.corr$param.variogram[ind],],
                   aniso.param = aniso.params[combs.corr$aniso[[ind]]],
                   param.outlier = param.outliers[combs.corr$param.outlier[[ind]]])
  save(data, file = paste0("Data/Correctionfactors/data_grid", 
                    combs.corr$gridsize[ind], "_vario",
                    combs.corr$variogram[ind], "_aniso", 
                    combs.corr$aniso[[ind]], "_", 
                    outlier.types[combs.corr$type[ind]], ".RData"))
}
save(seeds, file = "Data/Correctionfactors/seeds.RData")


##########################################################################
## Gaussian random fields: for simulation without outliers and consistency
##########################################################################

combs.kons <- subset(combs, combs$type == 2)
seeds <- rep(0, nrow(combs.kons))
for(ind in 1:nrow(combs.kons)){ 
  s <- sample(1:10000, 1)
  set.seed(s)
  seeds[ind] <- s
  data <- sim.data(gridsize = grids[combs.kons$gridsize[ind],],
                   dist = dists,
                   dist.outlier = dists.outlier,
                   variogram = variograms[combs.kons$variogram[[ind]]][[1]],
                   param.variogram = params.variogram[combs.kons$param.variogram[ind],],
                   aniso.param = aniso.params[combs.kons$aniso[[ind]]],
                   param.outlier = param.outliers[combs.kons$param.outlier[[ind]]])
  save(data, file = paste0("Data/Consistency/data_grid", 
                           combs.kons$gridsize[ind], "_vario",
                           combs.kons$variogram[ind], "_aniso", 
                           combs.kons$aniso[[ind]], "_", 
                           outlier.types[combs.kons$type[ind]], ".RData"))
}
save(seeds, file = "Data/Consistency/seeds.RData")

################################################################
## Gaussian random fields: for simulation with isolated outliers
################################################################

combs.iso <- subset(combs, combs$type == 3)
seeds <- rep(0, nrow(combs.iso))
for(ind in 1:nrow(combs.iso)){ 
  s <- sample(1:10000, 1)
  set.seed(s)
  seeds[ind] <- s
  data <- sim.data(gridsize = grids[combs.iso$gridsize[ind],],
                   dist = dists,
                   dist.outlier = dists.outlier,
                   variogram = variograms[combs.iso$variogram[[ind]]][[1]],
                   param.variogram = params.variogram[combs.iso$param.variogram[ind],],
                   aniso.param = aniso.params[combs.iso$aniso[[ind]]],
                   out.type = outlier.types[combs.iso$type[ind]],
                   amount = outlier.amount[combs.iso$amount[ind]],
                   param.outlier = unlist(param.outliers[combs.iso$param.outlier[[ind]]]))
  save(data, file = paste0("Data/Isolated/data_grid", 
                           combs.iso$gridsize[ind], "_vario",
                           combs.iso$variogram[ind], "_aniso", 
                           combs.iso$aniso[[ind]], "_", 
                           outlier.types[combs.iso$type[ind]], "_",
                           outlier.amount[combs.iso$amount[ind]], "_param",
                           combs.iso$param.outlier[[ind]], ".RData"))
}
save(seeds, file = "Data/Isolated/seeds.RData")


##############################################################
## Gaussian random fiels: for simulation with an outlier block
##############################################################

combs.block <- subset(combs, combs$type == 4)
seeds <- rep(0, nrow(combs.block))
for(ind in 1:nrow(combs.block)){ 
  s <- sample(1:10000, 1)
  set.seed(s)
  seeds[ind] <- s
  data <- sim.data(gridsize = grids[combs.block$gridsize[ind],],
                   dist = dists,
                   dist.outlier = dists.outlier,
                   variogram = variograms[combs.block$variogram[[ind]]][[1]],
                   param.variogram = params.variogram[combs.block$param.variogram[ind],],
                   aniso.param = aniso.params[combs.block$aniso[[ind]]],
                   out.type = outlier.types[combs.block$type[ind]],
                   amount = outlier.amount[combs.block$amount[ind]],
                   param.outlier = unlist(param.outliers[combs.block$param.outlier[[ind]]]),
                   blocktype = "quad")
  save(data, file = paste0("Data/Block/data_grid", 
                           combs.block$gridsize[ind], "_vario",
                           combs.block$variogram[ind], "_aniso", 
                           combs.block$aniso[[ind]], "_", 
                           outlier.types[combs.block$type[ind]], "_",
                           outlier.amount[combs.block$amount[ind]], "_param",
                           combs.block$param.outlier[[ind]], ".RData"))
}
save(seeds, file = "Data/Block/seeds.RData")


