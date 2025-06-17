################################################################################
## Parameter values for the simulations presented in the paper                ##
################################################################################

## grid sizes
grids <- rbind(c("nx" = 15, "ny" = 15),
               c("nx" = 25, "ny" = 25),
               c("nx" = 50, "ny" = 50),
               c("nx" = 60, "ny" = 60),
               c("nx" = 75, "ny" = 75))

## distribution of the outliers
dists.outlier <- rnorm

## variogram models
variograms <- list("spherical", "exponential", "gaussian")

## parameters of the variogram models
params.variogram <- rbind(c("var" = 0.4, "scale" = 5),
                          c("var" = 0.4, "scale" = 3))

## nugget effects
nugget <- c(0.1)

## Parameters for the anisotropy
aniso.params <- c("angle" = 3*pi/8, "ratio" = 2)

## type of outliers
outlier.types <- c("non", "non", "isolated", "block")
# first non is for correction factors
# second non for the simulation without outliers

## amount of outliers
outlier.amount <- c(0, 0.05, 0.10, 0.15, 0.25)

## parameters for the outliers
param.outliers <- list(NULL, c(1.5, 1), c(3,1), c(5, 1), c(0,4))

## number of lags
hmaxs <- c(7,7,5,5)


