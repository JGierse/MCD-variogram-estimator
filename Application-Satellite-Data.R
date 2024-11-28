################################################################################
## Simulation with the modified MCD estimators described in section 4.3       ##
################################################################################

##############################
## Load functions and packages
##############################

# Necessary packages for preparing the satellite data
library(gdalcubes)
library(magrittr)
library(xts)
library(magick)

# Necessary for the variogram estimation
library(tidyverse)
library(robustbase)

# Function for calculation of the MCD variogram estimators
source("Functions/MCD-variogram-estimator.R")
source("Functions/function-for-simulation.R")
source("Functions/Auxiliary-functions.R")

# graphs
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

####################################
## Preparation of the satellite data
####################################

## Load the satellite data
## Data can be downloaded from https://earthexplorer.usgs.gov/
## the used data here have the following dimensions in space and time:
# A data cube view object
# 
# Dimensions:
#   low              high count pixel_size
# t        2013-04-12        2019-11-30  2424        P1D
# y -1096876.57818711 -993876.578187112   206        500
# x -7370181.54841633 -7235181.54841633   270        500
# 
# SRS: "EPSG:3857"
# Temporal aggregation method: "first"
# Spatial resampling method: "near"
#IMAGE_DIR = "L8_cropped"  # Path of the data 

# load the structure of the data
col <- create_image_collection(list.files(IMAGE_DIR, recursive = TRUE, pattern=".tif", full.names  = TRUE), "L8_SR")

## Selection of pixels based on quality band
# only pixel classified as clear
L8.clear_mask = image_mask("PIXEL_QA", values=c(322, 386, 834, 898, 1346, 324, 388, 836, 900, 1348), invert = TRUE)

# Define what the data set should look like
v = cube_view(srs="EPSG:3857", extent=col, dx=500, dy=500, dt="P1D")

# Transform data into the required format and select the required spectral bands
comp_mask <- raster_cube(col, v, mask = L8.clear_mask) %>% select_bands((names(raster_cube(col,v))[2:8]))
comp <- raster_cube(col, v) %>% select_bands((names(raster_cube(col,v))[2:8]))
v.date <-cube_view(view = v, extent=list(left = -7331535-400*30, right = -7331535+400*30, 
                                         bottom = -999642-400*30, top = -999642 + 150*30, t0 = "2016-10-20", t1 = "2016-10-20"), dx=30, dy=30, dt="P1D") 
comp.date <- raster_cube(col, v.date) %>% select_bands(c("B02", "B03", "B04"))


# Function for calculating the NDVI for a subregion
NDVI_sub <- function(pixel_x, pixel_y, time, mask = NULL){
  delta = 30
  v.sub = cube_view(view = v, extent=list(left=pixel_x-(20 * delta), right=pixel_x+(40 * delta),
                                          bottom=pixel_y-(20 * delta),top=pixel_y+(40 * delta),
                                          t0 = time, t1 = time), dx=delta, dy=delta, dt="P1D")
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI")  -> tseries_plot

  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI") %>%
    as_array() %>% drop -> tseries

  return(list(tseries = tseries, tseries_plot = tseries_plot))
}

#############
## 60x60 grid 
#############

# all data in the subregion
delta = 30
v.sub1 = cube_view(view = v, extent=list(left=-7331535-20*delta, right=-7331535+40*delta, 
                                        bottom=-999642-20*delta,top=-999642+40*delta,
                                        t0 = "2016-10-20", t1 = "2016-10-20"), dx=delta, dy=delta, dt="P1D") 
sub1 <- raster_cube(col, v.sub1) %>% select_bands(c("B02", "B03", "B04"))

# only clear data in the subregion
sub.clear1 <- raster_cube(col, v.sub1, mask = L8.clear_mask) %>% select_bands(c("B02", "B03", "B04"))

########################################################
## Calculation: complete region (without missing values) 
########################################################

# calculate the pixlewise NDVI for the subregion
reg1 <-  NDVI_sub(-7331535, -999642, "2016-10-20")

# extract the NDVI Index
reg1.data = reg1$tseries

# rename the y-cordinates
rownames(reg1.data) <- sapply(60:1, function(x) paste0("y", x))

# add the coordinates
reg1.data <- cbind(sapply(60:1, function(x) paste0("y", x)), reg1.data)

# rename the x-coordinates
colnames(reg1.data) <- c("y", sapply(1:60, function(x) paste0("x", x)))

# data set with a row for each data point
reg1.data.long <- as_tibble(reg1.data) %>% pivot_longer(cols = starts_with("x"), values_to = "data") %>%  
  mutate(data = as.numeric(data), ID = row_number()) %>% rename(x = name) %>% 
  mutate(x = as.numeric(sub("x", "", x)), y = as.numeric(sub("y", "", y))) 

# estimate the standard deviation of the data with an robust estimator
sd.all <- mad(as.vector(reg1$tseries))

# standardaize the data with the standard deviation fo the complete data 
reg1.data.long <- reg1.data.long %>% mutate(data.sd = data/sd.all)

# Graphic of the NDVI of the subregion (Figure 9)
pdf("Graphs/NDVI.pdf")
plot(reg1$tseries_plot, key.pos = 4, col = rev(hcl.colors(10, palette = "Greens 2")))
dev.off()

# Variogram estimation 
varo.reg1 <- sim.it(reg1.data.long$data.sd, cbind(reg1.data.long$y, reg1.data.long$x), hmax = c(4,4,3,3), cp = 1)


################################################
## Calculation: only clear data (missing values) 
################################################

# calculate the pixlewise NDVI for the subregion with only clear data
reg1.clear <- NDVI_sub(-7331535, -999642, "2016-10-20", L8.clear_mask)

# extract the NDVI Index
reg1.clear.data = reg1.clear$tseries 

# rename the y-cordinates
rownames(reg1.clear.data) <- sapply(60:1, function(x) paste0("y", x)) 

# add the coordinates
reg1.clear.data <- cbind(sapply(60:1, function(x) paste0("y", x)), reg1.clear.data)

# rename the x-coordinates
colnames(reg1.clear.data) <- c("y", sapply(1:60, function(x) paste0("x", x)))

# data set with a row for each data point
reg1.clear.data.long <- as_tibble(reg1.clear.data) %>% pivot_longer(cols = starts_with("x"), values_to = "data") %>%  
  mutate(data = as.numeric(data), ID = row_number()) %>% rename(x = name) %>% 
  mutate(x = as.numeric(sub("x", "", x)), y = as.numeric(sub("y", "", y))) 

# standardaize the data with the standard deviation fo the complete data 
reg1.clear.data.long <- reg1.clear.data.long %>% mutate(data.sd = data/sd.all)

# amount of missing values
sum(is.na(reg1.clear.data.long$data))/(60*60)
#[1] 0.1730556

# Variogram estimation
varo.reg1.clear <- sim.it(reg1.clear.data.long$data.sd, cbind(reg1.clear.data.long$y, reg1.clear.data.long$x), hmax = c(4,4,3,3), cp = 1, missing = TRUE)


##########
### Graphs
##########

varos <- list(varo.reg1, varo.reg1.clear)

## lags
# S-N
lag.vec.SN <- cbind(rep(0, 4), 1:4)
lags.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x) %*% x))

# E-W
lag.vec.EW <- cbind(1:4, rep(0, 4))
lags.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x) %*% x))

# SE-NW
lag.vec.SENW <- cbind(1:3, -(1:3))
lags.SENW <- apply(lag.vec.SENW, 1, function(x) sqrt(t(x) %*% x))

# SW-NE
lag.vec.SWNE <- cbind(1:3, 1:3)
lags.SWNE <- apply(lag.vec.SWNE, 1, function(x) sqrt(t(x) %*% x))

lags <- list(lags.SN, lags.EW, lags.SENW, lags.SWNE)   

varo.long <- data.frame(lag = NA, direction = NA, estimator = NA, vario = NA, outlier = NA)
dic <- c("S-N", "E-W", "SW-NE", "SE-NW")
outlier <- c("with", "without")
est <- c("MCD.diff.re", "MCD.org.re", "Matheron", "Genton")
inde <- c(2,4,5,6)

# prepare for ggplot
for(d in 1:4){
  for(i in 1:2){
    for(e in 1:4){
      ls <- lags[[d]]
      varo.long <- rbind(varo.long, data.frame(lag = ls, direction = dic[d], estimator = est[e], vario = varos[[i]][[d]][,inde[e]], outlier = outlier[i])) 
    }
  }
}


# colours
cbbPalette1 <- c("#E69F00", "#56B4E9",  "#D55E00","#0072B2")

## direction S-N with and without Matheron (with outliers)
S_N <- filter(varo.long, direction == "S-N")
S_Ng <- ggplot(S_N, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator, linetype = outlier)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette1) +
  scale_shape_manual(values = c(3, 4, 17, 8)) + theme_minimal(base_size = 10)

S_N2 <- filter(S_N, estimator != "Matheron" & outlier == "with")
S_N2 <- rbind(S_N2, filter(S_N, outlier == "without"))
S_N2g <- ggplot(S_N2, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator, linetype = outlier)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette1) +
  scale_shape_manual(values = c(3, 4, 17, 8)) + theme_minimal(base_size = 10)
(S_Ng | S_N2g) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = 'bottom')
ggsave("Graphs/Anwendung_MCD_SN.pdf", width = 16.5, units = "cm") # Figure 10


## only robust estimators on the complete data to compare the different directions
cbbPalette2 <- c("#E69F00", "#D55E00","#0072B2")

S_N <- filter(varo.long, direction == "S-N" & outlier == "with" & estimator != "Matheron")
S_Ng <- ggplot(S_N, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette2) +
  theme_minimal(base_size = 10) + scale_shape_manual(values = c(3, 17, 8))


E_W <- filter(varo.long, direction == "E-W" & outlier == "with" & estimator != "Matheron")
E_Wg <- ggplot(E_W, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette2) +
  theme_minimal(base_size = 10)  + scale_shape_manual(values = c(3, 17, 8))

SW_NE <- filter(varo.long, direction == "SW-NE" & outlier == "with" & estimator != "Matheron")
SW_NEg <- ggplot(SW_NE, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette2) +
  theme_minimal(base_size = 10)  + scale_shape_manual(values = c(3, 17, 8))

SE_NW <- filter(varo.long, direction == "SE-NW" & outlier == "with" & estimator != "Matheron")
SE_NWg <- ggplot(SE_NW, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette2) +
  theme_minimal(base_size = 10)  + scale_shape_manual(values = c(3, 17, 8))


(S_Ng | E_Wg) / (SW_NEg | SE_NWg)  + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = 'bottom')
ggsave("Graphs/Anwendung_MCD.pdf", width = 16.5, units = "cm") # Figure 11
