################################################################################
#### Graphics and tables for the paper                                      ####
################################################################################

################
## Load packages
################

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(xtable)

##########
## Colours
##########

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73",  "#D55E00", "#CC79A7",  "#0072B2","#000000", "#F0E442")

##############
#### Chapter 4
##############

############
## Figure 1: True Variogram

sp.var <- function(h, a = 5, c = 1){
  if(h < a){res <- c*(((3*h)/(2*a)) - (1/2) * (h/a)^3)}
  if(h >= a){res <- c}
  if(h == 0){res <- 0}
  return(res)
}

# Anisotropie parameters
Ts <- diag(c(1, 1/2))
R <- matrix(c(cos(3*pi/8), -sin(3*pi/8), sin(3*pi/8), cos(3*pi/8)), nrow = 2)

# lag-vectors
# S-N
lag.vec.SN <- cbind(rep(0, 15), 1:15)
lag.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x)  %*% x))
lags.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
    
# E-W
lag.vec.EW <- cbind(1:15, rep(0, 15))
lag.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x)  %*% x))
lags.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))

# SW-NE
lag.vec.SWNE <- cbind(1:15, 1:15)
lag.SWNE <- apply(lag.vec.SWNE, 1, function(x) sqrt(t(x)  %*% x))
lags.SWNE <- apply(lag.vec.SWNE, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
    
# SE-NW
lag.vec.SENW <- cbind(1:15, -(1:15))
lag.SENW <- apply(lag.vec.SENW, 1, function(x) sqrt(t(x)  %*% x))
lags.SENW <- apply(lag.vec.SENW, 1, function(x) sqrt(t(x) %*% t(R) %*% t(Ts) %*% Ts %*% R %*% x))
  
# true model
true.var <- list()
true.var[[1]] <- sapply(lags.SN, function(l) 2*sp.var(l, a = 5, c = 1))
true.var[[2]] <- sapply(lags.EW, function(l) 2*sp.var(l, a = 5, c = 1))
true.var[[3]] <- sapply(lags.SWNE, function(l) 2*sp.var(l, a = 5, c = 1))
true.var[[4]] <- sapply(lags.SENW, function(l) 2*sp.var(l, a = 5, c = 1))


true.long <- data.frame("h" = c(lag.SN, lag.EW, lag.SWNE, lag.SENW), 
                        "vario" = c(true.var[[1]], true.var[[2]], true.var[[3]], true.var[[4]]),
                        "direction" = c(rep("S-N", 15), rep("E-W", 15), rep("SW-NE", 15), rep("SE-NW",15)))

plot.sph <- ggplot(true.long, aes(x = h, y = vario, col = direction, shape = direction)) + geom_line(linewidth = 0.25) + geom_point(size= 2) +
  xlab("||h||") + ylab(expression(2 * gamma(h))) + scale_colour_manual(values=cbbPalette) +
  theme_minimal(base_size = 18) + xlim(c(0,12))
plot.sph
ggsave("Graphs/sph_var_true.pdf")

#########################
## 4.1: Correctionfactors
#########################

###########
## Table 1: 15 X 15 grid, all directions, all estimators
load(file = "Correctionfactors/correctionfactors.RData")

xtable(round(res.corr[,,1, 1], 3))

####################################
## 4.3: Non-contaminated normal data
####################################

############
## Figure 2: Bias, all estimators, all directions, 15 x 15 grid
load("Consistency/estimation_bias.RData")

# 15 x 15 grid
bias <- kons.bias[, , , 1]

# prepare data for ggplot
lag.vec.1 <- cbind(1:7, rep(0, 7))
lags.1 <- apply(lag.vec.1, 1, function(x) sqrt(t(x)%*% x))
lags.1 <- c(lags.1, NA, NA)


lag.vec.2 <- cbind(1:5, -(1:5))
lags.2 <- apply(lag.vec.2, 1, function(x) sqrt(t(x) %*% x))
lags.2 <- c(lags.2, NA, NA, NA, NA)


est <- c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re", "Matheron", "Genton")

dic <- c("S-N", "E-W", "SW-NE", "SE-NW")


bias.long <- data.frame(estimator = NA, direction = NA, lag = NA, bias = NA)
for(s in 1:6){
  for(r in 1:4){
    for(l in 1:7){
      if(r == 1 | r == 2){
        lags <- lags.1
      } else{ lags <- lags.2}

      bias.long <- rbind(bias.long, data.frame(estimator = est[s], direction = dic[r], lag = lags[l], bias = bias[s, r, l]))
    }
  }
}


E_W <- filter(bias.long, direction == "E-W")
Bias_E_W <- ggplot(data = E_W, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +  geom_line(linewidth = 0.25) + 
  geom_point(size= 2) + ylab("Bias") + xlab("||h||") + scale_colour_manual(values=cbbPalette) +
        scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(-0.6, 0.15)
Bias_E_W


S_N <- filter(bias.long, direction == "S-N")
Bias_S_N <- ggplot(data = S_N, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +  geom_line(linewidth = 0.25) + 
  geom_point(size= 2) + ylab("Bias") + xlab("||h||")  + scale_colour_manual(values=cbbPalette) +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(-0.6, 0.15)
Bias_S_N

SW_NE<- filter(bias.long, direction == "SW-NE")
Bias_SW_NE <- ggplot(data = SW_NE, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +  geom_line(linewidth = 0.25) +  
  geom_point(size= 2) + ylab("Bias") + xlab("||h||")  + scale_colour_manual(values=cbbPalette) +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(-0.6, 0.15)
Bias_SW_NE


SE_NW<- filter(bias.long, direction == "SE-NW")
Bias_SE_NW <- ggplot(data = SE_NW, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +  geom_line(linewidth = 0.25) +  
  geom_point(size= 2) + ylab("Bias") + xlab("||h||")  + scale_colour_manual(values=cbbPalette) +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(-0.6, 0.15)
Bias_SE_NW


(Bias_E_W | Bias_S_N)/(Bias_SE_NW| Bias_SW_NE) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")& theme(legend.position = 'bottom')
ggsave("Graphs/Bias_oA.pdf", width = 16.5, units = "cm")


############
## Figure 3: rMSE, all estimators, all directions, 15 x 15 grid
load("Consistency/estimation_sqrtMSE.RData")

# variogram 1, gridsize 1
sqrtMSE <- kons.sqrtMSE[,,,1]

# prepare data for ggplot
lag.vec.1 <- cbind(1:7, rep(0, 7))
lags.1 <- apply(lag.vec.1, 1, function(x) sqrt(t(x)%*% x))
lags.1 <- c(lags.1, NA, NA)


lag.vec.2 <- cbind(1:5, -(1:5))
lags.2 <- apply(lag.vec.2, 1, function(x) sqrt(t(x) %*% x))
lags.2 <- c(lags.2, NA, NA, NA, NA)


est <- c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re", "Matheron", "Genton")

dic <- c("S-N", "E-W", "SW-NE", "SE-NW")


sqrtMSE.long <- data.frame(estimator = NA, direction = NA, lag = NA, sqrtMSE = NA)
for(s in 1:6){
  for(r in 1:4){
    for(l in 1:7){
      if(r == 1 | r == 2){
        lags <- lags.1
      } else{ lags <- lags.2}

      sqrtMSE.long <- rbind(sqrtMSE.long, data.frame(estimator = est[s], direction = dic[r], lag = lags[l], sqrtMSE = sqrtMSE[s, r, l]))
    }
  }
}


E_W <- filter(sqrtMSE.long, direction == "E-W")
sqrtMSE_E_W <- ggplot(data = E_W, mapping = aes(x = lag, y = sqrtMSE, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(sqrt("MSE"))) +
  xlab("||h||") + scale_colour_manual(values=cbbPalette) +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(0, 1.3)
sqrtMSE_E_W


S_N <- filter(sqrtMSE.long, direction == "S-N")
sqrtMSE_S_N <- ggplot(data = S_N, mapping = aes(x = lag, y = sqrtMSE, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(sqrt("MSE"))) +
  xlab("||h||") + scale_colour_manual(values=cbbPalette) +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(0, 1.3)
sqrtMSE_S_N

SE_NW <- filter(sqrtMSE.long, direction == "SE-NW")
sqrtMSE_SE_NW <- ggplot(data = SE_NW, mapping = aes(x = lag, y = sqrtMSE, col = estimator, shape = estimator), shape = estimator) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(sqrt("MSE"))) +
  xlab("||h||") + scale_colour_manual(values=cbbPalette) +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(0, 1.3)
sqrtMSE_SE_NW

SW_NE <- filter(sqrtMSE.long, direction == "SW-NE")
sqrtMSE_SW_NE <- ggplot(data = SW_NE, mapping = aes(x = lag, y = sqrtMSE, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(sqrt("MSE"))) +
  xlab("||h||") + scale_colour_manual(values=cbbPalette) +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(0, 1.3)
sqrtMSE_SW_NE


(sqrtMSE_E_W | sqrtMSE_S_N)/(sqrtMSE_SE_NW| sqrtMSE_SW_NE) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = 'bottom')
ggsave("Graphs/MSE_oA.pdf", width = 16.5, units = "cm")


###################
## 4.3: Consistency
###################

############
## Figure 4: Boxplots, 1 direction (E-W), 4 lags, 4 grids, incl. modified estimators
load(file = "Modified/estimation_errors.RData")


errors <- kons.est.error[,,c(1,4,7),,]


# prepare data for ggplot (for Figure 4 and Figure 5)
est <- c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re", "Matheron", "Genton",
         "MCD.diff.mod", "MCD.diff.mod.re", "MCD.org.mod", "MCD.org.mod.re")

it <- 1:1000

lags <- c(1, 4, 7)

grid <- c("15x15", "25x25", "50x50", "75x75")

dic <- c("S-N", "E-W")

errors_long <- data.frame(estimator = NA, iteration = NA, lag = NA, direction = NA, grid = NA, error = NA)
for(s in 1:10){
  for(r in 1:2){
    for(l in 1:3){
      for(g in 1:4){  
        errors_long <- rbind(errors_long, data.frame(estimator = est[s], iteration = it, direction = dic[r], lag = lags[l], grid = grid[g], error = errors[s, r, l,, g]))
      }
    }
  }
}

errors_long$estimator <- factor(errors_long$estimator, levels = c("Matheron", "Genton", "MCD.diff", "MCD.diff.mod", "MCD.diff.re", "MCD.diff.mod.re", "MCD.org", "MCD.org.mod", "MCD.org.re", "MCD.org.mod.re"))


# E-W
E_W_l1 <- filter(errors_long, direction == "E-W" & lag == 1)
boxplot_EW_l1 <- ggplot(E_W_l1, aes(x = estimator, y = error, fill = grid)) + geom_boxplot() + ylab("estimation error") + ggtitle(expression(paste("(1,0)", )^{T})) + 
  scale_colour_manual(values=cbbPalette) + theme_minimal(base_size = 10)

E_W_l4 <- filter(errors_long, direction == "E-W" & lag == 4)
boxplot_EW_l4 <- ggplot(E_W_l4, aes(x = estimator, y = error, fill = grid)) + geom_boxplot() + ylab("estimation error") + ggtitle(expression(paste("(4,0)", )^{T})) + 
  scale_colour_manual(values=cbbPalette) + theme_minimal(base_size = 10)

E_W_l7 <- filter(errors_long, direction == "E-W" & lag == 7)
boxplot_EW_l7 <- ggplot(E_W_l7, aes(x = estimator, y = error, fill = grid)) + geom_boxplot() + ylab("estimation error") + ggtitle(expression(paste("(7,0)", )^{T})) + 
  scale_colour_manual(values=cbbPalette) + theme_minimal(base_size = 10)

boxplot_EW_l1 / boxplot_EW_l4 / boxplot_EW_l7 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("Graphs/Boxplots_EW.pdf", width = 22, height = 22, units = "cm")


############
## Figure 5: Boxplots, 1 direction (S-N), 1 lag, 4 grids, incl. modified estimators

# S-N
S_N_l7 <- filter(errors_long, direction == "S-N" & lag == 7)
boxplot_SN_l7 <- ggplot(S_N_l7, aes(x = estimator, y = error, fill = grid)) + geom_boxplot() + ylab("estimation error") + ggtitle(expression(paste("(0,7)", )^{T})) + 
  scale_colour_manual(values=cbbPalette)  + theme_minimal(base_size = 10)

boxplot_SN_l7 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("Graphs/Boxplot_SN.pdf", width = 22, height = 10, units = "cm")



#########################################
## 4.4: Contamination with block outliers
#########################################

############
## Figure 6: Bias for E-W, S-N, 5% und 15%, 25%, N(3,1) 
load("Block/estimation_bias.RData")

# 15 x 15 grid, N(3,1),  5% und 15%
bias <- kons.bias[,1:2,,c(1,3, 4),2]

# prepare data for ggplot
lag.vec.1 <- cbind(1:7, rep(0, 7))
lags.1 <- apply(lag.vec.1, 1, function(x) sqrt(t(x)%*% x))
lags.1 <- c(lags.1, NA, NA)

est <- c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re", "Matheron", "Genton")

dic <- c("S-N", "E-W")

am <- c("0.05", "0.15", "0.25")

bias.long <- data.frame(estimator = NA, direction = NA, lag = NA, amount = NA, bias = NA)
for(s in 1:6){
  for(r in 1:2){
    for(l in 1:7){
      for(a in 1:3){
      lags <- lags.1

      bias.long <- rbind(bias.long, data.frame(estimator = est[s], direction = dic[r], lag = lags[l], amount = am[a], bias = bias[s, r, l, a]))
      }
    }
  }
}

E_W_0.05 <- filter(bias.long, direction == "E-W", amount == "0.05" | amount == "0")
Bias_E_W_0.05 <- ggplot(data = E_W_0.05, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") +
  labs(title = "a") +  scale_colour_manual(values=cbbPalette) + xlab("||h||") +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10)  + ylim(-0.5, 1.1)
Bias_E_W_0.05


S_N_0.05 <- filter(bias.long, direction == "S-N", amount == "0.05"| amount == "0")
Bias_S_N_0.05 <- ggplot(data = S_N_0.05, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") +
  labs(title = "b") + scale_colour_manual(values=cbbPalette) + xlab("||h||") +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10)  + ylim(-0.5, 1.1)
Bias_S_N_0.05

E_W_0.15 <- filter(bias.long, direction == "E-W", amount == "0.15"| amount == "0")
Bias_E_W_0.15 <- ggplot(data = E_W_0.15, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + scale_colour_manual(values=cbbPalette) + theme_minimal(base_size = 10)
Bias_E_W_0.15


S_N_0.15 <- filter(bias.long, direction == "S-N", amount == "0.15"| amount == "0")
Bias_S_N_0.15 <- ggplot(data = S_N_0.15, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10)
Bias_S_N_0.15


E_W_0.25 <- filter(bias.long, direction == "E-W", amount == "0.25"| amount == "0")
Bias_E_W_0.25 <- ggplot(data = E_W_0.25, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + scale_colour_manual(values=cbbPalette) + theme_minimal(base_size = 10)
Bias_E_W_0.25


S_N_0.25 <- filter(bias.long, direction == "S-N", amount == "0.25"| amount == "0")
Bias_S_N_0.25 <- ggplot(data = S_N_0.25, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) + 
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10)
Bias_S_N_0.25


(Bias_E_W_0.05 | Bias_S_N_0.05)/(Bias_E_W_0.15| Bias_S_N_0.15 )/(Bias_E_W_0.25| Bias_S_N_0.25 ) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("Graphs/Bias_BA.pdf", width = 16.5, height = 16.5, units = "cm")



############
## Figure 7: Bias, SW-NE & SE-NW, 10%, N(3, 1), N(5,1), N(0,4)

load("Block/estimation_bias.RData")

# 15 x 15 grid, N(3, 1), N(5,1), N(0,4), 10%
bias <- kons.bias[,3:4,,2,c(2, 3, 4)]

# prepare data for ggplot
lag.vec.1 <- cbind(1:5, 1:5)
lags.1 <- apply(lag.vec.1, 1, function(x) sqrt(t(x)%*% x))


est <- c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re", "Matheron", "Genton")

dic <- c("SW-NE", "SE-NW")

dis <- c("N(3,1)", "N(5,1)", "N(0,4)")

bias.long <- data.frame(estimator = NA, direction = NA, lag = NA, dist = NA, bias = NA)
for(s in 1:6){
  for(r in 1:2){
    for(l in 1:5){
      for(d in 1:3){
        lags <- lags.1

        bias.long <- rbind(bias.long, data.frame(estimator = est[s], direction = dic[r], lag = lags[l], dist = dis[d], bias = bias[s, r, l, d]))
      }
    }
  }
}

SW_NE_1 <- filter(bias.long, direction == "SW-NE", dist == "N(3,1)")
Bias_SW_NE_1 <- ggplot(data = SW_NE_1, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") +
  labs(title = "a") + scale_colour_manual(values=cbbPalette) + xlab("||h||") +
  scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) + theme_minimal(base_size = 10) + ylim(-0.5,2.5)
Bias_SW_NE_1

SW_NE_2 <- filter(bias.long, direction == "SW-NE", dist == "N(5,1)")
Bias_SW_NE_2 <- ggplot(data = SW_NE_2, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10) + ylim(-0.5,5.5)
Bias_SW_NE_2

SW_NE_3 <- filter(bias.long, direction == "SW-NE", dist == "N(0,4)")
Bias_SW_NE_3 <- ggplot(data = SW_NE_3, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") + 
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10) + ylim(-0.5,4)
Bias_SW_NE_3

SE_NW_1 <- filter(bias.long, direction == "SE-NW", dist == "N(3,1)")
Bias_SE_NW_1 <- ggplot(data = SE_NW_1, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  labs(title = "b") + scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10) + ylim(-0.5, 2.5)
Bias_SE_NW_1

SE_NW_2 <- filter(bias.long, direction == "SE-NW", dist == "N(5,1)")
Bias_SE_NW_2 <- ggplot(data = SE_NW_2, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10) + ylim(-0.5,5.5)
Bias_SE_NW_2

SE_NW_3 <- filter(bias.long, direction == "SE-NW", dist == "N(0,4)")
Bias_SE_NW_3 <- ggplot(data = SE_NW_3, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10) + ylim(-0.5,4)
Bias_SE_NW_3

(Bias_SW_NE_1 | Bias_SE_NW_1)/(Bias_SW_NE_2| Bias_SE_NW_2)/(Bias_SW_NE_3| Bias_SE_NW_3) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("Graphs/Bias_BA_2.pdf", width = 16.5, height = 16.5, units = "cm")


############################################
## 4.5: Contamination with isolated outliers
############################################

############
## Figure 8: Bias, E-W direction, 5 % , 15% , 25 %, N(3,1), N(5,1)

load("Isolated/estimation_bias.RData")

# 15 x 15 grid, N(3, 1), N(5,1), 5%, 15%, 25%, E-W
bias <- kons.bias[,2,,c(1, 3, 4),c(2, 3)]

# prepare data for ggplot
lag.vec.1 <- cbind(1:7, rep(0,7))
lags.1 <- apply(lag.vec.1, 1, function(x) sqrt(t(x)%*% x))

est <- c("MCD.diff", "MCD.diff.re", "MCD.org", "MCD.org.re", "Matheron", "Genton")

dis <- c("N(3,1)", "N(5,1)")

am <- c("0.05", "0.15", "0.25")

bias.long <- data.frame(estimator = NA, lag = NA, dist = NA, amount = NA, bias = NA)
for(s in 1:6){
  for(a in 1:3){
    for(l in 1:7){
      for(d in 1:2){
        lags <- lags.1

        bias.long <- rbind(bias.long, data.frame(estimator = est[s], lag = lags[l], dist = dis[d], amount = am[a], bias = bias[s, l, a, d]))
      }
    }
  }
}

E_W_d1_0.05 <- filter(bias.long,  dist == "N(3,1)", amount == "0.05")
Bias_E_W_d1_0.05 <- ggplot(data = E_W_d1_0.05, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  labs(title = "a") + scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10)
Bias_E_W_d1_0.05

E_W_d1_0.15 <- filter(bias.long,  dist == "N(3,1)", amount == "0.15")
Bias_E_W_d1_0.15 <- ggplot(data = E_W_d1_0.15, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10)
Bias_E_W_d1_0.15

E_W_d1_0.25 <- filter(bias.long,  dist == "N(3,1)", amount == "0.25")
Bias_E_W_d1_0.25 <- ggplot(data = E_W_d1_0.25, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10)
Bias_E_W_d1_0.25

E_W_d2_0.05 <- filter(bias.long,  dist == "N(5,1)", amount == "0.05")
Bias_E_W_d2_0.05 <- ggplot(data = E_W_d2_0.05, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  labs(title = "b") + scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10)
Bias_E_W_d2_0.05

E_W_d2_0.15 <- filter(bias.long,  dist == "N(5,1)", amount == "0.15")
Bias_E_W_d2_0.15 <- ggplot(data = E_W_d2_0.15, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10)
Bias_E_W_d2_0.15

E_W_d2_0.25 <- filter(bias.long,  dist == "N(5,1)", amount == "0.25")
Bias_E_W_d2_0.25 <- ggplot(data = E_W_d2_0.25, mapping = aes(x = lag, y = bias, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab("Bias") + xlab("||h||") +
  scale_colour_manual(values=cbbPalette) + scale_shape_manual(values = c(3, 4, 19, 17, 15, 8)) +
  theme_minimal(base_size = 10)
Bias_E_W_d2_0.25

(Bias_E_W_d1_0.05| Bias_E_W_d2_0.05)/(Bias_E_W_d1_0.15| Bias_E_W_d2_0.15)/(Bias_E_W_d1_0.25| Bias_E_W_d2_0.25) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("Graphs/Bias_IA.pdf", width = 16.5, height = 16.5, units = "cm")


###################################
## 5: Application to Satellite Data
###################################

## see document: Application-Satellite-Data

###########
## Appendix
###########

###########
## Table 2: 15 X 15 grid, all directions, all estimators, exponential and gaussian variogram
load(file = "Correctionfactors/correctionfactors.RData")

xtable(round(res.corr[,,1, 2], 3))
xtable(round(res.corr[,,1, 3], 3))

###########
## Table 2: Bias & rMSE multiplied with 10, 3 lags, 2 distributions, 2 amounts, direction E-W & S-N
##          GRF with outlier block

load("Block/estimation_bias.RData")
load("Block/estimation_sqrtMSE.RData")

# 15 x 15 grid, N(3,1) & N(0,4),  5% and 15%, S-N  & E-W, lags c(1,4,7)
bias <- kons.bias[,1:2,c(1,4,7),c(1,3),c(2,4)]
rMSE <- kons.sqrtMSE[,1:2,c(1,4,7),c(1,3),c(2,4)]

bias <- round(bias * 10, 2)
rMSE <- round(rMSE * 10, 2)

# modified estimators
load("Modified/estimation_bias_block.RData")
load("Modified/estimation_sqrtMSE_block.RData")

bias <- kons.bias[,,c(1,4,7),,]
rMSE <- kons.sqrtMSE[,,c(1,4,7),,]

bias <- round(bias * 10, 2)
rMSE <- round(rMSE * 10, 2)

