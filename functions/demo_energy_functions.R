#                           dynapop-oft 
#       A population dynamics model accounting for the individual process 
#               of resource foraging on heterogeneous agricultural landscapes.
#     Copyright (C) 2017  Yoann BOURHIS at bourhis.yoann@gmail.com
#       see the corresponding publication at ... for details
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# About this source file:
# --------------------------------------------------------------------
#   functions need to set the non-spatial (i.e. demo-energetic) parameters
#    of the population dynamics model
# _____________________________________________________________________


library(limSolve)
library(ggplot2)
library(GGally)
library(nloptr)


#==================================================================================
#---    evaluate adequacy between population hypothesis and temporal dynamics   ---
#==================================================================================
eval.fn.light <- function(par, props, init=1, tend=1){
  
  res <- RK.model(par=par, pf=props[1], pl=props[2], init=init, t.end=tend)
  
  peak <- 200 + order(res$youngs.density[200:800], decreasing=TRUE)[1]
  
  # hypotheses :
  - min(res$youngs.density[peak], 2) + # maximise (until certain point) the population peak
    abs(res$youngs.density[1] - res$youngs.density[1000]) + # minimise the difference of yound density betwen start and end of season
    res$adults.density[1000]  # minimise the density of adults at the end of the season
  
}

#=============================================================
#--- Runge Kutta 4 solving of demography, no space             -------------
#=============================================================
RK4. <- function(t, species, par, pf, pl){
  # tau, mu, gb, gf, gl
  dJ <- - par[3] * pl * species[2] - # laying / birth rate
    par[1] * species[1] - # stage development (-)
    par[2] * species[1] * t # death rate
  
  dA <- par[1] * species[1] -# stage development (-)
    par[2] * species[2] * t # death rate
  
  c(dJ, dA)
}

#===================================================================
#--- main function of the non-spatial population dynamics model  ---
#===================================================================
RK.model <- function(par, pf=.05, pl=.1, init=1, t.end=1){
  
  # Dimension parameters
  time.step <- .001
  time <- seq(from = time.step, to = t.end, by = time.step)
  unit.t <- 1 / time.step # duration in hours of a time unit (t=1)
  
  # Initial conditions
  species <- c(init,0)
  youngs.density <- species[1]
  adults.density <- species[2]
  
  # Results and timer lists
  clockGlobal <- proc.time()[3]
  
  for (t in 1:length(time)){
    
    # Reaction RK4
    #______________________
    k1 <- RK4.(t=time[t], species=species, par=par, pf=pf, pl=pl) * time.step 
    k2 <- RK4.(t=time[t]+time.step/2, species=species+k1/2, par=par, pf=pf, pl=pl) * time.step 
    k3 <- RK4.(t=time[t]+time.step/2, species=species+k2/2, par=par, pf=pf, pl=pl) * time.step 
    k4 <- RK4.(t=time[t]+time.step, species=species+k3, par=par, pf=pf, pl=pl) * time.step 
    derivatives <- (k1 + 2*k2 + 2*k3 + k4)/6
    species <- species + derivatives
    
    youngs.density <- c(youngs.density, species[1])
    adults.density <- c(adults.density, species[2])
  }
  
  list(youngs.density=youngs.density, 
       adults.density=adults.density, 
       time=proc.time()[3] - clockGlobal)
  
}



#=============================================================
#    ABOUT ENERGY PARAMETERS
#-----------       Identify Gb,Gf,Gl,Gm          -------------
#=============================================================
EnergyParEst <- function(Gb, pf, pl, c, D, dt, ut, ux, resol, 
                         abs.max, plot=FALSE, sample.size=50000){
  
  D <- D * (ux^2 / ut)
  c <- c * (ux / ut)
  
  t <- 0 ; area <- 0
  # area explored during one time step
  while (area < 2){
    t <- t + 1
    r <- sqrt(4 * D * t)
    
    area <- (pi*r^2 + r*c*t - min(pi * r^2 / 2, r*c*t/2)) / resol^2
    
    if (c*t < r) area <- area - c * t * 2 * r / 2
    
  }
  
  if (c==0){
    ratio <- c(1, 0)
  }else{
    ratio <- c((pi * r^2 / resol^2) / (c*t/resol), 1) / 
      ((pi * r^2 / resol^2) / (c*t/resol) + 1)  
  }
  
  min <- c(-Gb, -abs.max, -abs.max) * t
  max <- c(abs.max, Gb, Gb/40) * t # flying cost always outreach BMR when at full speed
  
  E <- matrix(ncol = 3, byrow = TRUE, data = c(area*pf, area*pl, 1)) # global foraging is dE=0
  
  F <- c(-Gb)
  
  G <- matrix(ncol = 3, byrow = TRUE, data = c(1, 0, 0, 
                                               0, 1, 0,
                                               0, 0, 1,
                                               -1, 0, 0,
                                               0, -1, 0,
                                               0, 0, -1,
                                               -pf, -pl, -1/(area))) # local total foraging is dE<0
  
  H <- c(min, -max, Gb)
  
  xs <- xsample(E = E, F = F, G = G, H = H, iter = sample.size, type = "mirror", burninlength = sample.size/5)
  
  data <- xs$X ; colnames(data) <- c("g.feed", "g.lay", "g.mvt")
  g.mvt <- data[,"g.mvt"]
  # divide mvt cost in adv cost and diff cost, diff cost is added to gb because homogeneous
  data <- cbind(data, data[,"g.mvt"] * ratio[1], g.mvt, Gb * t)
  data[,"g.mvt"] <- data[,"g.mvt"] * ratio[2]
  colnames(data) <- c("g.feed", "g.lay", "g.mvt.adv", 'g.mvt.dif', "g.mvt", "g.base")
  data <- data.frame(data)
  
  if (plot){
    p <- ggpairs(rbind(data[which(1:nrow(data) %% 10 == 0),], min, max), 
                 lower=list(continuous=wrap("points", alpha = 0.1)))
    return(list(cbind(data, area) / t, p))
  } 
  
  cbind(data, area) / t
}

#=============================================================
#--- Multi local (MLSL+Nelder Mead Simplex) optimisation algorithm
#--- looking for demogrpahic parameters conditionnaly to the 
#---  estimated energy parameters
#=============================================================
demoParEst <- function(data, props, min=c(0,0), max=c(100,500), print=0, 
                       maxit=1000, LDEI=FALSE, plot=TRUE){
  
  if (LDEI){
    minGl <- maxGl <- data[order(apply(apply(data, 2, abs), 1, sum))[1],2]  
  }else{
    minGl <- min(data[,"g.lay"])
    maxGl <- max(data[,"g.lay"])
  } 
  
  par.bound <- cbind(tau=c(min[1],max[1]), mu=c(min[2],max[2]), Gl = c(minGl, maxGl))
  
  local_opts <- list( "algorithm" = "NLOPT_LN_NELDERMEAD",
                      maxeval=maxit/10,
                      ftol_abs=1e-3)
  
  opts <- list("algorithm" = "NLOPT_GN_MLSL_LDS",
               "local_opts" = local_opts,
               maxeval=maxit,
               print_level=print)
  
  OPT <- nloptr(
    x0 = runif(min=par.bound[1,], max=par.bound[2,], ncol(par.bound)),
    eval_f=eval.fn.light,
    lb = par.bound[1,],
    ub = par.bound[2,],
    eval_grad_f=NULL,
    opts=opts,
    props=props,
    init=1,
    tend=1) 
  
  print(OPT)  
  if (plot){
    # png('outputs/temp.png')
    res <- RK.model(par=c(OPT$solution), pf=props[1], pl=props[2], init=1, t.end=1)
    plot(res$adults.density,type='l', ylim=c(0, max(res[[1]], res[[2]])), 
         xlab='time-steps', ylab='Pop. density (adults:black & young:red)')
    lines(res$youngs.density,col="red")
    # dev.off()
  }
  
  who <- which(abs(data[,2] - OPT$solution[3]) < abs(OPT$solution[3]/100))
  fitted <- c(OPT$solution[-3], data[who[order(data[who,"g.mvt"], decreasing=TRUE)[1]],])
  names(fitted)[1:2] <- c("tau", "mu")
  fitted
  
}

