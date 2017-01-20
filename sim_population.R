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
# ---------------------------------------------
#   simulate one season of population dynamics
#     on the r.map landscape resulting from
#       script gen_ladscape.R
# _____________________________________________


# GLOBAL (dicretisation) PARAMETERS
x.end <- y.end <- 1 # non-dimensionalised space extent
resol <- 10 # side dimension in meter of space cells
space.dim.x <- space.dim.y <- 800 # extent in meter of the spatial domain
dx <- x.end / (space.dim.x/resol) # non-dimensionalised cell side size
dy <- y.end / (space.dim.y/resol)
t.end <- 1 # non-dimensionalised time extent
t.resol <- 1 # duration in hour of a time step
time.dim <- 1000 # duration in hours of the simulation
time.step <- t.end / (time.dim/t.resol) # non-dimensionalised time step
e.max <- 30 # number of discretisation in energy dimension
# diffusion parameter :
diff <- .01


# actual LCP (proportions) after rasterisation
props <- c(sum(as.matrix(r.map)==2)/prod(ncells), # feed 
           sum(as.matrix(r.map)==1)/prod(ncells)) # lay

#  FIRST STEP :
#     IDENTIFYING A RELEVANT DEMO-ENERGY PARAMETERS SET
# _______________________________________________________
source('functions/demo_energy_functions.R')
g.base <- -0.00833 * 1000 # from E=1 to E=0 in 5 days (Quite wrong, but working)
abs.max <- -g.base * 48 # at least 1 hour/day feeding to maintain BMR

pdf(file='outputs/demoParPlot.pdf', width=20, height=15)
demo.sets <- EnergyParEst(Gb=g.base, abs.max=abs.max, c=1, D=diff, plot=TRUE, 
                          pf=landscape.res$metrics.res['LCP_f'], 
                          pl=landscape.res$metrics.res['LCP_l'], 
                          dt=time.step, resol=resol, 
                          ut=time.dim/t.end, ux=space.dim.x/x.end)

demo.sets[[1]][order(apply(apply(demo.sets[[1]], 2, abs), 1, sum))[1],]
demo.sets[[2]]

fitted <- demoParEst(data=demo.sets[[1]], 
                     props=landscape.res$metrics.res[c('LCP_f', 'LCP_l')],
                     #min=c(1,.1), 
                     max=c(100,100), 
                     print=0, maxit=1000, LDEI=FALSE, plot=TRUE)
print(fitted)

# convert par names to fit FP or adv xp.
fitted[['g.base']] <- fitted[['g.base']] + fitted[['g.mvt.dif']]
fitted[['g.mvt']] <- fitted[['g.mvt.adv']]
print(t(fitted))

parametersDF <- data.frame(matrix(fitted, nrow=1))
colnames(parametersDF) <- names(fitted)
plot(tableGrob(round(unlist(parametersDF), 5)))
dev.off()

#  SECOND STEP :
#     SPATIAL SIMULATION !
# _______________________________________________________
source('functions/RAD_function.R')
source('functions/schemes_functions.R')

pars <- c(resol=resol, t.end=t.end, time.step=time.step, e.max=e.max, cfl=.8, e.init=1/3, 
          dif=diff, max.speed=42, K=1, fitted, f.l.ratio=1, 
          save.freq=25,
          percept.dist=100, 
          s.lay=5,
          s.feed=5)

simulation <- RAD.model(parameters=pars, maps=list(r.init, r.map), plot=TRUE)
