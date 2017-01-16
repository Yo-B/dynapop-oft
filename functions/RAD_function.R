#                           dynapop-oft 
#       A population dynamics model accounting for the individual process 
#               of resource foraging on heterogeneous agricultural landscapes.
#     Copyright (C) 2017  Yoann BOURHIS at bourhis.yoann@gmail.com
#       see the associated publication at ... for details
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
#   Main function of the population dynamics model
#   relies on schemes_functions.R 
# _____________________________________________________________________

RAD.model <- function(parameters=NULL, maps=maps, pad=c(0,0), verb=1, plot=FALSE, 
                      exportOnly=FALSE, save.freq=10){
  
  if(!is.null(parameters)){
    mapply(assign, names(parameters), parameters, pos=1)
  }
  par.demo <- c(-g.lay, tau, mu, K)
  
  print(parameters)
  
  if (!exists('percept.dist')){ # perception distance can be resource-dependent
    percept.dist <- c(percept.dist.lay, percept.dist.feed)
  }
  
  # Extracting maps
  matrices <- maps
  init.map <- padding.map(matrix(matrices[[1]]@data@values, 
                                 nrow=matrices[[1]]@nrows, byrow=TRUE), pad[1], pad[2])
  resources.map <- padding.map(matrix(matrices[[2]]@data@values, 
                                      nrow=matrices[[2]]@nrows, byrow=TRUE), pad[1], pad[2])
  # slightly bluring laying sites helps the resolution (regularising the map)
  blur.laying <- as.matrix(blur(as.im(resources.map==1), sigma=1))
  blur.laying[blur.laying<0] <- 0
  blur.laying <- blur.laying * (sum(resources.map==1)/sum(blur.laying))
  # # initial condition can be heterogeneous, if so some slight blur secure stability at resolution start
  # blur.init <- as.matrix(blur(as.im(init.map), sigma=1))
  # blur.init[blur.init<0.01] <- 0
  # blur.init <- blur.init * (sum(init.map==1)/sum(blur.init))
  # # try homogeneous initial conditions
  # blur.init <- blur.init * 0 + .1
  
  # Dimension parameters
  ny <- dim(resources.map)[2] ; nx <- dim(resources.map)[1]
  dx <- 1/nx ; dy <- 1/ny
  dE <- 1/e.max
  time <- seq(from = time.step, to = t.end, by = time.step)
  unit.t <- 1 / time.step # duration in hours of a time unit (t=1)
  
  # Advective fields (and CFL approved time step)
  
  fields <- gen.adv.fields(d=percept.dist, dx=dx, dy=dy, Nx=nx, Ny=ny,
                           dt=time.step, resol=resol, map=resources.map)
  
  adv.cubes <- gen.adv.cubes(fields, Emax=e.max, max.speed=max.speed, norm=TRUE,
                             dt=time.step, resol=resol, dx=dx, dy=dy,
                             func.poly=function(x) exp(-s.lay*x), 
                             func.edges=function(x) f.l.ratio * exp(-s.feed*x)) ##################
  
  
  e.start <- round(e.max * e.init) # where new adults are spawning in E dimension
  
  # Energy LS is divided into energy.ls.cte, depending on land-covers only
  energy.ls.cte <- resources.map * 0 + g.base # accounting for blur laying
  energy.ls.cte[resources.map==2] <- g.feed + g.base
  # and energy.ls.mvt accounting for mvt cost and varying along E
  energy.ls.mvt<- array(unlist(lapply(1:e.max, FUN=function(x) energy.ls.cte * 0 + 
                                        g.mvt * sqrt(((adv.cubes[[1]][,,x]+
                                                         adv.cubes[[2]][,,x])/2)^2 +
                                                       ((adv.cubes[[3]][,,x]+
                                                           adv.cubes[[4]][,,x])/2)^2))), 
                        dim=c(dim(resources.map), e.max))
  energy.ls.cte <- array(rep(energy.ls.cte, e.max), dim=c(nx,ny,e.max))
  # finite volumes upwind cubes
  e.ls.p <- abind(energy.ls.mvt[,,-e.max], matrix(0,nx,ny), along=3) + 
    abind(energy.ls.cte[,,-e.max], matrix(0,nx,ny), along=3) # accumulating boundary
  e.ls.m <- abind(energy.ls.mvt[,,1], energy.ls.mvt[,,-e.max], along=3) + energy.ls.cte # dirichlet condition
  # laying cost evolve depending on field occupation, see in loop
  
  # adjusting CFL compliant sub time-steps
  it.adv.XY <- sapply(1:e.max, FUN=function(x)
    ceiling(time.step / (min(dx,dy)/max(max(abs(adv.cubes[[1]][,,x])),
                                        max(abs(adv.cubes[[2]][,,x])),
                                        max(abs(adv.cubes[[3]][,,x])),
                                        max(abs(adv.cubes[[4]][,,x])))*cfl))) # depend on slices
  dt.adv.XY <- time.step / it.adv.XY
  
  # compute cfl according to varying laying cost
  max.adv.E <- max(abs(cbind(e.ls.m + rep(blur.laying * g.lay, e.max), 
                             e.ls.p + c(rep(blur.laying * g.lay, e.max-1), rep(0, nx*ny)))))
  
  it.adv.E <- ceiling(max.adv.E*time.step/(1/e.max))
  dt.adv.E <- time.step / pmax(it.adv.E, 1)
  
  
  # Diffusion map, possibly heterogeneous
  DifMap <- resources.map * 0 + dif # in m2/h
  DifMap <- DifMap * unit.t * dx / resol # in dx^2/dt
  # ADI matrices
  matrix.list <- gen.adi.matrix(dt = time.step, D=DifMap, BC='neumann')[1:4]
  
  
  # Initial conditions
  species <- c(init.map, init.map * 0)
  species <- array(species, dim = c(nx, ny, 2))
  youngs.density <- species[,,1]
  adults.density <- species[,,2]
  species.cube <- array(0, dim=c(nx,ny,e.max)) ####################
  species[,,2] <- apply(species.cube, FUN=sum, c(1,2))
  
  # Results and timer lists
  RES <- apply(species, FUN = sum, 3)
  timeTotal <- NULL
  clockGlobal <- proc.time()[3]
  
  # debug option
  if(exportOnly) return(list(parameters=parameters, adv.cubes=adv.cubes, 
                             energy.ls.cte=energy.ls.cte, energy.ls.mvt=energy.ls.mvt, 
                             fields=ifelse(exists('fields'), fields, NA), 
                             it.adv.E=it.adv.E, dt.adv.E=dt.adv.E, dt.adv.E=dt.adv.E, 
                             it.adv.E=it.adv.E, it.adv.XY=it.adv.XY, dt.adv.XY=dt.adv.XY,
                             kern.poly=ifelse(exists('fields'), fields[[length(fields)-1]], NA), 
                             kern.edges=ifelse(exists('fields'), fields[[length(fields)]], NA),
                             DifMap=ifelse(exists('DifMap'), DifMap, NA)))
  
  # delete intermediate heavy objects from memory
  rm(list=c('fields', 'resources.map', 'DifMap', 'energy.ls.cte', 'energy.ls.mvt',
            'init.map')) ; gc()
  
  # display init status
  if (verb > 0) {
    cat('\nsum(Y) = ', sum(species[,,1]), '   min(Y) = ', min(species[,,1]), '   max(Y) = ', max(species[,,1]))
    cat('\nsum(A) = ', sum(species[,,2]), '   min(A) = ', min(species[,,2]), '   max(A) = ', max(species[,,2]))
    cat('\nsum(A+Y) = ', sum(species[,,1])+sum(species[,,2]))
    cat('\ntime-step ', 0,' for date ', 0)
    if (plot){
      par(mfrow=c(2,2))
      plot(maps[[2]]) ; hist(apply(species.cube, 3, sum), e.max)
      plot(raster(species[,,1])) ; plot(raster(species[,,2])) 
    }
  }
  df <- data.frame(time=c(0, time), Y=c(sum(species[,,1]), rep(NA, length(time))), A=c(0, rep(NA, length(time))))
  status <- 'ok'
  #=============================================================
  #-----------           Resolution loop           -------------
  #=============================================================
  cat('\nStarting resolution ...\n')
  
  for (t in 1:length(time)){
    
    clockTotal <- proc.time()[3]  
    temp <- sum(species.cube)
    
    if(sum(is.nan(species)) >= 1){
      cat('\nExiting simulation with status NaNError')
      status <- 'NaNError'
      break
    }
    
    # Reaction RK4
    #______________________
    k1 <- RK4(t=time[t], species=species, par.demo=par.demo, 
              laying=blur.laying, Nx=nx, Ny=ny) * time.step 
    k2 <- RK4(t=time[t]+time.step/2, species=species+k1/2, par.demo=par.demo, 
              laying=blur.laying, Nx=nx, Ny=ny) * time.step 
    k3 <- RK4(t=time[t]+time.step/2, species=species+k2/2, par.demo=par.demo, 
              laying=blur.laying, Nx=nx, Ny=ny) * time.step 
    k4 <- RK4(t=time[t]+time.step, species=species+k3, par.demo=par.demo, 
              laying=blur.laying, Nx=nx, Ny=ny) * time.step
    derivatives <- (k1 + 2*k2 + 2*k3 + k4)/6
    evo <- derivatives[,,2]/species[,,2] + 1
    species <- species + derivatives
    # dispatching densities across the energy dimension
    evo[is.na(evo)] <- 1 
    evo[evo==Inf] <- 1
    species.cube <- array(apply(species.cube, 3, function(x) x * pmin(evo,1)), dim=c(nx,ny,e.max))
    species.cube[,,e.start] <- species.cube[,,e.start] + pmax(derivatives[,,2],0)    
    
    t.rk <- proc.time()[3] - clockTotal
    t.dif <- proc.time()[3]
    
    dRK4 <- temp-sum(species.cube) # check flows
    dRK4 <- ifelse(is.na(dRK4), 'Na(N)ERROR\n', 'OK')
    temp <- sum(species.cube)
    # Diffusion 
    # _____________________
    species.cube <- diffusion3d(cube=species.cube, Emax=e.max, Nx=nx, Ny=ny, 
                                matrix.list=matrix.list) 
    
    t.dif <- proc.time()[3] - t.dif
    t.advX <- proc.time()[3]
    
    dD <- temp-sum(species.cube)
    dD <- ifelse(round(dD, 14) != 0, paste('Diff.Warning', dD, '\n'), 'OK')
    temp <- sum(species.cube)
    
    # Advection in Space
    # _____________________   
    species.cube <- adv.upw.2d.3d(cube=species.cube, nx=nx, ny=ny, adv.p.x=adv.cubes[[1]], 
                                  adv.m.x=adv.cubes[[2]], adv.p.y=adv.cubes[[3]], 
                                  adv.m.y=adv.cubes[[4]], dx=dx, dy=dy, 
                                  it.cfl=it.adv.XY, Emax=e.max, dt.cfl=dt.adv.XY)
    
    t.advX <- proc.time()[3] - t.advX
    t.advE <- proc.time()[3]
    
    dA <- temp-sum(species.cube)
    dA <- ifelse(round(dA, 15) != 0, paste('AdvX.Warning', dA, '\n'), 'OK')
    temp <- sum(species.cube)
    
    # Advection in Energy
    # _____________________
    
    # update flows to account for laying cost dependency on field occupation
    e.ls.m.temp <- e.ls.m + rep(blur.laying * g.lay * (1 - species[,,1] / K), e.max)
    e.ls.p.temp <- e.ls.p + c(rep(blur.laying * g.lay * (1 - species[,,1] / K), e.max-1), rep(0, nx*ny))
    # update cfl compliant time steps
    max.adv.E <- max(abs(cbind(e.ls.m + rep(blur.laying * g.lay, e.max), 
                               e.ls.p + c(rep(blur.laying * g.lay, e.max-1), rep(0, nx*ny)))))
    it.adv.E <- ceiling(max.adv.E*time.step/(1/e.max))
    dt.adv.E <- time.step / pmax(it.adv.E, 1)
    
    species.cube <- adv.upw.energy(cube=species.cube, nx=nx, ny=ny, adv.p.e=e.ls.p.temp, 
                                   adv.m.e=e.ls.m.temp, de=dE, it.cfl=it.adv.E, 
                                   dt.cfl=dt.adv.E, Emax=e.max)
    
    
    t.advE <- proc.time()[3] - t.advE
    dADVE <- temp-sum(species.cube)
    dADVE <- ifelse(dADVE < 0, 'ERROR Adv.E. Flows\n', 'OK')
    
    # flatenning back the cube
    species[,,2] <- apply(species.cube, FUN=sum, c(1,2))
    
    # predictions & timer concatenation
    RES <- cbind(RES, apply(species, FUN=mean, 3))
    timeTotal <- cbind(timeTotal, c(t.rk=t.rk, t.diff=t.dif, t.advE=t.advE, 
                                    t.advX=t.advX, t.total=proc.time()[3] - clockTotal))
    
    df[t+1, 2:3] <- c(sum(species[,,1]), sum(species[,,2]))
    # display current status
    if (t %% verb == 0) {
      cat('\nDerivatives checks - ', 'dRK4', dRK4, ' ; dDif', dD, ' ; dAdvX', dA, ' ; dAdvE', dADVE)
      cat('\nsum(Y) = ', sum(species[,,1]), '   min(Y) = ', min(species[,,1]), '   max(Y) = ', max(species[,,1]))
      cat('\nsum(A) = ', sum(species[,,2]), '   min(A) = ', min(species[,,2]), '   max(A) = ', max(species[,,2]))
      cat('\nsum(A+Y) = ', sum(species[,,1])+sum(species[,,2]))
      cat('\nComp. times - ', 'RK4', t.rk, 'diff', t.dif, 'advX', t.advX, 'advE', t.advE, 'it.advE', it.adv.E)
      cat('\ntime-step ', t,' for date ', time[t], ' computed in ', timeTotal[length(timeTotal)],'s\n')
      if (plot){
        par(mfrow=c(2,2), mar=c(1,1.75,1,1))
        plot(1:e.max, apply(species.cube, 3, sum), type='b')
        plot(raster(species[,,1]))
        plot(raster(apply(species.cube, FUN=sum, c(1,2))))
        # plot(raster(species.cube[,,30]))
        plot(x=df$time, y=df$Y, type='l', ylim=c(0,max(df[,-1],na.rm=TRUE)), col='red')
        lines(x=df$time, y=df$A, type='l')
        
      }
    }
    
    if (t %% save.freq == 0) {
      youngs.density <- abind(youngs.density, species[,,1], along=3)
      adults.density <- abind(adults.density, species[,,2], along=3)
    }
    
  }
  
  cat('\nSimulation complete in ', proc.time()[3] - clockGlobal,'s, representing ', length(time), 
      ' time steps')
  
  
  list(youngs.density=youngs.density, 
       adults.density=adults.density, 
       times=timeTotal,
       df=df, parameters=parameters,
       status=status,
       time=proc.time()[3] - clockGlobal)
  
}
