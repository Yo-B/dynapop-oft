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
#   Numerical scheme functions called in RAD_functions.R
#     for the resolution of the opulation dynamics model
# _____________________________________________________________________



library(stats)
library(raster)
library(Matrix)
library(abind)
library(plyr)
library(compiler)
library(spatstat)


#=============================================================
#--- Solving the reaction terme with Runge-Kutta 4 scheme ----
#=============================================================

RK4 <- function(t, species, laying, par.demo, Nx, Ny){
  
  dY <- par.demo[1] * laying * species[,,2] * # laying / birth rate (+) 
    (1 - species[,,1] / par.demo[4]) - # carrying capacity
    par.demo[2] * species[,,1] - # stage development (-)
    par.demo[3] * species[,,1] * t # death rate (-)
    
  dA <- par.demo[2] * species[,,1] - # stage development (+)
    par.demo[3] * species[,,2] * t # death rate (-)
  
  array(c(dY, dA), dim = c(Nx, Ny, 2))
}

#=============================================================
#--------  Generating the 4 matrices for ADI scheme  ---------
#=============================================================
gen.adi.matrix <- function(dt = dt, # time step
                           D=DifMap, # Diffusion coefficient mapping (matrix w/ same dim as pop. matrix)
                           BC='neumann'){ # 'neumann' 'dirichlet' boundary condition for the diffusive part
  
  
  NxA <- nrow(D) ; dx <- 1/NxA
  NyA <- ncol(D) ; dy <- 1/NyA
  
  I <- bandSparse(n=(NxA)*(NyA),k=0, 
                  diagonal=list(rep(1,(NxA)*(NyA))))
  
  # multiplying coefficient for the principal diagonal boundary-related cells
  if (BC=='dirichlet'){
    bc <- 2
  }else if (BC=='neumann'){
    bc <- 1
  } 
  
  #     1st matrix generation
  # ____________________________________________________
  
  a <- b <- c <- matrix(0, NyA, NxA)
  for (j in 1:(NyA)){ 
    for (i in 2:(NxA-1)){   
      a[j,i] <- -dt/dx**2*(0.5*(D[i+1,j]+D[i,j]) + 0.5*(D[i-1,j]+D[i,j]))
      
      b[j,i] <- dt/dx**2*(0.5*(D[i+1,j]+D[i,j]))
      
      c[j,i] <- dt/dx**2*(0.5*(D[i-1,j]+D[i,j]))
      
    }
    b[j,1] <- dt/dx**2*(0.5*(D[2,j]+D[1,j]))
    c[j,NxA] <- dt/dx**2*(0.5*(D[NxA-1,j]+D[NxA,j]))
    
    a[j,1] <- -dt/dx**2*(0.5*(D[2,j]+D[1,j])) * bc
    a[j,NxA] <- -dt/dx**2*(0.5*(D[NxA-1,j]+D[NxA,j])) * bc
    
  }
  
  a <- as.vector(t(a)) ; b <- as.vector(t(b)) ; c <- as.vector(t(c))
  Axx <- bandSparse(n=(NxA)*(NyA),k=c(0,1,-1), 
                    diagonal=list(a,b[-length(b)], c[-1]))
  
  
  #     2nd matrix generation 
  # ____________________________________________________
  
  a <- b <- c <- matrix(0, NxA, NyA)
  for (i in 1:(NxA)){ 
    for (j in 2:(NyA-1)){   
      a[i,j] <- -dt/dy**2*(0.5*(D[i,j+1]+D[i,j]) + 0.5*(D[i,j-1]+D[i,j]))
      
      b[i,j] <- dt/dy**2*(0.5*(D[i,j+1]+D[i,j]))
      
      c[i,j] <- dt/dy**2*(0.5*(D[i,j-1]+D[i,j]))
      
    }
    b[i,1] <- dt/dy**2*(0.5*(D[i,2]+D[i,1]))
    c[i,NyA] <- dt/dy**2*(0.5*(D[i,NyA-1]+D[i,NyA]))
    
    a[i,1] <- -dt/dy**2*(0.5*(D[i,2]+D[i,1])) * bc
    a[i,NyA] <- -dt/dy**2*(0.5*(D[i,NyA-1]+D[i,NyA])) * bc
  }
  
  a <- as.vector(a) ; b <- as.vector(b) ; c <- as.vector(c)
  Ayy <- bandSparse(n=(NxA)*(NyA),k=c(0, NxA, -(NxA)), 
                    diagonal=list(a, 
                                  b[-(((NxA)*(NyA-1)+1):((NxA)*(NyA)))], 
                                  c[-(1:(NxA))]))
  
  Axx1 <- (I-0.5*Axx) 
  Axx2 <- (I+0.5*Axx) 
  Ayy1 <- (I-0.5*Ayy)
  Ayy2 <- (I+0.5*Ayy)
  
  #     Advective matrices generation for fokker-planck
  # ____________________________________________________
  mu.p.x <- -(rbind(D[-1,], D[NxA,]) - D) / dx
  mu.m.x <- -(D - rbind(D[1,], D[-NxA,])) / dx
  mu.p.y <- -(cbind(D[,-1], D[,NyA]) - D) / dy
  mu.m.y <- -(D - cbind(D[,1], D[,-NyA])) / dy
  
  list(Axx1=Axx1, Axx2=Axx2, Ayy1=Ayy1, Ayy2=Ayy2, mu.p.x=mu.p.x,
       mu.m.x=mu.m.x, mu.p.y=mu.p.y, mu.m.y=mu.m.y)
}

#=============================================================
#----------  Diffusion by ADI scheme on a square    ----------
#=============================================================
diffusion2d <- function(square, Nx=Nx, Ny=Ny, Axx1=Axx1, Axx2=Axx2, Ayy1=Ayy1, Ayy2=Ayy2){
  u1 <- as.numeric(square)
  
  u2 <- solve(Axx1 , Ayy2 %*% u1)
  u3 <- solve(Ayy1 , Axx2 %*% u2)
  
  square <- matrix(u3, Nx, Ny)
  
  square
}

#=============================================================
#----------  Diffusion by ADI scheme on the cube    ----------
#=============================================================
diffusion3d <- function(cube, Nx=Nx, Ny=Ny, Emax, matrix.list){
  # linearise matrix and concatenate them by column (3D -> 2D)
  u1 <- matrix(cube, Nx*Ny , ncol=Emax) # account for padding
  
  u2 <- solve(matrix.list['Axx1'][[1]] , matrix.list['Ayy2'][[1]] %*% u1) # this works !!
  u3 <- solve(matrix.list['Ayy1'][[1]] , matrix.list['Axx2'][[1]] %*% u2)
  
  cube <- array(u3, dim=c(Nx, Ny, Emax))
  
  cube
}

diffusion3d.fp <- function(cube, Nx=Nx, Ny=Ny, Emax, matrix.list){
  
  for (i in 1:Emax){
    u1 <- matrix(cube[,,i], Nx*Ny , ncol=1)
    u2 <- solve(matrix.list[[i]]['Axx1'][[1]], matrix.list[[i]]['Ayy2'][[1]] %*% u1)
    u3 <- solve(matrix.list[[i]]['Ayy1'][[1]], matrix.list[[i]]['Axx2'][[1]] %*% u2)
    cube[,,i] <- array(u3, dim=c(Nx, Ny, 1))
  }
  
  cube
}


#=============================================================
#----------  Space upwind on cube by layer    ----------
#=============================================================

adv.upw.2d.3d <- function(cube, nx, ny, adv.p.x, adv.m.x, adv.p.y, adv.m.y, dx, dy, it.cfl, dt.cfl, Emax){
  
  for (i in 1:Emax){
    
    for (j in 1:it.cfl[i]){
      F.p.x <- pmin(adv.p.x[,,i], 0) * rbind(cube[-1,,i], 0) + pmax(adv.p.x[,,i], 0) * cube[,,i]
      F.m.x <- pmin(adv.m.x[,,i], 0) * cube[,,i] + pmax(adv.m.x[,,i], 0) * rbind(0, cube[-nx, ,i]) 
      F.p.y <- pmin(adv.p.y[,,i], 0) * cbind(cube[,-1,i], 0) + pmax(adv.p.y[,,i], 0) * cube[,,i]
      F.m.y <- pmin(adv.m.y[,,i], 0) * cube[,,i] + pmax(adv.m.y[,,i], 0) * cbind(0, cube[, -ny,i]) 
      
      cube[,,i] <- cube[,,i] - dt.cfl[i]/dx*(F.p.x - F.m.x) - dt.cfl[i]/dy*(F.p.y - F.m.y)
      
    }
    
  }
  
  cube
}



#=============================================================
#-------  Energy flows by upwind scheme on the cube    -------
#=============================================================

adv.upw.energy <- function(cube, nx, ny, adv.p.e, adv.m.e, de, it.cfl, dt.cfl, Emax){
  # finite volume resolution of energy advection
  for (i in 1:it.cfl){
    
    F.p.e <- pmin(adv.p.e, 0) * abind(cube[,,-1], matrix(0,nx,ny), along=3) + pmax(adv.p.e, 0) * cube
    F.m.e <- pmin(adv.m.e, 0) * cube + pmax(adv.m.e, 0) * abind(matrix(0,nx,ny), cube[,,-Emax], along=3) 
    
    cube <- cube - dt.cfl/de*(F.p.e - F.m.e)
  }
  cube
}

#=============================================================
#-----------     Generating Advection fields     -------------
#=============================================================

gen.adv.fields <- function(d = 200, map = resources.map, resol = 15,
                               dt = 0.01, dx = .01, dy = .01, Nx=Nx, Ny=Ny){ 
  
  mask.poly <- map==1 * 1.0
  mask.edges <- map==2 * 1.0
  
  nc <- nr <- max(dim(map))
  if (nc %% 2 == 0){ # need odd numbers
    nc <- nr <- nc + 1
  }
  m <- matrix(nc=nc, nr=nr)
  col <- rep(1:nc, nc); row <- rep(1:nr, each=nr)
  x <- col - ceiling(nc/2); y <- row - ceiling(nr/2)
  m[cbind(row, col)] <- sqrt((x*resol)^2+(y*resol)^2)
  
  proj <- CRS('+proj=longlat +datum=WGS84')
  
  # if only one value of d is provided, it applies on both resources
  m1 <- 1 - m / d[1]
  m1[m1 < 0] <- 0
  kern.poly <- as.matrix(focal(x=raster(mask.poly*1.0, crs=proj), w=m1, 
                               pad=TRUE, padValue=0))
  m2 <- 1 - m / d[length(d)]
  m2[m2 < 0] <- 0
  kern.edges <- as.matrix(focal(x=raster(mask.edges*1.0, crs=proj), w=m2, 
                                pad=TRUE, padValue=0))
  # negative potential functions, normalised between 0 and -1
  kern.poly <- - kern.poly / max(kern.poly) 
  kern.edges <- - kern.edges / max(kern.edges) 
  
  # compute vectors from kernels
  edges.p.x <- -(rbind(kern.edges[-1,], kern.edges[Nx,]) - kern.edges) / dx
  edges.m.x <- -(kern.edges - rbind(kern.edges[1,], kern.edges[-Nx,])) / dx
  edges.p.y <- -(cbind(kern.edges[,-1], kern.edges[,Ny]) - kern.edges) / dy
  edges.m.y <- -(kern.edges - cbind(kern.edges[,1], kern.edges[,-Ny])) / dy
  
  poly.p.x <- -(rbind(kern.poly[-1,], kern.poly[Nx,]) - kern.poly) / dx
  poly.m.x <- -(kern.poly - rbind(kern.poly[1,], kern.poly[-Nx,])) / dx
  poly.p.y <- -(cbind(kern.poly[,-1], kern.poly[,Ny]) - kern.poly) / dy
  poly.m.y <- -(kern.poly - cbind(kern.poly[,1], kern.poly[,-Ny])) / dy
  
  tmp <- max(c(sqrt(edges.p.x^2+edges.p.y^2), sqrt(edges.m.x^2+edges.m.y^2)))
  edges.p.x <- edges.p.x / tmp
  edges.p.y <- edges.p.y / tmp
  edges.m.x <- edges.m.x / tmp
  edges.m.y <- edges.m.y / tmp
  
  tmp <- max(c(sqrt(poly.p.x^2+poly.p.y^2), sqrt(poly.m.x^2+poly.m.y^2)))
  poly.p.x <- poly.p.x / tmp
  poly.p.y <- poly.p.y / tmp
  poly.m.x <- poly.m.x / tmp
  poly.m.y <- poly.m.y / tmp
  
  
  return(list(poly.p.x, poly.m.x, poly.p.y, poly.m.y,
              edges.p.x, edges.m.x, edges.p.y, edges.m.y,
              kern.poly, kern.edges))
}

#=============================================================
#-----------     Melting fields in the cube      -------------
#=============================================================

gen.adv.cubes <- function(fields, Emax, max.speed = 40, resol = 15,
                              dt = 0.01, dx = .01, dy = .01, norm=TRUE,
                              func.poly=function(x) 1-x,
                              func.edges=function(x) 1-x){
  # (maxspeed) : m / dt -> 1000m/j -> 40 m/h (Finch & Skinner, 1975)
  # for a 1h time step, with 1T = 1000h ..
  # therefore precision is not enhanced by decreasing the time step
  
  max.adv.x <- max.speed * (1/dt) * dx / resol
  max.adv.y <- max.speed * (1/dt) * dy / resol
  
  adv.poly.x.p <- fields[[1]]; adv.poly.x.m <- fields[[2]]
  adv.poly.y.p <- fields[[3]]; adv.poly.y.m <- fields[[4]]
  adv.edges.x.p <- fields[[5]]; adv.edges.x.m <- fields[[6]]
  adv.edges.y.p <- fields[[7]]; adv.edges.y.m <- fields[[8]]
  
  adv.cube.x.p <- adv.cube.x.m <- adv.cube.y.p <- adv.cube.y.m <- array(0, dim=c(dim(adv.poly.x.m), Emax))
  
  x <- seq(0, 1, length.out=Emax)
  for (i in 1:Emax){
    adv.cube.x.m[,,i] <- func.poly(x[Emax-i+1]) * adv.poly.x.m + func.edges(x[i]) * adv.edges.x.m
    adv.cube.x.p[,,i] <- func.poly(x[Emax-i+1]) * adv.poly.x.p + func.edges(x[i]) * adv.edges.x.p
    adv.cube.y.m[,,i] <- func.poly(x[Emax-i+1]) * adv.poly.y.m + func.edges(x[i]) * adv.edges.y.m
    adv.cube.y.p[,,i] <- func.poly(x[Emax-i+1]) * adv.poly.y.p + func.edges(x[i]) * adv.edges.y.p
  }
  
  if (norm){ # normalise by ùax speed, do not if fokker-planck
    mx.spd <- max(c(sqrt(adv.cube.x.m^2+adv.cube.y.m^2), sqrt(adv.cube.x.p^2+adv.cube.y.p^2)))
    adv.cube.x.m <- adv.cube.x.m / mx.spd * max.adv.x
    adv.cube.y.m <- adv.cube.y.m / mx.spd * max.adv.y
    adv.cube.x.p <- adv.cube.x.p / mx.spd * max.adv.x
    adv.cube.y.p <- adv.cube.y.p / mx.spd * max.adv.y
  }

  return(list(adv.cube.x.p, adv.cube.x.m, adv.cube.y.p, adv.cube.y.m))
}


#=============================================================
#---------------      Useful Functions       -----------------
#=============================================================

# pad matrices with 0 fringes
# ----------------------------
padding.map <- function(map, vpad, hpad){
  if (hpad > 0){
    map <- cbind(matrix(0, nrow=nrow(map), ncol=hpad),
                 map,
                 matrix(0, nrow=nrow(map), ncol=hpad))
  }
  if (vpad > 0){
    map <- rbind(matrix(0, ncol=ncol(map), nrow=vpad),
                 map,
                 matrix(0, ncol=ncol(map), nrow=vpad))
  }
  map
}

#=================================================================
#-------- Compile all functions called in the RAD loop -----------
#=================================================================

# might be efficient in some cases
RK4 <- cmpfun(RK4)
diffusion2d <- cmpfun(diffusion2d)
diffusion3d <- cmpfun(diffusion3d)
diffusion3d.fp <- cmpfun(diffusion3d.fp)
adv.upw.2d.3d <- cmpfun(adv.upw.2d.3d)
adv.upw.energy <- cmpfun(adv.upw.energy)

