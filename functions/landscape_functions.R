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
#   functions used by script oneLS.R
#   includes low level functions for spatial object manipulation
#   and GA function for controled landscapes generation
# _____________________________________________________________________

require(spatstat)
require(maptools)
require(deldir)
require(rgeos)



#==================================================================
#---    use deldir to get SP from voronoi, from Melen Leclerc  ----
#==================================================================
voronoi = function(layer,xmin,xmax,ymin,ymax) {
  crds = layer@coords     #extraction de la matrice des coordonn?es de l'objet layer
  z = deldir(crds[,1], crds[,2],rw=c(xmin,xmax,ymin,ymax)) # triangulation en fonction des points propos?s
  w = tile.list(z)      #pour chaque point-centre extrait le polygone-pavage
  polys = vector(mode='list', length=length(w))
  for (i in seq(along=polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP = SpatialPolygons(polys)
  SP@bbox[,1]<-c(xmin,ymin)
  SP@bbox[,2]<-c(xmax,ymax)
  return (SP)
}




#=============================================================
#---        Get spatialines from spatial polygons         ----
#=============================================================
as.spl <- function(sppol){
  yscale <- c(0, max(sppol@bbox)*10, 0, max(sppol@bbox)*10)
  line.list <- c()
  line.mat <- rep(max(sppol@bbox)*20,4)
  id <- 0
  for (i in seq(sppol)){
    a <- sppol@polygons[[i]]@Polygons[[1]]@coords
    for (j in 1:(nrow(a)-1)){     
      
      if (sum(apply(matrix(round(line.mat, 12) %in% 
                           round(cbind(a[-nrow(a),], a[-1,])[j,] + yscale, 12), ncol=4), FUN=sum, 1) == 4) == 0){
        id <- id + 1
        line.mat <- rbind(line.mat, cbind(a[-nrow(a),], a[-1,])[j,] + yscale)
        line.list <- c(line.list, Lines(list(Line(matrix(cbind(a[-nrow(a),], a[-1,])[j,],
                                                         nrow=2, byrow=TRUE))), ID=id))
      }
    }
  }
  SpatialLines(line.list)
}


#=============================================================
#---  get mid points (ppp) of segments of a spatial lines ----
#=============================================================
get.mid <- function(spl, window){
  mid.coord <- t(sapply(seq(spl), 
                        FUN=function(x) spl@lines[[x]]@Lines[[1]]@coords[1,] + 
                          diff(spl@lines[[x]]@Lines[[1]]@coords)/2))
  ppp(mid.coord[,1], mid.coord[,2], window=window) #owin(c(xlim[1], xlim[2]), c(ylim[1], ylim[2]))
}


#=====================================================================
#---    Metropolis towards homogeneous distance between seeds     ----
#=====================================================================
f.gibbs.points <- function(n=200, iterations=1000, ylim=c(0,1), xlim=c(0,1), disp=TRUE){
  
  coord <- cbind(runif(n,xlim[1],xlim[2]),
                 runif(n,ylim[1],ylim[2]))
  
  f.eval <- function(coord){
    ppp.temp <- ppp(coord[,1], coord[,2], window=owin(c(xlim[1], xlim[2]), c(ylim[1], ylim[2])))
    m <- pairdist(ppp.temp, periodic=TRUE, method="C")
    sum(1/m[lower.tri(m)]^2)
  }
  
  old.eval <- f.eval(coord)
  eval.vec <- old.eval
  
  for (i in 1:iterations){
    
    old.coord <- coord
    coord[sample(1:n, 1),] <- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
    new.eval <- f.eval(coord)
    
    eval.vec <- c(eval.vec, new.eval)
    
    if (new.eval > old.eval){ # rejection
      coord <- old.coord
    } else { # acceptance
      old.eval <- new.eval
    }
    
  }
  
  spp.new <- SpatialPoints(coord)
  if (disp==TRUE){
    par(mfrow=c(1,3))
    plot(spp.new)
    plot(voronoi(spp.new, xlim[1], xlim[2], ylim[1], ylim[2]))
    plot(eval.vec, type='l')
    par(mfrow=c(1,1))
  }
  
  ppp(coord[,1], coord[,2], window=owin(c(xlim[1], xlim[2]), c(ylim[1], ylim[2])))
}


#=============================================================
#---        controlled generation the tesselation
#---        structural basis of the landscape    
#=============================================================
init.struct <- function(n, iterations, xlim, ylim, prop, edges.bbox=FALSE){
  
  points <- f.gibbs.points(n=n, iterations=iterations, ylim=ylim, xlim=xlim, disp=FALSE)
  
  window <- owin(c(xlim[1], xlim[2]), c(ylim[1], ylim[2]))
  # generating structures from the point pattern
  spp <- SpatialPoints(cbind(points$x, points$y))
  sppol <- voronoi(spp, window$xrange[1], window$xrange[2], 
                   window$yrange[1], window$yrange[2])
  spl <- as.spl(sppol)
  psp <- as.psp.SpatialLines(spl)
  psp$window <- window
  ppp.edges <- get.mid(spl, window)
  ppp.poly <- ppp(x=spp@coords[,1], y=spp@coords[,2],  window=window)
  
  t(sapply(1:length(spl), FUN=function(x) spl@lines[[x]]@Lines[[1]]@coords))
  
  # edges not fully on bbox
  tabu.edges <- !(ppp.edges$x %in% ppp.edges$window$xrange + 
                    ppp.edges$y %in% ppp.edges$window$yrange) == 1
  # edges w/o any point on bbox
  tabu.edges.hard <- apply(t(sapply(1:length(spl), FUN=function(x) spl@lines[[x]]@Lines[[1]]@coords)), 1, 
                           FUN=function(x) sum(x[1:2] %in% ppp.edges$window$xrange, 
                                               x[3:4] %in% ppp.edges$window$yrange)) == 0
  # pol w/o any point on the bbox
  tabu.pol <- unlist(lapply(sapply(1:length(sppol), FUN=function(x) 
    unlist(sppol@polygons[[x]]@Polygons[[1]]@coords)), 
    FUN=function(x) sum(x[,1] %in% ppp.edges$window$xrange) + 
      sum(sum(x[,2] %in% ppp.edges$window$yrange)))) == 0
  
  list(ppp.poly=ppp.poly, ppp.edges=ppp.edges, sppol=sppol, spl=spl, 
       tabu.e=tabu.edges, tabu.e.h=tabu.edges.hard, tabu.p=tabu.pol)
  
}

#=============================================================
#---    computation of metrics on one landscape
#=============================================================
compute.metrics <- function(x, population=population, poly.mat.dist, edges.mat.dist,
                            p.e.mat.dist, p.e.mat.inter, areas, lengths, total.area,
                            resol, ncells, perimeters){
  
  # compute ENN_AW (Euclidean Nearest Neighbour Distance, Area Weighted)
  # Leitao, Measuring Landscapes, 2006
  m1 <- as.matrix(poly.mat.dist[population[[1]][x,], population[[1]][x,]])
  diag(m1) <- Inf
  m1 <- apply(m1, 1, min) * areas[population[[1]][x,]] / sum(areas[population[[1]][x,]])
  if(length(m1)==0) m1 <- Inf
  # ENN Length weighted for the segments
  m2 <- as.matrix(edges.mat.dist[population[[2]][x,], population[[2]][x,]])
  diag(m2) <- Inf
  m2 <- apply(m2, 1, min) * lengths[population[[2]][x,]] / sum(lengths[population[[2]][x,]])
  if(length(m2)==0) m2 <- Inf
  # ENN area weighted for the segments to polys
  m3 <- as.matrix(p.e.mat.dist[population[[2]][x,], population[[1]][x,]])
  m3 <- apply(m3, 2, min) * areas[population[[1]][x,]] / sum(areas[population[[1]][x,]])
  if(length(m3)==0) m3 <- Inf
  # ENN Length and area weighted for the segments and polys RECIPROQUELY
  m3r <- as.matrix(p.e.mat.dist[population[[2]][x,], population[[1]][x,]])
  m3r <- c(apply(m3r, 1, min) * lengths[population[[2]][x,]] / sum(lengths[population[[2]][x,]]), 
           apply(m3r, 2, min) * areas[population[[1]][x,]] / sum(areas[population[[1]][x,]]))
  
  # Mean distance Length and area weighted for the segments and polys
  md <- as.matrix(p.e.mat.dist[population[[2]][x,], population[[1]][x,]])
  md <- apply(md, 2, mean) * areas[population[[1]][x,]] / sum(areas[population[[1]][x,]])
  if(length(md)==0) md <- Inf
  # Mean distance Length and area weighted for the segments and polys RECIPROQUELY
  mdr <- as.matrix(p.e.mat.dist[population[[2]][x,], population[[1]][x,]])
  mdr <- c(apply(mdr, 1, mean) * lengths[population[[2]][x,]] / sum(lengths[population[[2]][x,]]), 
           apply(mdr, 2, mean) * areas[population[[1]][x,]] / sum(areas[population[[1]][x,]]))
  mdr <- na.omit(mdr)
  if(length(mdr)==0) mdr <- Inf
  
  # compute interface length between both resources, normalised by edges length
  m4 <- sum(unlist(apply(as.matrix(p.e.mat.inter[population[[2]][x,], population[[1]][x,]]), 1, unique))) / 
    sum(lengths[population[[2]][x,]])
  # normalised by crops perimeters, then Area weighted
  pIL <- apply(t(t(p.e.mat.inter[population[[2]][x,], population[[1]][x,]]) / 
                   perimeters[population[[1]][x,]]), 2, sum)
  m5 <- sum(pIL * areas[population[[1]][x,]] / sum(areas[population[[1]][x,]]))
  
  # compute proportions
  pl <- sum(population[[1]][x,] * areas) / total.area
  pf <- sum(ceiling(population[[2]][x,] * lengths / resol))/ncells
  
  return(c('ENN_l'=sum(m1), 'ENN_f'=sum(m2), 'ENN_fl'=sum(m3), 'ENN_flR'=sum(m3r), 
           'IL_f'=m4, 'IL_l'=m5, 'LCP_l'=pl, 'LCP_f'=pf, 'sd_l'=sd(m1), 'sd_f'=sd(m2), 
           'sd_fl'=sd(m3), 'MD_fl'=sum(md), 'MD_flR'=sum(mdr)))
}

#=============================================================
#---    GENETIC ALGORITHM, GA
#=============================================================

#   evolution functions
# _______________________________________
init.pop <- function(N, props, dim, tabu.e, tabu.p){ # random generation of population of landscapes (marks)
  poly.p <- rep(0, dim[1]*N) 
  poly.p[1:round(props[1]*dim[1]*N+sum(!tabu.p)*props[1])] <- 1
  poly.p <- matrix(sample(poly.p), nrow=N, byrow=TRUE)
  poly.p[,!tabu.p] <- 0  
  
  poly.e <- rep(0, dim[2]*N) 
  poly.e[1:round(props[2]*dim[2]*N+sum(!tabu.e)*props[2])] <- 1
  poly.e <- matrix(sample(poly.e), nrow=N, byrow=TRUE)
  poly.e[,!tabu.e] <- 0
  
  list(poly.p=(poly.p==1), poly.e=(poly.e==1))
}

eval.func <- function(x, targets, scale){ # evaluation func : difference between current metrics and targeted ones
  
  if(anyNA(x)) return(-Inf)
  
  -sum(((x-targets)*scale)^2, na.rm=TRUE)
  
}

crossover.func <- function(parents, max.cross = 5, len){ # multiple point crossing over between two parents
  crosspoint <- sort(unique(sample(1:len, sample(1:max.cross,1))))
  tab <- rbind(c(1, crosspoint+1), c(crosspoint,len))
  child <- c()
  for (i in 1:(length(crosspoint)+1)){
    child <- c(child, unlist(parents[[i %% 2 + 1]][tab[1, i]:tab[2, i]]))
  }
  child
}


mutate.func <- function(child, tau=.01, len, tabu){ # mutation func., randomly switch some marks
  who <- (1:len)[as.logical((runif(len) < tau)*tabu)]
  child[who] <- !child[who]
  child
}


#
#   GA main function
# _______________________________________
ls.gen.GA <- function(n=100, # number of polygons
                      ls.metrics=c(30,30,30,.1), #ENN_AW poly, edges, poly-edges, Interface lengths
                      props=c(.1,.04), # surface proportion of each resource
                      sd.max=NULL, # max sd values ENN_AW sould not outreach (specif 1 number for the 3 ENN or a vector of 3 values)
                      scale=rep(1,6), # weigths of each metrics in convergeance, needed to scale metrics with very different units
                      xlim=c(0,1000),  ylim=c(0,1000), # space dimensions
                      iterations=100, # number of iterations in Gibbs process disposing poly. seeds
                      N=100, # population size
                      max.gen=300, # number of generations = stopping conditions
                      max.gen.fail=10, # number of gen. without improvement of best fit. = stopping condition (useless)
                      r=.5, # proportion of the population renewed at each generation
                      tau=.00025, # mutation rate for new borns
                      crs.ovr=2, # maximum crossover point number during recombination
                      tournament=FALSE, # is selection of reproducing individuals a tournament ? (FALSE=fitness biased sampling)
                      k.tourn=0, # size of the tournaments
                      p.tourn=.2, # geometricly decreasing probability for the biased sampling inside tournaments
                      p.trunc=.2, # geom. decreasing prob. for the biased sampling at truncation (ending of generation), 0=only the best
                      verbose=10, # how often (in gen.) should it display status ?
                      plot=TRUE){ # display graphic results ?
  
  
  bias.tourn <- pmax(p.tourn*(1-p.tourn)^(0:(k.tourn-1)), 1e-20)
  bias.pop <- p.trunc*(1-p.trunc)^(0:(N-1)) 
  # bias.pop must be > 1e-35 otherwise "sample()" does not account for it
  if (bias.pop[floor(round(N*(1-r),5))] < 1e-30) {
    bias.pop[floor(round(N*(1-r),5))] <- 1
  }
  
  # needed to approximate surface taken by edges from their lengths
  resol <- 10
  ncells <- prod(diff(cbind(xlim, ylim))/resol)
  
  window <- owin(c(xlim[1], xlim[2]), c(ylim[1], ylim[2]))
  structures <- init.struct(n=n, xlim=xlim, ylim=ylim, iterations=iterations, prop=props)
  
  areas <- sapply(1:length(structures[[3]]), FUN=function(x) structures[[3]]@polygons[[x]]@area)
  total.area <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1])
  lengths <- sapply(1:length(structures[[4]]), 
                    FUN=function(x) sqrt(sum(diff(structures[[4]]@lines[[x]]@Lines[[1]]@coords)^2)))
  total.length <- sum(lengths) - xlim[2] * 2 - ylim[2] * 2
  
  edges.nb <- length(structures[[4]])
  
  population <- init.pop(N=N, dim=c(structures[[1]]$n, structures[[2]]$n), 
                         props = c(.1,.2), tabu.e=structures[['tabu.e.h']], 
                         tabu.p=structures[['tabu.p']])
  
  
  poly.mat.dist <- gDistance(structures[['sppol']], structures[['sppol']], byid=TRUE)
  edges.mat.dist <- gDistance(structures[['spl']], structures[['spl']], byid=TRUE)
  p.e.mat.dist <- gDistance(structures[['sppol']], structures[['spl']], byid=TRUE)
  # compute interface lengths matrix
  pts.coords <- sapply(1:length(structures[['spl']]), FUN=function(x) t(structures[['spl']]@lines[[x]]@Lines[[1]]@coords))
  pts.coords1 <- SpatialPoints(t(pts.coords)[,1:2])
  pts.coords2 <- SpatialPoints(t(pts.coords)[,3:4])
  p.e.mat.inter <- lengths * ((gDistance(structures[['sppol']], pts.coords1, byid=TRUE) == 0) * 
                                (gDistance(structures[['sppol']], pts.coords2, byid=TRUE) == 0))
  getPeri <- function(polygons){
    pts <- nrow(polygons@Polygons[[1]]@coords)
    sum(sqrt(apply((polygons@Polygons[[1]]@coords[-1,] - polygons@Polygons[[1]]@coords[-pts,])^2, 1, sum)))
  }
  perimeters <- sapply(1:length(structures[['sppol']]), FUN=function(x) getPeri(structures[['sppol']]@polygons[[x]]))
  
  renew <- floor(round(N*(1-r),5))
  fit.ls <- matrix(c(-Inf, NA, -Inf, rep(-Inf, N*5)), nrow=1)
  fail <- 0
  gen <- 1
  best.fit.ever <- -Inf
  
  while(fail < max.gen.fail && gen < max.gen){
    # compute landscape metrics
    metrics <- sapply(1:N, FUN=compute.metrics, population=population,
                      poly.mat.dist=poly.mat.dist, edges.mat.dist=edges.mat.dist, 
                      p.e.mat.dist=p.e.mat.dist, p.e.mat.inter=p.e.mat.inter,
                      areas=areas, lengths=lengths, total.area=total.area,
                      resol=resol, ncells=ncells, perimeters=perimeters)
    # compute fitness from metrics
    #     fitness <- apply(metrics[c(1:3, 5:7),], 2, FUN=eval.func, targets=c(ls.metrics, props), scale=scale)
    fitness <- apply(metrics[c('ENN_l', 'ENN_f', 'ENN_fl', 'IL_f', 'LCP_l', 'LCP_f'),], 2, 
                     FUN=eval.func, targets=c(ls.metrics, props), scale=scale)
    # add penalty for high sd if sd.max is specified
    if (!is.null(sd.max)) {
      sd.cost <- apply(metrics[c('sd_l', 'sd_f', 'sd_fl'),], 2, FUN=function(x) sum(pmax(0, x-sd.max) * scale[1:3]))
      fitness[!is.na(sd.cost)] <- fitness[!is.na(sd.cost)] + sd.cost[!is.na(sd.cost)]
    }
    # scaling fitness between 0 & 1, with log scaling to avoid every one but the noobest to have fit.p=1 (which makes no selection)
    fit.p <- fitness 
    fit.p[fit.p==-Inf] <- min(fit.p[fit.p > -Inf])
    fit.p <- -log(-fit.p)
    fit.p <- (fit.p-min(fit.p)) / (max(fit.p)-min(fit.p))
    
    if (gen %% verbose == 0){
      cat('\ngen', gen, '/', max.gen ,' ; mean fit ', mean(fitness), ' ; best fit', 
          max(fitness), ' ; fail ', fail , '/', max.gen.fail)
    }
    fit.ls <- rbind(fit.ls, c(mean(fitness[!is.infinite(fitness)]), sd(fitness[!is.infinite(fitness)]), 
                              max(fitness[!is.infinite(fitness)]), fitness, 
                              metrics[c('ENN_fl'),], metrics[c('IL_f'),], metrics[c('LCP_l'),], metrics[c('LCP_f'),]))
    
    if (fit.ls[nrow(fit.ls),3] > best.fit.ever){ # conserve the best landscape independantly of the generation
      best.fit.ever <- fit.ls[nrow(fit.ls),3]
      who <- order(fitness, decreasing=TRUE)[1]
      best.map.ever <- list(population[[1]][who,], population[[2]][who,])
      best.metrics.ever <- metrics[,who]
    }
    
    # stopping conditions
    if (fit.ls[nrow(fit.ls),3] <= fit.ls[nrow(fit.ls)-1,3]){
      fail <- fail + 1
    }else{
      fail <- 0
    }
    
    new.borns <- list(c(), c())
    
    for (j in 1:(N-renew)){
      # selection
      if (tournament){ # or tournament selection
        happy.coople <- sample(which(order(fitness, decreasing=TRUE) %in% sample(1:N, k.tourn)), 2, prob=bias.tourn)
      }else{ # or fitness proportionate selection
        happy.coople <- sample(1:N, 2, prob=fit.p)
      }
      # cross-over
      child.poly <- crossover.func(parents=apply(population[[1]][happy.coople,], 1, FUN=function(x) list(x)), 
                                   len=n, max.cross=crs.ovr)
      child.edges <- crossover.func(parents=apply(population[[2]][happy.coople,], 1, FUN=function(x) list(x)), 
                                    len=edges.nb, max.cross=edges.nb*tau)
      # mutation
      child.poly <- mutate.func(child.poly, tau, len=n, tabu=structures[['tabu.p']])
      child.edges <- mutate.func(child.edges, tau, len=edges.nb, tabu=structures[['tabu.e.h']])
      
      new.borns[[1]] <- rbind(new.borns[[1]], child.poly)
      new.borns[[2]] <- rbind(new.borns[[2]], child.edges)
      
    }
    # ordering by fitness
    population[[1]] <- population[[1]][order(fitness, decreasing=TRUE),]
    population[[2]] <- population[[2]][order(fitness, decreasing=TRUE),]
    # truncation and merging of the new borns to population
    if(p.trunc == 0){
      survivors <- 1:floor(round(N*(1-r),5)) # strict fitness ordinated truncation
    }else{
      survivors <- sort(sample(1:N, floor(round(N*(1-r),5)), prob=bias.pop)) # fitness biased sampling truncation
    }
    population[[1]] <- rbind(population[[1]][survivors,], new.borns[[1]])
    population[[2]] <- rbind(population[[2]][survivors,], new.borns[[2]])
    
    gen <- gen + 1
  }
  
  metrics <- metrics[,order(fitness, decreasing = TRUE)]
  fitness <- fitness[order(fitness, decreasing = TRUE)]
  
  
  structures[[1]]$marks <- population[[1]][1,]
  structures[[2]]$marks <- population[[2]][1,]
  
  
  spp <- SpatialPoints(cbind(structures[[1]]$x, structures[[1]]$y))
  sppol <- voronoi(spp, window$xrange[1], window$xrange[2], 
                   window$yrange[1], window$yrange[2])
  spl <- as.spl(sppol)
  psp <- as.psp.SpatialLines(spl)
  spdf <- SpatialPolygonsDataFrame(sppol, data=data.frame(marks=population[[1]][1,]))
  sldf <- SpatialLinesDataFrame(spl, data=data.frame(marks=population[[2]][1,]))
  
  
  if (plot){
    par(mfrow=c(3,5), mar=c(0,0,3,0))
    for (x in seq(1,15,by=1)){ 
      title <- paste(as.character(c(round(fitness[x],1), round(metrics[1:2,x],1),
                                    round(metrics[c( 'ENN_fl'),x],1)),
                                  '\n', round(metrics[5:7,x],2)), collapse=' ')
      plot(spdf, col=population[[1]][x,], main=title)
      plot(psp, col=c('grey', 'red')[population[[2]][x,]+1], 
           lwd=population[[2]][x,]+2, add=TRUE) 
    }
    
    par(mfrow=c(2,3), mar=c(5, 4, 4, 2) + 0.1)
    
    pal <- colorRampPalette(c("red","yellow","springgreen","royalblue"))(N)
    new.ls <- t(apply(fit.ls[,4:(N+3)], 1, FUN=function(x) x[order(x, decreasing=TRUE)]))
    plot(x=1:nrow(new.ls[,]), log(-new.ls[,1]), # !is.infinite(new.ls[,1])
         ylim=c(min(log(-as.numeric(new.ls)[!is.infinite(as.numeric(new.ls))])), 
                max(log(-as.numeric(new.ls)[!is.infinite(as.numeric(new.ls))]))),
         xlab='generations', ylab='log(-fitness)', cex=.1, col=pal[1], alpha=.5)
    for (i in 2:round(N/1)) lines(x=1:nrow(new.ls), log(-new.ls[,i]), type='p', cex=.1, col=pal[i], alpha=.5)
    
    new.ls <- fit.ls[,-(1:3)]
    plot(x=1:nrow(new.ls), new.ls[,(N+1)], cex=.1, ylab=dist.metric, xlab='generations')
    for (i in 2:round(N/1)) lines(x=1:nrow(new.ls), new.ls[,N+i], type='p', cex=.1)
    
    plot(x=1:nrow(new.ls), new.ls[,(3*N+1)], cex=.1, ylab='LCP (E. red, P. black)', xlab='generations', 
         ylim=c(min(min(new.ls[-1,(3*N+1):(3*N+round(N/1))]), min(new.ls[-1,(4*N+1):(4*N+round(N/1))])),
                max(max(new.ls[-1,(3*N+1):(3*N+round(N/1))]), max(new.ls[-1,(4*N+1):(4*N+round(N/1))]))))
    for (i in 2:round(N/1)) lines(x=1:nrow(new.ls), new.ls[,3*N+i], type='p', cex=.1)
    for (i in 1:round(N/1)) lines(x=1:nrow(new.ls), new.ls[,4*N+i], type='p', cex=.1, col="red")
    
    plot(as.numeric(new.ls[,(N+1):(2*N)]), as.numeric(new.ls[,(2*N+1):(3*N)]),
         ylab='Interface length', xlab= 'ENN_fl', cex=.1)
    
    plot(x=2:nrow(fit.ls), log(-fit.ls[-1,1] + fit.ls[-1,2]/sqrt(N)), 
         type='l', ylim=c(min(c(log(-fit.ls[!is.infinite(fit.ls[,1]),1] - fit.ls[!is.infinite(fit.ls[,1]),2]/sqrt(N)), 
                                log(-fit.ls[!is.infinite(fit.ls[,1]),3]))),
                          max(log(-fit.ls[!is.infinite(fit.ls[,1]),1] + fit.ls[!is.infinite(fit.ls[,1]),2]/sqrt(N)))),
         xlab='generations', ylab='log(-fitness)')
    lines(x=2:nrow(fit.ls), log(-fit.ls[-1,1] - fit.ls[-1,2]/sqrt(N)))
    lines(x=2:nrow(fit.ls), log(-fit.ls[-1,1]), col='red')
    lines(x=2:nrow(fit.ls), log(-fit.ls[-1,3]), col='green')
    
    par(mar=c(0,0,2,0))
    title <- paste(as.character(c(round(best.fit.ever,4), round(best.metrics.ever[c('ENN_l', 'ENN_f', 'ENN_fl')],3), 
                                  '\n', round(best.metrics.ever[c('IL_f','LCP_l','LCP_f')],3))), collapse=' ')
    plot(spdf, col=best.map.ever[[1]], main=title)
    plot(psp, col=c('grey', 'red')[best.map.ever[[2]]+1], 
         lwd=best.map.ever[[2]]+2, add=TRUE) 
    
  }
  
  return(list(spdf=spdf, sldf=sldf, spl=spl, psp=psp, population=population, metrics=metrics, fit.ls=fit.ls,
              map.res=best.map.ever, fit.res=best.fit.ever, metrics.res=best.metrics.ever,
              poly.mat.dist=poly.mat.dist, edges.mat.dist=edges.mat.dist, p.e.mat.dist=p.e.mat.dist, 
              p.e.mat.inter=p.e.mat.inter, areas=areas, lengths=lengths, perimeters=perimeters, 
              stop.gen=gen, fitness=fitness))
  
}



