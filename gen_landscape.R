#                           dynapop-oft 
#       A population dynamics model accounting for the individual process 
#               of resource foraging on heterogeneous agricultural landscapes.
#     Copyright (C) 2017  Yoann BOURHIS at bourhis.yoann@gmail.com
#       see the corresponding publication for details at
#  https://www.sciencedirect.com/science/article/pii/S0304380016302885?via%3Dihub
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
#   design a landscape with 2 kinds of resource (linear and polygonal)
#   according to precise aggregation and proportion metrics
#   with genetic algorithm
#   plot it with ggplot
# _____________________________________________________________________

xlim <- ylim <- c(0,1000) # size of the landscape in meters
n <- xlim[2] * ylim[2] / (1000 * 1000) * 100 # polygon density = 100 / km2
iterations <- 100 # regularity of the polygonal seeds pattern

source('functions/landscape_functions.R')

# 6 metrics can be controled
target.metrics <- c(ENN_l = 100, # intra aggregation for laying sites (polygonal)
                    ENN_f = 45, # intra aggregation for feeding sites (linear)
                    ENN_fl = 80, # inter aggreg between feeding and laying resources (ISOLATION METRIC)
                    IL = .1, # interface length between laying and feeding resources (CONTIGUITY METRIC)
                    LCP_l = .1, # proportion of laying resource
                    LCP_f = .04) # proportion of feeding resource
# some metrics can be relaxed like the first two here (intra ressource agreggation),
# the other need a scale factor to have equivalent impact on the GA convergence
scale <- c(0, # ENN_l
           0, # ENN_f
           1, # ENN_fl , used as the scaling unit here
           pmax(1, target.metrics[3] / pmax(target.metrics[4], .1)), # interface length IL
           pmax(1, target.metrics[3] / target.metrics[5]), # LCP_f (more constraint is applied on proportions)
           pmax(1, target.metrics[3] / target.metrics[6])) # LCP_l

# Genetic algorithm parameters, see function ls.gen.GA in functionsLSGA.R for detailed meanings of the parameters
# increase N, max.gen and max.gen.fail to increase the computational effort to reach the target metrics
# modifying the other parameters could also improve convergence, but in a less straitforward way
N <- 300
max.gen <- 500
max.gen.fail <- 25 
r <- .5
tau <- .001
crs.ovr <- 3
tournament <- TRUE
k.tourn <- pmax(3, round(N/5))
p.tourn <- .2
p.trunc <- .3
plot <- FALSE ##
verbose <- 1

# run the genetic algorithm
landscape.res <- ls.gen.GA(n=n, ls.metrics=target.metrics[1:4], props=target.metrics[5:6], 
                           sd.max=target.metrics[7], xlim=xlim, scale=scale, ylim=ylim, N=N, 
                           max.gen=max.gen, max.gen.fail=max.gen.fail, r=r, tau=tau, 
                           crs.ovr=crs.ovr, tournament=tournament, k.tourn=k.tourn, 
                           iterations=iterations, p.tourn=p.tourn, p.trunc=p.trunc, 
                           plot=plot, verbose=verbose)
# resulting metrics of the best metric-compliant landscape
landscape.res$metrics.res[c(1:3,6:8)]
# distance from target metrics (in percent)
(landscape.res$metrics.res[c(1:3,6:8)]-target.metrics) / target.metrics * 100
# in this example, only the four last metrics are controled

# plotting the best landscape, i.e. the more in adequacy with target metrics
# the object is composed of a spdf object (polygonal resource) and sldf object (linear resource)
# for ggplot it is converted here in dataframe
require(ggplot2)
df <- NULL
df <- rbind(df, c(landscape.res$metrics.res[c('ENN_pe', 'IL_p', 'LCP_p', 'LCP_e')], 
                  ENN=paste('ENN =', landscape.res$metrics.res[c('ENN_pe')]), 
                  IL=paste('IL =', landscape.res$metrics.res[c('IL_p')])))
data <- landscape.res$sldf@data@.Data[[1]] # here marques are land-cover types
data <- data.frame(marques = data, length=landscape.res$sldf@data@.Data[[1]], ID = 1:length(data))
edges <- spChFIDs(landscape.res$sldf, as.character(row.names(data)))
edges <- SpatialLinesDataFrame(edges, data = data, match.ID = "ID")
data <- landscape.res$spdf@data@.Data[[1]] # here marques are land-cover types
data <- data.frame(marques = data, ID = 1:length(data))
poly <- spChFIDs(landscape.res$spdf, as.character(row.names(data)))
poly <- SpatialPolygonsDataFrame(poly, data = data, match.ID = "ID") # reconstruct the spdf
poly <- gBuffer(poly, byid=TRUE, width=0)
poly.id <- fortify(poly, region="ID")
poly.marques <- fortify(poly, region="marques")
edges.id <- fortify(edges, region="ID")
edges.id.marques <- c()
edges.id.marques <- c(edges.id.marques, rep(edges$marques, each=2))


leg.txt=c('matrix', 'crops')
colorScale=c('white', "#2c7fb8")

p1 <- ggplot() + 
  geom_path(data = poly.id, aes(y=lat, x=long, group = group), alpha=.1) +
  geom_polygon(data=poly.marques, aes(x = long, y = lat, group = group, fill = id)) +
  geom_path(data=poly.id, aes(x = long, y = lat, group=group), color="lightgray") +
  coord_equal() + 
  theme_bw() + 
  scale_fill_manual('Land-covers', values = colorScale, labels = leg.txt) + # defined legend based on predifined groups
  theme(legend.title = element_text(size=14), legend.text = element_text(size = 13), 
        legend.key.size = unit(1, "cm")) + labs(x = 'Distance from Origin in East Direction (m)', 
                                                y = 'Distance from Origin in North Direction (m)') + # classical informations for plots
  scale_x_continuous(breaks = round(seq(min(poly.id$long), 
                                        max(poly.id$long), by = 250),1)) +
  scale_y_continuous(breaks = round(seq(min(poly.id$lat), 
                                        max(poly.id$lat), by = 250),1)) +
  geom_path(data=edges.id[edges.id.marques==TRUE,],  
            aes(x = long, y = lat, group=id, color=rep('1',length(id))), size=2) +
  scale_color_manual('', label=c('field\nbanks'), values=c("#7fcdbb")) +
  guides(colour = guide_legend(order = 2), fill = guide_legend(order = 1)) 
p1

# the use of this landscape in the population dynamics model requires its conversion in a raster object
# i.e. a matrix of 0, 1 and 2 coding respectively for "no resource", laying resource and feeding resource
require(raster)
resol <- 10 # 10 meter sided pixels
ncbb <- landscape.res$spdf@bbox
xymin <- ncbb[,1]
r_xy <- apply(ncbb, 1, diff)
ncells <- ceiling(as.integer(r_xy / resol))
rasterTemplate <- raster(nrows = ncells[[2]], ncols=ncells[[1]],
                         xmn = ncbb[,1][[1]]-.000001, # for the edges on windows' borders to be accounted for ...
                         xmx=ncbb[, 1][[1]] + ncells[[1]] * resol+.000001,
                         ymn = ncbb[,1][[2]]-.000001,
                         ymx=ncbb[, 1][[2]] + ncells[[2]] * resol+.000001)

mod.spdf <-  landscape.res$spdf
mod.spdf$marks <-  landscape.res$spdf$marks
r.map <- rasterize(mod.spdf,  rasterTemplate, "marks")
r.edges <- rasterize( landscape.res$sldf[ landscape.res$sldf@data@.Data[[1]]==1,], rasterTemplate, "marks")
r.init <- r.map * 0 + target.metrics['LCP_l'] # homogeneous initial conditions for population dynamics model
r.map@data@values[r.edges@data@values==1] <- 2
val <- getValues(r.map)
xy <- as.data.frame(xyFromCell(r.map,1:ncell(r.map)))
xy <- cbind(xy,val)
colnames(xy) <- c('Longitude', 'Latitude', 'Density')


p2  <-  ggplot() + theme_bw() +
  geom_raster(data=xy, aes(x=Longitude, y=Latitude, fill=factor(Density), alpha=factor(Density))) + coord_equal() + 
  scale_fill_manual('Land-covers', values = c(colorScale, "#7fcdbb"), labels = c(leg.txt,'field banks')) + 
  scale_alpha_manual('Land-covers', values = c(.5, 1, 1), labels = c(leg.txt,'field banks')) + 
  theme(legend.title = element_text(size=14), legend.text = element_text(size = 13), 
        legend.key.size = unit(1, "cm")) + labs(x = 'Distance from Origin in East Direction (m)', 
                                                y = 'Distance from Origin in North Direction (m)') + # classical informations for plots
  scale_x_continuous(breaks = round(seq(min(poly.id$long), 
                                        max(poly.id$long), by = 250),1)) +
  scale_y_continuous(breaks = round(seq(min(poly.id$lat), 
                                        max(poly.id$lat), by = 250),1))
p2

# make a combined plot
require(grid)
require(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plot <- grid.arrange(p1 + theme(legend.position = "none") + ggtitle('A') +
                       theme(plot.title=element_text(hjust=0.01, size=20, face="bold", margin = margin(b = -20))), 
                     p2 + theme(legend.position = "none")  + labs(y='') + ggtitle('B') +
                       theme(plot.title=element_text(hjust=0.01, size=20, face="bold", margin = margin(b = -20))), 
                     g_legend(p1), ncol=3, widths=c(1,1, .3))


# export plot as a pdf
pdf('outputs/oneLS.pdf', width=12, height=5)
grid.newpage() 
grid.draw(plot)
dev.off()

