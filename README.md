# dynapop-oft

The compagnon paper to this repository can be found in [Ecological Modelling](https://www.sciencedirect.com/science/article/pii/S0304380016302885).

A population dynamics model accounting for the individual process of resource foraging on heterogeneous agricultural landscapes.

This is a GNU R program with 2 scripts loading 4 source files containing functions. Both scripts need to be run sequentially, __gen_landscape.R__ before __sim_population.R__.


1. The first script is __gen_landscape.R__ which generates a landscape from specified landscape metrics with a genetic algorithm  
      * the generated landscape will be plotted in vectorial and raster format in the _output/_ folder  
      * the landscape object is _landscape.res_ in which you can find  
          * _landscape.res$sldf_ (spatial lines data frame) of the linear resource presence/absence  
          * _landscape.res$spdf_ (spatial polygons data frame) of the polygonal resource presence/absence  
      * the landscape is converted into a raster/matrix called _r.map_ which is used as substrate in the population dynamics model  
      
      
2. The second script is __sim_population.R__ which simulates a population dynamics on the previously generated landscape raster/matrix  
     * the population demographic and energetic parameters are first fitted  
     * then the spatial simulation can be done  
     * the resulting object _simulation_ contains densities for young and adult stages for all the time and space discretisations  
     * during computation, status will be displayed at each time step, as well as 4 plots:  
          a. the summed densities (over _x,y_) for each energy discretisation  
          b. the densities over _x,y_ of the young stage  
          c. the densities over _x,y_ of the adult stage  
          d. record of the temporal dynamics for both stages
            

### Required R packages:

_ggplot2,
raster,
grid,
gridExtra,
limSolve,
GGally,
nloptr,
stats,
Matrix,
abind,
plyr,
compiler,
spatstat,
maptools,
deldir,
rgeos_

Running smoothly with R version 3.3.2 (2016-10-31) -- "_Sincere Pumpkin Patch_"
