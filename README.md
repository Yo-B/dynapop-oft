# dynapop-oft
A population dynamics model accounting for the individual process of resource foraging on heterogeneous agricultural landscapes.

This is a GNU R program with 2 scripts loading 4 source files containing functions.


1) The first script is gen_landscape.R which generates a landscape from specified landscape metrics with a genetic algorithm

      the generated landscape will be plotted in vectorial and raster format in the output/ folder
      
      the landscape object is landscape.res in which you can find
      
          landscape.res$sldf (spatial lines data frame) of the linear resource presence/absence
          
          landscape.res$spdf (spatial polygons data frame) of the polygonal resource presence/absence
          
      the landscape is converted into a raster/matrix called r.map which is used as substrate in the population dynamics model
      
      
2) The sencond script is sim_population.R which simulates a population dynamics on the previously generated landscape raster/matrix

      the population demographic and energetic parameters are first fitted
      
      then the spatial simulation can be done
      
      the resulting object "simulation" contains densities for young and adult stages for all the time and space discretisations
      
      while computing it can display status at each time step as well as 4 plots :
      
            i) the summed densities (over x,y) for each energy discretisation
            
            ii) the densities over x,y of the young stage
            
            iii) the densities over x,y of the adult stage
            
            iv) record of the temporal dynamics for both stages
            
