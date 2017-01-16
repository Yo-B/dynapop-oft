# dynapop-oft
A population dynamics model accounting for the individual process of resource foraging on heterogeneous agricultural landscapes.

This is a GNU R program with 2 scripts loading 4 source files fo functions.

1) The first script is gen_landscape.R which generates a landscape from specified landscape metrics with a genetic algorithm
      the generated landscape will be plotted in vectorial and raster format in the output/ folder
      the landscape object is landscape.res in which you can find 
          landscape.res$sldf (spatial lines data frame) of the linear resource presence/absence
          landscape.res$spdf (spatial polygons data frame) of the polygonal resource presence/absence
      the landscape is converted into a raster/matrix called r.map which is used as substrate in the population dynamics model
      
2) The sencond script is sim_population.R
