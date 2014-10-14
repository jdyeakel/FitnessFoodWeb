rm(list=c(ls()))

library(RColorBrewer)
library(Rcpp)
library(igraph)
library(beanplot)
library(lattice)
source("src/filled_contour.r")
source("src/smooth_pal.R")



#Import the RData for a given SDP solution matrix
#Uniform Habitat
load("Hab0.Rdata")

#Patchy Habitat
load("Hab1.Rdata")



##########################
# Set Initial Conditions #
##########################


############################
# Begin Forward Iterations #
############################


