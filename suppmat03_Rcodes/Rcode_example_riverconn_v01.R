#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
#remove everything in the working environment, without a warning!!
rm(list=ls())

#https://www.sciencedirect.com/science/article/pii/S1364815222001748
#https://github.com/damianobaldan/riverconn.git

#https://damianobaldan.github.io/riverconn_tutorial/
# install package if required
if(!require(ggnetwork)){
  install.packages("ggnetwork")
}
if(!require(RANN)){
  install.packages("RANN")
}
if(!require(elevatr)){
  install.packages("elevatr")
}
if(!require(corrmorant)){
  library("devtools")
  devtools::install_github('r-link/corrmorant', quiet =F)
  #install.packages("corrmorant", repos='http://cran.us.r-project.org')
}
library(corrmorant)
if(!require("riverconn")){
  if (!require('devtools')) install.packages('devtools')
  library("devtools")
  # install from GitHub
  devtools::install_github('damianobaldan/riverconn', build_vignettes = T)
}
# get libraries
library(tidyverse)
library(sf)
library(raster)
library(ggspatial)
library(viridis)
library(igraph)
library(riverconn)
library(elevatr)
library(gridExtra)
library(ggnetwork)
library(lwgeom)
library(gridExtra)
library(corrmorant)
library(RANN)
library(ggpubr)
library(cowplot)


# river network shapefile can be retrieved from HydroSHEDS, a global hydrographic atlas
# https://www.hydrosheds.org/products

shape_river <- st_read("Ebro_shape_corrected/Ebro_rivers.shp")
shape_basin <- st_read("Ebro_catchment/Ebro_mask.shp")
shape_dams <- st_read("Ebro_dams/Embalses/Embalses.shp")
