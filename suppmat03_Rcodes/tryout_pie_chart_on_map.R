#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#remove everything in the working environment, without a warning!!
rm(list=ls())

#define path for input dir
wd01 <- "/suppmat01_inp_files"
wd05 <- "/suppmat05_out_files"
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
wd_out01 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis"
wd00_wd01 <- paste0(wd00,wd01)
wd00_wd05 <- paste0(wd00,wd05)
#df_hap_loc04
## install all packages required
if(!require(mapdata)){
  install.packages("mapdata")
  library(mapdata)
}

if(!require(sp)){
  install.packages("sp")
  library(sp)
}

if(!require(ggmap)){
  install.packages("ggmap")
  library(ggmap)
}

if(!require(mapproj)){
  install.packages("mapproj")
  library(mapproj)
}

if(!require(rworldmap)){
  install.packages("rworldmap")
  library(rworldmap)
}

if(!require(akima)){
  install.packages("akima")
  library(akima)
}

## import all required libraries from packages
library(sp)
library(ggmap)
library(mapproj)
library(rworldmap)
library(maps)
library(mapdata)
library(akima)

## install the package 'scales', which will allow you to make points on your plot more transparent
if(!require(scales)){
  install.packages("scales")
  library(scales)
}
library("scales")


if(!require(fields)){
  install.packages("fields")
  library(fields)
}
library(fields)

## install the package 'marmap', which will allow you to plot bathymetric maps
if(!require(marmap)){
  install.packages("marmap")
  library(marmap)
}
library(marmap)

#get the package that enables the function 'subplot'
if(!require(TeachingDemos)){
  install.packages("TeachingDemos")
  library(TeachingDemos)
}
library(TeachingDemos)
#_______________make map for graphical abstract with average eDNA levels for all 6 species as histograms on trawl positions ______________________
# figname01 <- "Fig04_v07_map_haplotype_pie.pdf"
# pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# # set to save plot as pdf file with dimensions 8.26 to 2.9
# # 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# # and 210 mm and 74.25 mm matches 1/4 of a A4 page
# pdf(pthfignm01,width=(1.6*8.2677),height=(1.6*2.9232))
library(marmap)
# get bathymetric map
bathymap <- marmap::getNOAA.bathy(lon1 = 4, lon2 = 17,
                        lat1 = 52.88, lat2 = 60, resolution = 2)
# paste path and filename for input file together
loc04fl <- paste(wd00_wd05,"/loc04.csv",sep="")
#read in csv
df_hap_loc04 <- read.table(loc04fl, sep=",")
#change column names
colnames(df_hap_loc04) <- df_hap_loc04[1,]
# drop first row with column names
df_hap_loc04 <- df_hap_loc04[-1,]
#copy data frame
df_hl05 <- df_hap_loc04
# Creating a custom palette of blues for the bathymetric background
blues <- c("midnightblue", "lightsteelblue3","lightsteelblue2", "white")
# define a colour for transparent points
transp_col <- rgb(0, 0, 0, 0)
# get number of columns
nclhl05 <- ncol(df_hl05)
# get column 2 to the fifth last column
df_hl06 <- df_hl05[,c(2:(nclhl05-4))]
#make it a 'matrix array'
ma <- as.matrix(df_hl06)
if (!exists("zpfl"))
{
# Download the shapefile.
download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip" 
              , destfile=paste0(wd_out01,"/world_shape_file.zip"))
  } #else {flnmshpfl <- paste0(wd00_wd05,"/TM_WORLD_BORDERS_SIMPL-0.3.shp")}
# Make a full path to this file
zpfl <- paste0(wd_out01,"/world_shape_file.zip")
# Unzip this file. You can do it with R (as below),
unzip(zpfl, exdir = wd_out01)
# get a list of the files inside the zip file
lzpfl <- unzip(zpfl,list=T)
#grep in the names of files inside the list of files
flnmshpfl <- lzpfl$Name[grepl("shp",lzpfl$Name)]
#  -- > You now have 4 files. One of these files is a .shp file! (TM_WORLD_BORDERS_SIMPL-0.3.shp)
# Read this shape file with the rgdal library. 
library(rgdal)
mspdf <- rgdal::readOGR( 
  dsn= paste0(wd_out01,"/",flnmshpfl) , 
  layer="TM_WORLD_BORDERS_SIMPL-0.3",
  verbose=FALSE
)
#
library(shapefiles)
sfl<- gsub(".shp","",paste0(wd_out01,"/",flnmshpfl))
wlshpfl <- shapefiles::read.shapefile(sfl)
# change the row names
rownames(ma) <- df_hl05$smplloca
# get package library
library(mapplots)
# define plot limits
xlim = c(4, 16)
ylim = c(53, 59)
## remove column
df_hl05$rws2 <- NULL
# re arrange from wide to long
df_hl06 <- reshape2::melt(df_hl05,                                 # Apply melt function
                   id.vars = c("smplloca", "dec_lat","dec_lon","rws"))
# copy by renaming columns
df_hl06$Hpt <- df_hl06$variable
df_hl06$Cnt <- df_hl06$value
#remove columns
df_hl06$variable <-  NULL
df_hl06$value <-  NULL
df_hl06$Cnt <- as.numeric(df_hl06$Cnt)
# make xyz list 
xyz.Hpt <- make.xyz(df_hl06$dec_lon,df_hl06$dec_lat,df_hl06$Cnt,df_hl06$Hpt)
# count number of haplotypes
nHpt <- length(unique(df_hl06$Hpt))
# make a colour range to reflect the number of haplotypes
colra <- rainbow(nHpt)
# try making a different colour range
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
# make a colour ramp
colfunc <- colorRampPalette(cbbPalette2)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# use the colour function
cl03 <- colfunc(nHpt)
# Replace the colour range with rainbow colours
colra <- cl03
# read in map
world3 <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
# make it a "SpatialPolygonsDataFrame"
w3sp<- as(world3, 'Spatial')
#shapefiles::convert.to.shapefile(w3sp, field = )
# end plot
dev.off()


library(shapefiles)
#______________________________________________________________
#___________start  plot map halpotype v 07____________________
#______________________________________________________________
# define name for output plot file
figname01 <- "Fig04_v07_map_haplotype_pie.pdf"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
pdf(pthfignm01,width=(1.6*2.9),height=(0.4*8.26))

#add extra space to the right of the plot
par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 3), xpd=FALSE)
#reset this parameter
par(mfrow = c(1, 1)) 

#shapefiles::convert.to.shapefile(w3sp, field = )
mapplots::basemap(xlim, ylim, main = NA, bg="white", las=1,
                  
                  #surpress tickmarks
                  xaxt="n", yaxt="n")
mtext("a", side=3, adj=0, line=0.4, cex=1.4, font=2)
#add the coastline
draw.shape(wlshpfl, col="azure3")
#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 0.6, col=scales::alpha(colra,0.6))
# get unique haplotypes
Hpts <- unique(df_hl06$Hpt)
# bind colors and unique haplotypes in a data frame
df_colra <- as.data.frame(cbind(as.character(unique(df_hl06$Hpt)),colra))
# assign new column names to data frame
colnames(df_colra) <- c("HptNo","hexcol")
#exclude haplotypes that have a count of 0
df_hl07  <- df_hl06[!df_hl06$Cnt==0,]
# match haplotype to get hex color back to data frame
df_hl07$hexcol <- df_colra$hexcol[match(df_hl07$Hpt,df_colra$HptNo)]
# add a legend to the plot
# legend.pie(5,58,labels=c(df_hl07$Hpt), radius=0.9, bty="n", col=df_hl07$hexcol,
#            cex=0.8, label.dist=1.3)
# add a legend to the plot
# legend.pie(5,58,labels=c(Hpts), radius=0.9, bty="n", col=colra,
#            cex=0.8, label.dist=1.3)
# define a vector that holds labels for tick marks
# and use this on the axis instead
longwE <- paste(5*(0:4),"ºE",sep="")
axis(1, at=c(5*(0:4)), labels=longwE,cex.axis = 0.8)
lattwN <- paste((53:59),"ºN",sep="")
axis(2, at=c(53:59), labels=lattwN, las=1,cex.axis = 0.8)
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(3,54,z=legend.z,round=0,maxradius=(1*0.6),bty="n",txt.cex=0.6)
text(3,55.5,"individuals",cex=0.8) 
#add legend to plot
legend("topright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.4, pt.cex = 1.2, ncol=3,
       x.intersp = 0.6,y.intersp = 0.7)
# end plot
dev.off()
#______________________________________________________________
#___________end plot map halpotype v 07____________________
#______________________________________________________________


#______________________________________________________________
#___________start plot map halpotype v 07____________________
#______________________________________________________________
# define name for outpu plot file
figname01 <- "Fig04_v08_map_haplotype_pie.pdf"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
pdf(pthfignm01,width=(1.6*2.9),height=(0.4*8.26))
#add extra space to the right of the plot
par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 3), xpd=FALSE)
#reset this parameter
par(mfrow = c(1, 1)) 
# begin plot, with defined borders
plot(NA,NA, xlim = c(4, 16), ylim = c(53, 59),
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     xaxt="n", yaxt="n",
     xlab="Longitude", ylab="Latitude")
     #xlab = NA, ylab = NA)
# #plot bathymetry
# plot(bathymap,  add=F, lwd = c(0, 0), col=transp_col, image = TRUE, bpal = c(alpha(blues, 0.4)),
#      xlim = c(4, 16), ylim = c(53, 59),
#      cex=2.8, vfont = c("sans serif", "plain")
# )
# #plot coastline
# plot(bathymap, add=TRUE, lwd = c(0.8, 1.6), lty = c(1, 1),
#      xlim = c(4, 16), ylim = c(53, 59),
#      deep = c(-500, 0), shallow = c(-20, 0), step = c(20, 0),
#      cex=2.8, vfont = c("sans serif", "plain"),
#      col = c("#00009B", "black"), drawlabels = c(TRUE, TRUE))
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
    xlim = c(4, 16), ylim = c(53, 59),
    xlab = "Longitude", ylab = "Latitude",
    col="azure3", bg=transp_col)
#par(new=T)
mtext("a", side=3, adj=0, line=0.4, cex=1.2, font=2)
#dev.off()
mapplots::draw.pie(df_hl05$dec_lon,
                   df_hl05$dec_lon,
                   ma,
                   radius = df_hl05$rws, add = TRUE,
                   col=c(colra)
                   )
longwE <- paste(5*(0:4),"ºE",sep="")
axis(1, at=c(5*(0:4)), labels=longwE,cex.axis = 0.8)
lattwN <- paste((53:59),"ºN",sep="")
axis(2, at=c(53:59), labels=lattwN, las=1,cex.axis = 0.8)

#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 0.6, col=scales::alpha(colra,0.7))
#class(mspdf)
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(5,54,z=legend.z,round=0,maxradius=(1*0.6),bty="n",txt.cex=0.6)
text(5,55.5,"individuals",cex=0.8) 
#add legend to plot
legend("topright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.4, pt.cex = 1.2, ncol=3,
       x.intersp = 0.6,y.intersp = 0.7)
# end plot
dev.off()
#______________________________________________________________
#___________end plot map halpotype v 08____________________
#______________________________________________________________

#https://stackoverflow.com/questions/21676721/r-plot-circular-histograms-rose-diagrams-on-map

if(!require(shapefiles)){
  install.packages("shapefiles")
  library(shapefiles)
}

if(!require(mapplots)){
  install.packages("mapplots")
  library(mapplots)
}

library(mapplots)
library(shapefiles)

