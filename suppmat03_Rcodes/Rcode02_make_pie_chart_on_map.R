#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#remove everything in the working environment, without a warning!!
#rm(list=ls())
#define path for input dir
# wd01 <- "/suppmat01_inp_files"
# wd05 <- "/suppmat05_out_files"
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
wd_out01 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis"
wd00_wd01 <- paste0(wd00,wd01)
wd00_wd05 <- paste0(wd00,wd05)
#https://stackoverflow.com/questions/21676721/r-plot-circular-histograms-rose-diagrams-on-map
if(!require(shapefiles)){
  install.packages("shapefiles")
  library(shapefiles)
}
if(!require(mapplots)){
  install.packages("mapplots")
  library(mapplots)
}
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
## install the package 'scales', which will allow you to make points on your plot more transparent
if(!require(scales)){
  install.packages("scales")
  library(scales)
}
## install  package
library("scales")
## install  package
if(!require(fields)){
  install.packages("fields")
  library(fields)
}
# ## install the package 'marmap', which will allow you to plot bathymetric maps
# if(!require(marmap)){
#   #install.packages("marmap", repos='http://cran.us.r-project.org')
#   library(marmap)
# }
#get the package that enables the function 'subplot'
if(!require(TeachingDemos)){
  install.packages("TeachingDemos")
  library(TeachingDemos)
}
## import all required libraries from packages
library(sp)
library(ggmap)
library(mapproj)
library(rworldmap)
library(maps)
library(rgdal)
library(mapdata)
library(akima)
library(fields)
library(TeachingDemos)
library(marmap)
library(mapplots)
library(shapefiles)

#_______________make map for graphical abstract with average eDNA levels for all 6 species as histograms on trawl positions ______________________
# figname01 <- "Fig04_v07_map_haplotype_pie.pdf"
# pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# # set to save plot as pdf file with dimensions 8.26 to 2.9
# # 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# # and 210 mm and 74.25 mm matches 1/4 of a A4 page
# pdf(pthfignm01,width=(1.6*8.2677),height=(1.6*2.9232))

# get bathymetric map
# bathymap <- marmap::getNOAA.bathy(lon1 = 4, lon2 = 17,
#                         lat1 = 52.88, lat2 = 60, resolution = 2)
# paste path and filename for input file together
#loc04fl <- paste(wd00_wd05,"/loc04.csv",sep="")
#read in csv
#df_hap_loc04 <- read.table(loc04fl, sep=",")
#df_hap_loc06
#change column names
# colnames(df_hap_loc04) <- df_hap_loc04[1,]
# # drop first row with column names
# df_hap_loc04 <- df_hap_loc04[-1,]
#copy data frame
df_hl05 <- df_hap_loc08
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
#______________________________________________________________
#___________start  plot map halpotype v 07____________________
#______________________________________________________________
# define name for output plot file
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
# define a vector that holds labels for tick marks
# and use this on the axis instead
#https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
longwE <- paste(5*(0:4),"\u00B0E",sep="")
axis(1, at=c(5*(0:4)), labels=longwE,cex.axis = 0.6)
lattwN <- paste((53:59),"\u00B0N",sep="")
axis(2, at=c(53:59), labels=lattwN, las=1,cex.axis = 0.6)
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(3,54,z=legend.z,round=0,maxradius=(1*0.6),bty="n",txt.cex=0.6)
text(3,55.5,"samples",cex=0.8) 
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
#___________start plot map halpotype v 08____________________
#______________________________________________________________
# define name for outpu plot file
figname01 <- "Fig04_v09_map_haplotype_pie.pdf"
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
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
    xlim = c(4, 16), ylim = c(53, 59),
    xlab = "Longitude", ylab = "Latitude",
    col="azure3", bg=transp_col)
#add new title
mtext("a", side=3, adj=0, line=0.4, cex=1.2, font=2)
#dev.off()
mapplots::draw.pie(df_hl05$dec_lon,
                   df_hl05$dec_lon,
                   ma,
                   radius = df_hl05$rws, add = TRUE,
                   col=c(colra)
                   )
#https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
longwE <- paste(5*(0:4),"\u00B0E",sep="")
axis(1, at=c(5*(0:4)), labels=longwE,cex.axis = 0.6)
lattwN <- paste((53:59),"\u00B0N",sep="")
axis(2, at=c(53:59), labels=lattwN, las=1,cex.axis = 0.6)
#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 0.6, col=scales::alpha(colra,0.7))
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(5,54,z=legend.z,round=0,maxradius=(1*0.6),bty="n",txt.cex=0.6)
text(5,55.5,"samples",cex=0.8) 
#add legend to plot
legend("topright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.3, pt.cex = 0.8, ncol=3,
       x.intersp = 0.6,y.intersp = 0.7)
# end plot
dev.off()
#______________________________________________________________
#___________end plot map halpotype v 08____________________
#______________________________________________________________

#_______________________________________________________________________________
# start - make plot for haplotypes in Europe
#_______________________________________________________________________________
# paste path and filename for input file together
loc06fl <- paste(wd00_wd05,"/loc06.csv",sep="")
#read in csv
df_hap_loc06 <- read.table(loc06fl, sep=",")
#change column names
colnames(df_hap_loc06) <- df_hap_loc06[1,]
# drop first row with column names
df_hap_loc06 <- df_hap_loc06[-1,]

#copy data frame
df_hl07 <- df_hap_loc06
# get number of columns
nclhl07 <- ncol(df_hl07)
colnames(df_hl07)[grepl("ov.aL",colnames(df_hl07))] <- "smplloca"

# get column 2 to the 7th last column
df_hl08 <- df_hl07[,c(2:(nclhl07-6))]
rownames(df_hap_loc06)
#df_hap_loc06$dec_loc2

# make all columns numeric
df_hl08 <- as.data.frame(lapply(df_hl08,as.numeric))
#make it a 'matrix array'
ma <- as.matrix(df_hl08)


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
sfl<- gsub(".shp","",paste0(wd_out01,"/",flnmshpfl))
wlshpfl <- shapefiles::read.shapefile(sfl)
#
#df_hl07$smplloca <- df_hl07$dec_loc2
# change the row names
rownames(ma) <- df_hl07$smplloca
# define plot limits
xlim = c(-12, 16)
ylim = c(34.8, 58.0)
## remove column
df_hl07$rws2 <- NULL
df_hl07$dec_loc2 <- NULL
df_hl07$dec_loc3 <- NULL
# re arrange from wide to long
df_hl08 <- reshape2::melt(df_hl07,                                 # Apply melt function
                          id.vars = c("smplloca", "dec_lat","dec_lon","rws"))
#
df_hl08$dec_lat <- as.numeric(df_hl08$dec_lat)
df_hl08$dec_lon <- as.numeric(df_hl08$dec_lon)
df_hl08$rws <- as.numeric(df_hl08$rws)
# copy by renaming columns
df_hl08$Hpt <- df_hl08$variable
df_hl08$Cnt <- df_hl08$value
#remove columns
df_hl08$variable <-  NULL
df_hl08$value <-  NULL
df_hl08$Cnt <- as.numeric(df_hl08$Cnt)
# make xyz list 
xyz.Hpt <- mapplots::make.xyz(df_hl08$dec_lon,df_hl08$dec_lat,df_hl08$Cnt,df_hl08$Hpt)
# count number of haplotypes
nHpt <- length(unique(df_hl08$Hpt))
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
# get library for package
library(shapefiles)

#______________________________________________________________
#___________start plot map halpotype v 09____________________
#______________________________________________________________
# define name for outpu plot file
figname01 <- "Fig04_v10_map_haplotype_pie.pdf"
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
plot(NA,NA, xlim = c(-12, 22),
     ylim = c(37, 58.0),
     # use las to turn the tick labels to be horizontal
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     xaxt="n", yaxt="n",
     xlab="Longitude", ylab="Latitude")
#xlab = NA, ylab = NA)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          xlim = c(-12, 36),
          ylim = c(28, 58.0),
          xlab = "Longitude", ylab = "Latitude",
          col="azure3", bg=transp_col)
#par(new=T)
mtext("b", side=3, adj=0, line=0.4, cex=1.2, font=2)
#dev.off()

mapplots::draw.pie(df_hl08$dec_lat,
                   df_hl08$dec_lon,
                   ma,
                   radius = df_hl08$rws, add = TRUE,
                   col=c(colra)
)
# https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
longwE <- paste(seq(-12, 16,4),"\u00B0E",sep="")
longwW <- paste(c(12,8,4),"\u00B0W",sep="")
longw0 <- paste("0","\u00B0",sep="")
longwE <- paste(c(4,8,12),"\u00B0E",sep="")
longwW0E<- c(longwW,longw0,longwE)
axis(1, at=c(-12,-8,-4,0,4,8,12), labels=longwW0E,cex.axis = 0.6)
lattwN <- paste(round(seq(34.8, 58.0,4),0),"\u00B0N",sep="")
axis(2, at=c(round(seq(34.8, 58.0,4),0)), labels=lattwN, las=1,cex.axis = 0.6)

#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 1.6, col=scales::alpha(colra,0.7))
#class(mspdf)
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(-10,48,z=legend.z,round=0,maxradius=(1.6),bty="n",txt.cex=0.6)
text(-10,46.5,"samples",cex=0.8) 
#add legend to plot
legend("bottomright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.3, pt.cex = 1.2, ncol=3,
       x.intersp = 0.6,y.intersp = 0.7)
# end plot
dev.off()


#______________________________________________________________
#___________end plot map halpotype v 09____________________
#______________________________________________________________



#______________________________________________________________
#___________start plot map halpotype v 10____________________
#______________________________________________________________
# define name for outpu plot file
figname01 <- "Fig04_v11_map_haplotype_pie.pdf"
figname01 <- "Fig04_v11_map_haplotype_pie.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 3600, height = 4800, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 2), xpd=FALSE)
par(oma=c(0, 0, 0, 0))
#reset this parameter
par(mfrow = c(2, 1)) 


# begin plot, with defined borders
plot(NA,NA, xlim = c(4, 16), ylim = c(53, 59),
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     xaxt="n", yaxt="n",
     xlab="Longitude", ylab="Latitude",
     cex.lab = 0.7)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          xlim = c(4, 16), ylim = c(53, 59),
          xlab = "Longitude", ylab = "Latitude",
          col="azure3", bg=transp_col)
#add new title
mtext("a", side=3, adj=0, line=0.4, cex=1.2, font=2)
#dev.off()
mapplots::draw.pie(df_hl05$dec_lon,
                   df_hl05$dec_lon,
                   ma,
                   radius = df_hl05$rws, add = TRUE,
                   col=c(colra)
)
#https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
longwE <- paste(5*(0:4),"\u00B0E",sep="")
axis(1, at=c(5*(0:4)), labels=longwE,cex.axis = 0.6)
lattwN <- paste((53:59),"\u00B0N",sep="")
axis(2, at=c(53:59), labels=lattwN, las=1,cex.axis = 0.6)
#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 0.6, col=scales::alpha(colra,0.7))
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(5,54,z=legend.z,round=0,maxradius=(1*0.6),bty="n",txt.cex=0.6)
text(5,55,"samples",cex=0.6) 
#add legend to plot
legend("topright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.8, pt.cex = 1.2, ncol=3,
       x.intersp = 0.6,y.intersp = 0.9)

# Next - start plot for Europe
# begin plot, with defined borders
plot(NA,NA, xlim = c(-12, 16),
     ylim = c(34.8, 58.0),
     # use las to turn the tick labels to be horizontal
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     xaxt="n", yaxt="n",
     xlab="Longitude", ylab="Latitude",
     cex.lab = 0.7)
#xlab = NA, ylab = NA)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          # xlim = c(-12, 16),
          # ylim = c(34.8, 58.0),
          xlim = c(-48, 48),
          ylim = c(20, 70.0),
          #xlab = "Longitude", ylab = "Latitude",
          col="azure3", bg=transp_col)
#par(new=T)
mtext("b", side=3, adj=0, line=0.4, cex=1.2, font=2)
#dev.off()

mapplots::draw.pie(df_hl08$dec_lat,
                   df_hl08$dec_lon,
                   ma,
                   radius = df_hl08$rws, add = TRUE,
                   col=c(colra)
)
# https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
#longwE <- paste(seq(-12, 16,4),"\u00B0E",sep="")
longwW <- paste(c(16,12,8,4),"\u00B0W",sep="")
longwW <- paste(c(20,15,10,5),"\u00B0W",sep="")
longw0 <- paste("0","\u00B0",sep="")
longwE <- paste(c(4,8,12,16),"\u00B0E",sep="")
longwE <- paste(c(seq(5,20,5)),"\u00B0E",sep="")
longwW0E<- c(longwW,longw0,longwE)
#axis(1, at=c(-16,-12,-8,4,0,4,8,12,16), labels=longwW0E,cex.axis = 0.6)

axis(1, at=c(seq(-20,20, by=5)), labels=longwW0E,cex.axis = 0.6)
lattwN <- paste(round(seq(34.8, 58.0,4),0),"\u00B0N",sep="")
axis(2, at=c(round(seq(34.8, 58.0,4),0)), labels=lattwN, las=1,cex.axis = 0.6)

#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 1.6, col=scales::alpha(colra,0.7))
#class(mspdf)
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(-16,55,z=legend.z*1.875,round=0,maxradius=(1.6*1.875),bty="n",txt.cex=0.6)
text(-16,51,"samples",cex=0.6) 
#add legend to plot
legend("bottomleft", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.8, pt.cex = 1.2, ncol=3,
       x.intersp = 0.6,y.intersp = 0.9)

# end plot
dev.off()


#______________________________________________________________
#___________end plot map halpotype v 10____________________
#______________________________________________________________



#______________________________________________________________
#___________start plot map halpotype v 10____________________
#______________________________________________________________
# define name for outpu plot file
figname01 <- "Fig04_v11_map_haplotype_pie.pdf"
figname01 <- "Fig04_v12_map_haplotype_pie.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 3600, height = 4800, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 2), xpd=FALSE)
par(oma=c(0, 0, 0, 0))
#reset this parameter
par(mfrow = c(2, 1)) 


# begin plot, with defined borders
plot(NA,NA, xlim = c(4, 16), ylim = c(53, 59),
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     xaxt="n", yaxt="n",
     xlab="Longitude", ylab="Latitude",
     cex.lab = 0.7)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          xlim = c(4, 16), ylim = c(53, 59),
          xlab = "Longitude", ylab = "Latitude",
          col="azure3", bg=transp_col)
#add new title
mtext("a", side=3, adj=0, line=0.4, cex=1.2, font=2)
#dev.off()
mapplots::draw.pie(df_hl05$dec_lon,
                   df_hl05$dec_lon,
                   ma,
                   radius = df_hl05$rws, add = TRUE,
                   col=c(colra)
)
#https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
longwE <- paste(5*(0:4),"\u00B0E",sep="")
axis(1, at=c(5*(0:4)), labels=longwE,cex.axis = 0.6)
lattwN <- paste((53:59),"\u00B0N",sep="")
axis(2, at=c(53:59), labels=lattwN, las=1,cex.axis = 0.6)
#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 0.6, col=scales::alpha(colra,0.7))
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(5,55,z=legend.z,round=0,maxradius=(1*0.6),bty="n",txt.cex=0.6)
text(5,55.86,"samples",cex=0.6) 
#add legend to plot
legend("topright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.8, pt.cex = 1.2, ncol=4,
       x.intersp = 0.6,y.intersp = 0.9)

# Next - start plot for Europe
# begin plot, with defined borders
plot(NA,NA, xlim = c(-80, 60),
     ylim = c(6, 60.0),
     # use las to turn the tick labels to be horizontal
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     #xaxt="n", yaxt="n",
     xlab="Longitude", ylab="Latitude",
     cex.lab = 0.7)
#xlab = NA, ylab = NA)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          # xlim = c(-12, 16),
          # ylim = c(34.8, 58.0),
          xlim = c(-120, 80),
          ylim = c(-10, 90.0),
          #xlab = "Longitude", ylab = "Latitude",
          col="azure3", bg=transp_col)
#par(new=T)
mtext("b", side=3, adj=0, line=0.4, cex=1.2, font=2)
#dev.off()

mapplots::draw.pie(df_hl08$dec_lat,
                   df_hl08$dec_lon,
                   ma,
                   radius = df_hl08$rws, add = TRUE,
                   col=c(colra)
)
# https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
#longwE <- paste(seq(-12, 16,4),"\u00B0E",sep="")
longwW <- paste(c(16,12,8,4),"\u00B0W",sep="")
longwW <- paste(c(20,15,10,5),"\u00B0W",sep="")
longw0 <- paste("0","\u00B0",sep="")
longwE <- paste(c(4,8,12,16),"\u00B0E",sep="")
longwE <- paste(c(seq(5,20,5)),"\u00B0E",sep="")
longwW0E<- c(longwW,longw0,longwE)
#axis(1, at=c(-16,-12,-8,4,0,4,8,12,16), labels=longwW0E,cex.axis = 0.6)

# axis(1, at=c(seq(-80,60, by=10)), labels=longwW0E,cex.axis = 0.6)
# lattwN <- paste(round(seq(6, 6,10),0),"\u00B0N",sep="")
# axis(2, at=c(round(seq(6, 60,10),0)), labels=lattwN, las=1,cex.axis = 0.6)

#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 1.6, col=scales::alpha(colra,0.7))
#class(mspdf)
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(-40,40,z=legend.z*1.875,round=0,maxradius=(1.6*1.875),bty="n",txt.cex=0.6)
text(-40,44,"samples",cex=0.6) 
#add legend to plot
legend("bottomright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.8, pt.cex = 1.2, ncol=4,
       x.intersp = 0.6,y.intersp = 0.9)

# end plot
dev.off()


#______________________________________________________________
#___________end plot map halpotype v 10____________________
#______________________________________________________________

#______________________________________________________________
#___________start plot map halpotype v 13____________________
#______________________________________________________________

#______________________________________________________________
# define name for output plot file
figname01 <- "Fig04_v13_map_haplotype_pie.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 2), xpd=FALSE)
par(oma=c(0, 0, 0, 0))
#reset this parameter
par(mfrow = c(1, 1)) 


# begin plot, with defined borders
plot(NA,NA, xlim = c(4, 16), ylim = c(53, 59),
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     xaxt="n", yaxt="n",
     xlab="Longitude", ylab="Latitude",
     cex.lab = 1.2)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          xlim = c(4, 16), ylim = c(53, 59),
          xlab = "Longitude", ylab = "Latitude",
          col="azure3", bg=transp_col)
#add new title
#mtext("a", side=3, adj=0, line=0.4, cex=1.2, font=2)
#dev.off()
mapplots::draw.pie(df_hl05$dec_lon,
                   df_hl05$dec_lon,
                   ma,
                   radius = df_hl05$rws, add = TRUE,
                   col=c(colra)
)
#https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
longwE <- paste(5*(0:4),"\u00B0E",sep="")
axis(1, at=c(5*(0:4)), labels=longwE,cex.axis = 1.1)
lattwN <- paste((53:59),"\u00B0N",sep="")
axis(2, at=c(53:59), labels=lattwN, las=1,cex.axis = 1.1)
#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 0.6, col=scales::alpha(colra,0.7))
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(5,54.4,z=legend.z,round=0,maxradius=(1*0.6),bty="n",txt.cex=1.1)
text(5,55.2,"samples",cex=1.1) 
#add legend to plot
legend("topright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.8, pt.cex = 1.2, ncol=3,
       x.intersp = 0.6,y.intersp = 0.9)

# end plot
dev.off()
#______________________________________________________________
#___________end plot map halpotype v 13____________________
#______________________________________________________________



# define name for outpu plot file
figname01 <- "Fig05_v01_map_haplotype_pie.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 2), xpd=FALSE)
par(oma=c(0, 0, 0, 0))
#reset this parameter
par(mfrow = c(1, 1)) 
# Next - start plot for Europe and the Atlantic Ocean
# begin plot, with defined borders
plot(NA,NA, xlim = c(-80, 60),
     ylim = c(6, 60.0),
     # use las to turn the tick labels to be horizontal
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     xaxt="n", yaxt="n",
     xlab="Longitude", ylab="Latitude",
     cex.lab = 1.2)
#xlab = NA, ylab = NA)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          # xlim = c(-12, 16),
          # ylim = c(34.8, 58.0),
          xlim = c(-120, 80),
          ylim = c(-10, 90.0),
          #xlab = "Longitude", ylab = "Latitude",
          col="azure3", bg=transp_col)
#mtext("a", side=3, adj=0, line=0.4, cex=1.2, font=2)
mapplots::draw.pie(df_hl08$dec_lat,
                   df_hl08$dec_lon,
                   ma,
                   radius = df_hl08$rws, add = TRUE,
                   col=c(colra)
)
# https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
#longwE <- paste(seq(-12, 16,4),"\u00B0E",sep="")
longwW <- paste(c(80,60,40,20),"\u00B0W",sep="")
longw0 <- paste("0","\u00B0",sep="")
#longwE <- paste(c(4,8,12,16),"\u00B0E",sep="")
longwE <- paste(c(seq(20,60,20)),"\u00B0E",sep="")
longwW0E<- c(longwW,longw0,longwE)
#axis(1, at=c(-16,-12,-8,4,0,4,8,12,16), labels=longwW0E,cex.axis = 0.6)

axis(1, at=c(seq(-80,60, by=20)), labels=longwW0E,cex.axis = 0.8)
lattwN <- paste(round(seq(10, 60,10),0),"\u00B0N",sep="")
axis(2, at=c(round(seq(6, 60,10),0)), labels=lattwN, las=1,cex.axis = 0.8)

#add pies to map
draw.pie(xyz.Hpt$x, xyz.Hpt$y, xyz.Hpt$z, radius = 1.6, col=scales::alpha(colra,0.7))
#class(mspdf)
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt$z,na.rm=TRUE))/10^0,0)
legend.bubble(-40,40,z=legend.z*1.875,round=0,maxradius=(1.6*1.875),bty="n",txt.cex=0.8)
text(-40,44,"samples",cex=0.8) 
#add legend to plot
legend("bottomright", inset=c(0.03, 0), legend=c(Hpts), pch=c(rep(22,nHpt)), 
       bg="white",title="Haplotype", pt.bg=colra, cex=0.8, pt.cex = 1.2, ncol=4,
       x.intersp = 0.6,y.intersp = 0.9)

# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 


#______________________________________________________________
#___________end plot map halpotype v 13____________________
#______________________________________________________________





#_______________________________________________________________________________
# end - make plot for haplotypes in Europe
#_______________________________________________________________________________
