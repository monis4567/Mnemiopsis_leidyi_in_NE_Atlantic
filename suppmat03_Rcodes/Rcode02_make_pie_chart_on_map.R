#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#remove everything in the working environment, without a warning!!
rm(list=ls())
#define path for input dir
# wd01 <- "/suppmat01_inp_files"
# wd05 <- "/suppmat05_out_files"
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
wd_out01 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis"
#define path for input dir
wd01 <- "/suppmat01_inp_files"
wd05 <- "/suppmat05_out_files"
wd00_wd01 <- paste0(wd00,wd01)
wd00_wd05 <- paste0(wd00,wd05)
#https://stackoverflow.com/questions/21676721/r-plot-circular-histograms-rose-diagrams-on-map
if(!require(shapefiles)){
  install.packages("shapefiles")
}
library(shapefiles)
if(!require(mapplots)){
  install.packages("mapplots")
}
library(mapplots)
## install all packages required
if(!require(mapdata)){
  install.packages("mapdata")
}
library(mapdata)

if(!require(sp)){
  install.packages("sp")
}
library(sp)

if(!require(ggmap)){
  install.packages("ggmap")
}
library(ggmap)

if(!require(mapproj)){
  install.packages("mapproj")
}
library(mapproj)

if(!require(rworldmap)){
  install.packages("rworldmap")
}
library(rworldmap)

if(!require(akima)){
  install.packages("akima")
}
library(akima)

## install the package 'scales', which will allow you to make points on your plot more transparent
if(!require(scales)){
  install.packages("scales")
}
library(scales)
## install  package
if(!require(fields)){
  install.packages("fields")
}
library(fields)

# ## install the package 'marmap', which will allow you to plot bathymetric maps
# if(!require(marmap)){
#   #install.packages("marmap", repos='http://cran.us.r-project.org')
#   library(marmap)
# }
#get the package that enables the function 'subplot'
if(!require(TeachingDemos)){
  install.packages("TeachingDemos")

}
library(TeachingDemos)

if(!require(rgdal)){
  install_version("rgdal","1.6-7")
}
library(rgdal)
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
loc04fl <- paste(wd00_wd05,"/loc04.csv",sep="")
#read in csv
df_hap_loc04 <- read.table(loc04fl, sep=",")
#change column names
 colnames(df_hap_loc04) <- df_hap_loc04[1,]
# # drop first row with column names
 df_hap_loc04 <- df_hap_loc04[-1,]
 df_hap_loc08 <- df_hap_loc04
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
# change column name on first element of column names
colnames(df_hl05)[1] <- "smplloca"
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
# change roman numerals to arabian numerals
df_hl06$Hpt <- as.numeric(df_hl06$Hpt)
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
#dev.off()



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
df_hl07$lnNm <- NULL
# make all columns numeric
df_hl08 <- as.data.frame(lapply(df_hl08,as.numeric))
#make it a 'matrix array'
ma <- as.matrix(df_hl08)
# check if map exists , if not then get it
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
#df_hl08$Cnt[is.na(df_hl08$Cnt)]
df_hl08$Cnt <- as.numeric(df_hl08$Cnt)
# replace any NAs with zero
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
#dev.off()
# get library for package
library(shapefiles)


#_______________________________________________________________________________
# start - subset DNAbin object and plot only haplotypes 
# for Denmark- Germany collected in 2017-2018
#_______________________________________________________________________________

# Read in package libraries
library(pegas)
library(phangorn)
library(ape)
library(scales) # to be able to use alpha for setting opacity in colours
#define input working directories
wd01 <- "suppmat01_inp_files"
wd05 <- "suppmat05_out_files"
setwd(wd00)
#define paths for in and output directories
wd00_wd01 <- paste0(wd00,"/",wd01)
wd00_wd05 <- paste0(wd00,"/",wd05)
#read in the data  - first define the input file name
inf01 <- "algn_Mnelei_18s_10.aligned.fasta.fas"
inf01 <- "algn_Mnelei_v13.fa"
inf01 <- "algn_Mnelei_v15.fa"
inf01 <- "algn_Mnelei_18s_17.fas.aligned.fasta"
inf02 <- "df_clo2.csv" 
inf03 <- "df_cll_colf_smpl_loc.csv" 
inf04 <- "df_lN02.csv"
inf05 <- "df_dd02.csv"
inf06 <- "df_cfy3.csv"
#  then paste directory path and file name together
wd00_wd05_inf01 <- paste0(wd00_wd05,"/",inf01)
# paste together path for file with table with colors for locations.
wd00_wd05_inf02 <- paste0(wd00_wd05,"/",inf02)
wd00_wd05_inf03 <- paste0(wd00_wd05,"/",inf03)
wd00_wd05_inf04 <- paste0(wd00_wd05,"/",inf04)
wd00_wd05_inf05 <- paste0(wd00_wd05,"/",inf05)
wd00_wd05_inf06 <- paste0(wd00_wd05,"/",inf06)
# read in file with locations and matching colors for sample locations
df_clo2 <- read.csv(wd00_wd05_inf02)
df_cll <- read.csv(wd00_wd05_inf03) 
df_lN02 <- read.csv(wd00_wd05_inf04) 
df_dd02 <- read.csv(wd00_wd05_inf05) 
df_cfy3 <- read.csv(wd00_wd05_inf06) 
# add long names to  data frame with abbreviated collection location names
df_cll$lngNm <- c("BalticSea","CaspianSea","CentralWAtlantic","Mediterranean","NEAtlantic","AtlanticOcean:NWAtlantic","FynBogense,","FynKerteminde","NWGermanyWaddenSeaBussumHaupstr","NWGermanyNSeaHelgolandRds","NGermanyKielFjord","NGermanyMecklenburgerBuchtWismarBucht","NJyllandLimfjord","JyllandMariagerfjord","SamsoeBallen","SjaellandSkovshoved")
# mathc to get latitude and longitude for Danish and German sampling locations
df_cll$dec_lat <- df_clo2$dec_lat[match(df_cll$lngNm,df_clo2$locality3)]
df_cll$dec_lon <- df_clo2$dec_lon[match(df_cll$lngNm,df_clo2$locality3)]
# also get the decimal latitude and longitude for the non- Danish and German sampling locations
df_cll$dec_lon[is.na(df_cll$dec_lon)] <- df_lN02$declon2[match(df_cll$lngNm[is.na(df_cll$dec_lon)], df_lN02$location)]
df_cll$dec_lat[is.na(df_cll$dec_lat)] <- df_lN02$declat2[match(df_cll$lngNm[is.na(df_cll$dec_lat)], df_lN02$location)]
# remove first column
df_dd02 <- df_dd02[,-1]
# then read in the fasta file which makes it a DNAbinxn object
dnb_pip4 <- read.dna(file=wd00_wd05_inf01,format="fasta")
#make the DNAbin object a genind object
geni_pip4 <- adegenet::DNAbin2genind(dnb_pip4)
#make the genind object a dataframe
df_pip4 <- adegenet::genind2df(geni_pip4)

# subset the DNAbin object and make a haplotype network only on Danish- German
# samples from 2017-2018
#_______________________________________________________________________________
# make the dnabin object a phydata object
phyd_pip4 <- phangorn::as.phyDat(dnb_pip4)
# strsplit names and turn into characters
# and rowbind the nested lists to a dataframe
df_pp4Nm <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(names(phyd_pip4)),
                                  "_")))
# modify columns names
colnames(df_pp4Nm) <- c("smplNm","loc","yer")
# get unique location names
uloc <- unique(df_pp4Nm$loc)
# order the unique location names
uloc <- uloc[order(uloc)]
# make a vector that categorizes the overall sampling locations for :
# North East Atlantic, Caspian Sea, Central Western Atlantic and Mediterranean
ovloc <- c("NEA", "CS", "CWA", "NEA", "NEA", "NEA", "NEA", "NEA", "NEA", 
           "NEA", "NEA", "M", "NEA", "NWA", "NEA", "NEA")
# combine in to a data frame
df_ovl <- as.data.frame(cbind(uloc, ovloc))
# subset to only comprise the NEA overall location names
df_ovlNEA <- subset(df_ovl, ovloc=="NEA")
# get unique NEA locatoins and all full sample names
NEAloc <- df_ovlNEA$uloc
Nmsloc <- names(phyd_pip4)
# use grepl instead of dplyr, to get all NEA full sample names
# but first collapse and paste with '|' sign to be able to grep for
# multiple chracteristics
# https://stackoverflow.com/questions/45357806/dplyr-select-and-starts-with-on-multiple-values-in-a-variable-list?noredirect=1&lq=1
NmslocNEA <- Nmsloc[grepl(paste(NEAloc, collapse="|"),Nmsloc)]
# now subset the phydat object by the full names that match bein NEA samples
phyd_pip5 <- subset(phyd_pip4, NmslocNEA)
# limit to only include samples collected in 2017 and 2018
Nmswy17_18 <- NmslocNEA[grepl(paste(c("2017","2018"), collapse="|"),NmslocNEA)]
# subset the phy dat object - I am replacing the previous object
phyd_pip5 <- subset(phyd_pip4, Nmswy17_18)
# make the phy-dat object a DNAbin object instead
dnb_pip5 <- as.DNAbin(phyd_pip5)
# make the DNAbin object a distance matrix 
dst_pip5 <- ape::dist.dna(dnb_pip5, model= "raw")
# make it a matrix
mtx_pip5 <- as.matrix(dst_pip5)
# make it a data frame
df_pip5 <- as.data.frame(mtx_pip5)
# split the row name string by a character 
lpip5 <- strsplit(as.character(row.names(mtx_pip5)), "_")
# get the 2 nd element per vector in the list - this holds the abbreviated 
# location name
grp.loc5 <- sapply(lpip5, "[[", 2)
#make haplotype object
ht5 <- pegas::haplotype(dnb_pip5)
#which polyps belong to which haplotype?
ht5_indices<-attr(ht5, "index")
# get number of labels 
hlab5 <- length(labels(ht5))
#assign labels
names(ht5_indices)<-c(1:hlab5)
# make roman numerals arabic numerals instead
labs.ht5 <- as.numeric(as.roman(labels(ht5)))
# use arabian numerals instead of roman numerals for haplotypes
ht5 <- haplotype(dnb_pip5,labels=c(labs.ht5))
# make haplonet object
hN5 <- pegas::haploNet(ht5)
#prepare hpt table
ind.hap5 <- with(
  utils::stack(setNames(attr(ht5, "index"), rownames(ht5))),
  table(hap=ind, pop=rownames(dnb_pip5)[values]))
#make it a dataframe
df_ihpt05 <- as.data.frame(ind.hap5)
# make haplotype labels arabic numerals instead of roman numerals 
df_ihpt05$hap.ab <- as.numeric(as.roman(as.character(df_ihpt05$hap)))
# split the string with sequence name which holds the location name
lbf5  <- strsplit(as.character(df_ihpt05$pop), "_")
# copy table
df_ihpt06 <- df_ihpt05
# get the second element from this string split
df_ihpt05$pop.loc <- sapply(lbf5, "[[", 2)
df_ihpt06$pop.yer <- sapply(lbf5, "[[", 3)
#remove any zero occurences
df_ihpt05 <- df_ihpt05[df_ihpt05$Freq >= 1,]	
df_ihpt06 <- df_ihpt06[df_ihpt06$Freq >= 1,]	
# make a table of locations and haplotype numbers
hsl5 <- table(df_ihpt05$hap.ab, df_ihpt05$pop.loc)
hsl6 <- table(df_ihpt06$hap.ab, df_ihpt06$pop.yer)
# get latitude and longtiude
df_ihpt05$dec_lat <-  df_cll$dec_lat[match(df_ihpt05$pop.loc,df_cll$coll_loc)]
df_ihpt05$dec_lon <-  df_cll$dec_lon[match(df_ihpt05$pop.loc,df_cll$coll_loc)]

#______________________________________________________________
#___________start plot map halpotype v 14____________________
#______________________________________________________________
# make xyz list 
xyz.Hpt5 <- mapplots::make.xyz(df_ihpt05$dec_lon,df_ihpt05$dec_lat,df_ihpt05$Freq,df_ihpt05$hap.ab)
length(xyz.Hpt5$x)==length(xyz.Hpt5$y) & length(xyz.Hpt5$y)==nrow(xyz.Hpt5$z)
#sort the haplotypes increasingly by number
sHpt5 <- unique(df_ihpt05$hap.ab)[order(unique(df_ihpt05$hap.ab))]
# count number of haplotypes
nHpt5 <- length(unique(df_ihpt05$hap.ab))
# make a colour range to reflect the number of haplotypes
colra5 <- rainbow(nHpt5)
# try making a different colour range
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
# make a colour ramp
colfunc <- colorRampPalette(cbbPalette2)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# use the colour function
cl05 <- colfunc(nHpt5)
# Replace the colour range with rainbow colours
colra5 <- cl05
# define name for output plot file
figname01 <- "Fig04_v14_map_haplotype_pie.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
jpeg(pthfignm01,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))
#add extra space to the right of the plot
par(mar=c(6, 6, 2, 2), xpd=FALSE)
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
     #xlab="Longitude", ylab="Latitude",
     xlab="", ylab="",
     cex.lab = 1.2)
title(xlab = "Longitude", ylab = "Latitude", line = 5)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          xlim = c(4, 16), ylim = c(53, 59),
          xlab = "Longitude", ylab = "Latitude",
          col="azure3", bg=transp_col)
#https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
longwE <- paste(5*(0:4),"\u00B0E",sep="")
axis(1, at=c(5*(0:4)), labels=longwE,cex.axis = 1.1)
lattwN <- paste((53:59),"\u00B0N",sep="")
axis(2, at=c(53:59), labels=lattwN, las=1,cex.axis = 1.1)
#add pies to map
# the 'draw.pie' function requires the largest pie diameter as input value for
# the radius of the largest pie
lrgst_pie <- 0.4*sqrt(max(rowSums(xyz.Hpt5$z))/pi)
#add pies to map
draw.pie(xyz.Hpt5$x, xyz.Hpt5$y, sqrt(xyz.Hpt5$z/pi), 
         radius = lrgst_pie, col=scales::alpha(colra5,0.7))
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt5$z,na.rm=TRUE)),0)
legend.bubble(5,54.4,legend.z,round=0,
              maxradius=lrgst_pie,bty="n",txt.cex=1.1)
text(5,55.2,"samples",cex=1.1) 
#add legend to plot
legend("topright", inset=c(0.03, 0), legend=c(sHpt5), pch=c(rep(22,nHpt5)), 
       bg="white",title="Haplotype", pt.bg=colra5, cex=0.96, pt.cex = 1.2, ncol=3,
       x.intersp = 0.6,y.intersp = 0.9)

# end plot
dev.off()
#______________________________________________________________
#___________end plot map halpotype v 14____________________
#______________________________________________________________
#define columns to keep
cke <- c("hap","Freq","pop.loc")
df_i05 <- df_ihpt05[cke]
# summarize haplotypes
library(dplyr)
df_i05 <- df_i05 %>% 
  dplyr::group_by(hap,pop.loc) %>% 
  dplyr::summarise(Freq2 = sum(Freq))
library(dplyr)
# use tidyr::pivot_wider to get the data frame re arranged
df_i05  <- df_i05 %>% tidyr::pivot_wider(names_from = "pop.loc",
                              values_from = "Freq2")
# convert to data frame
df_i05 <- as.data.frame(df_i05)
# make all values characters
df_i05 <- sapply(df_i05, as.character) 
df_i05[is.na(df_i05)] <- "0"
# make it a data frame
df_i05 <- as.data.frame(df_i05)
# make the entire data frame numeric
df_i05[] <- lapply(df_i05, as.numeric)
# sum up all columns
tot.hp <- colSums(df_i05)
# modify the total of haplotype nuumbers
tot.hp[1] <- "Total"
# bind the row with the total
df_i05 <- rbind(df_i05,tot.hp)
# replcae 0 with nothing
df_i05[(df_i05==0)] <- ""
# write out a csv file
write.csv(df_i05,
          file = paste0(wd00_wd05,"/Table02_v01_haplotypes_freq_DK_Germany.csv"))
#_______________________________________________________________________________
# end - subset DNAbin object and plot only haplotypes 
# for Denmark- Germany collected in 2017-2018
#_______________________________________________________________________________


#_______________________________________________________________________________
# start - plot the DNAbin object with all  haplotypes 
# for both Denmark- Germany samples collected in 2017-2018 and 
# NCBI GenBank samples
#_______________________________________________________________________________


# make the phy-dat object a DNAbin object instead
dnb_pip7 <- as.DNAbin(phyd_pip4)
# make the DNAbin object a distance matrix 
dst_pip7 <- ape::dist.dna(dnb_pip7, model= "raw")
# make it a matrix
mtx_pip7 <- as.matrix(dst_pip7)
# make it a data frame
df_pip7 <- as.data.frame(mtx_pip5)
# split the row name string by a character 
lpip7 <- strsplit(as.character(row.names(mtx_pip7)), "_")
# get the 2 nd element per vector in the list - this holds the abbreviated 
# location name
grp.loc7 <- sapply(lpip7, "[[", 2)
#make haplotype object
ht7 <- pegas::haplotype(dnb_pip7)

#which polyps belong to which haplotype?
ht7_indices<-attr(ht7, "index")
# get number of labels 
hlab7 <- length(labels(ht7))
#assign labels
names(ht7_indices)<-c(1:hlab7)
# make roman numerals arabic numerals instead
labs.ht7 <- as.numeric(as.roman(labels(ht7)))
# use arabian numerals instead of roman numerals for haplotypes
ht7 <- haplotype(dnb_pip7,labels=c(labs.ht7))
# make haplonet object
hN7 <- pegas::haploNet(ht7)
#prepare hpt table
ind.hap7 <- with(
  utils::stack(setNames(attr(ht7, "index"), rownames(ht7))),
  table(hap=ind, pop=rownames(dnb_pip7)[values]))
#make it a dataframe
df_ihpt07 <- as.data.frame(ind.hap7)
# make haplotype labels arabic numerals instead of roman numerals 
df_ihpt07$hap.ab <- as.numeric(as.roman(as.character(df_ihpt07$hap)))
# split the string with sequence name which holds the location name
lbf7  <- strsplit(as.character(df_ihpt07$pop), "_")
# copy table
df_ihpt08 <- df_ihpt07
# get the second element from this string split
df_ihpt07$pop.loc <- sapply(lbf7, "[[", 2)
df_ihpt08$pop.yer <- sapply(lbf7, "[[", 3)
#remove any zero occurences
df_ihpt07 <- df_ihpt07[df_ihpt07$Freq >= 1,]	
df_ihpt08 <- df_ihpt08[df_ihpt08$Freq >= 1,]	
# make a table of locations and haplotype numbers
hsl7 <- table(df_ihpt07$hap.ab, df_ihpt07$pop.loc)
hsl8 <- table(df_ihpt08$hap.ab, df_ihpt08$pop.yer)
# get latitude and longtiude
df_ihpt07$dec_lat <-  df_cll$dec_lat[match(df_ihpt07$pop.loc,df_cll$coll_loc)]
df_ihpt07$dec_lon <-  df_cll$dec_lon[match(df_ihpt07$pop.loc,df_cll$coll_loc)]


#xyz.Hpt7
# make xyz list 
xyz.Hpt7 <- mapplots::make.xyz(df_ihpt07$dec_lon,df_ihpt07$dec_lat,df_ihpt07$Freq,df_ihpt07$hap.ab)
length(xyz.Hpt7$x)==length(xyz.Hpt7$y) & length(xyz.Hpt7$y)==nrow(xyz.Hpt7$z)
#sort the haplotypes increasingly by number
sHpt7 <- unique(df_ihpt07$hap.ab)[order(unique(df_ihpt07$hap.ab))]
# count number of haplotypes
nHpt7 <- length(unique(df_ihpt07$hap.ab))
# make a colour range to reflect the number of haplotypes
colra7 <- rainbow(nHpt7)
# try making a different colour range
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
# make a colour ramp
colfunc <- colorRampPalette(cbbPalette2)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# use the colour function
cl07 <- colfunc(nHpt7)
# Replace the colour range with rainbow colours
colra7 <- cl07
# define name for output plot file
figname01 <- "Fig04_v15_map_haplotype_pie.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
jpeg(pthfignm01,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))
#add extra space to the right of the plot
par(mar=c(6, 6, 2, 2), xpd=FALSE)
par(oma=c(0, 0, 0, 0))
#reset this parameter
par(mfrow = c(1, 1)) 
# begin plot, with defined borders
plot(NA,NA, xlim = c(-80, 60),
     ylim = c(6, 60.0),
     # use las to turn the tick labels to be horizontal
     las=1,
     #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
     asp=1.6,
     #surpress tickmarks
     xaxt="n", yaxt="n",
     #xlab="Longitude", ylab="Latitude",
     xlab="", ylab="",
     #line= 4,
     cex.lab = 1.2)
title(xlab = "Longitude", ylab = "Latitude", line = 5)
# add internal map inside coastline
maps::map('worldHires', add=TRUE, fill=TRUE, 
          xlim = c(-120, 80),
          ylim = c(-10, 90.0),
          #xlab = "Longitude", ylab = "Latitude",
          
          col="azure3", bg=transp_col)

#https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
longwW <- paste(c(80,60,40,20),"\u00B0W",sep="")
longw0 <- paste("0","\u00B0",sep="")
#longwE <- paste(c(4,8,12,16),"\u00B0E",sep="")
longwE <- paste(c(seq(20,60,20)),"\u00B0E",sep="")
longwW0E<- c(longwW,longw0,longwE)


axis(1, at=c(seq(-80,60, by=20)), labels=longwW0E,cex.axis = 1.0)
lattwN <- paste(round(seq(10, 60,10),0),"\u00B0N",sep="")
axis(2, at=c(round(seq(6, 60,10),0)), labels=lattwN, las=1,cex.axis = 1.0)
# the 'draw.pie' function requires the largest pie diameter as input value for
# the radius of the largest pie
lrgst_pie <- sqrt(max(rowSums(xyz.Hpt7$z))/pi)
#add pies to map
draw.pie(xyz.Hpt7$x, xyz.Hpt7$y, sqrt(xyz.Hpt7$z/pi), 
         radius = lrgst_pie, col=scales::alpha(colra7,0.7))
#_______add histogram bars to positions on the map
# add a legend
legend.z <- round(max(rowSums(xyz.Hpt7$z,na.rm=TRUE)),0)
#write.csv2(as.data.frame(as.matrix(xyz.Hpt7$z)),paste0(wd00_wd05,"/xyz.Hpt7.csv"))
# also adjust the circle size for the legend , so that it is based on the
# square root of the max size divided by pi
legend.bubble(-40,40,z=legend.z,round=0,
              maxradius=lrgst_pie,bty="n",txt.cex=1.1)
text(-40,44.2,"samples",cex=1.1) 
#add legend to plot
legend("bottomright", inset=c(0.03, 0), legend=c(sHpt7), pch=c(rep(22,nHpt7)),
       bg="white",title="Haplotype", pt.bg=colra7, cex=0.96, pt.cex = 1.2, ncol=6,
       x.intersp = 0.6,y.intersp = 0.9)

# end plot
dev.off()
#__________________________________________________________________
#define columns to keep
cke <- c("hap","Freq","pop.loc")
df_i07 <- df_ihpt07[cke]
# summarize haplotypes
library(dplyr)
df_i07 <- df_i07 %>% 
  dplyr::group_by(hap,pop.loc) %>% 
  dplyr::summarise(Freq2 = sum(Freq))
library(dplyr)
# use tidyr::pivot_wider to get the data frame re arranged
df_i07  <- df_i07 %>% tidyr::pivot_wider(names_from = "pop.loc",
                                         values_from = "Freq2")
# convert to data frame
df_i07 <- as.data.frame(df_i07)
# make all values characters
df_i07 <- sapply(df_i07, as.character) 
df_i07[is.na(df_i07)] <- "0"
# make it a data frame
df_i07 <- as.data.frame(df_i07)
# make the entire data frame numeric
df_i07[] <- lapply(df_i07, as.numeric)
# sum up all columns
tot.hp <- colSums(df_i07)
# modify the total of haplotype nuumbers
tot.hp[1] <- "Total"
# bind the row with the total
df_i07 <- rbind(df_i07,tot.hp)
# replcae 0 with nothing
df_i07[(df_i07==0)] <- ""
# write out a csv file
write.csv(df_i07,
          file = paste0(wd00_wd05,"/Table02_v02_haplotypes_freq_global.csv"))

#_______________________________________________________________________________
# end - plot the DNAbin object with all  haplotypes 
# for both Denmark- Germany samples collected in 2017-2018 and 
# NCBI GenBank samples
#_______________________________________________________________________________

