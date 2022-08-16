#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#remove everything in the working environment, without a warning!!
rm(list=ls())
# Read in package libraries
library(pegas)
library(ape)
packageVersion("pegas")
#[1] ‘1.1’
packageVersion("ape")
# [1] ‘5.6.2’
#install.packages("pegas")
R.Version()
# $platform
# [1] "x86_64-pc-linux-gnu"
# 
# $arch
# [1] "x86_64"
# 
# $os
# [1] "linux-gnu"
# 
# $system
# [1] "x86_64, linux-gnu"
# 
# $status
# [1] ""
# 
# $major
# [1] "4"
# 
# $minor
# [1] "1.2"
# 
# $year
# [1] "2021"
# 
# $month
# [1] "11"
# 
# $day
# [1] "01"
# 
# $`svn rev`
# [1] "81115"
# 
# $language
# [1] "R"
# 
# $version.string
# [1] "R version 4.1.2 (2021-11-01)"
# 
# $nickname
# [1] "Bird Hippie"

# Or the development version from GitHub:
# install.packages("devtools")
#devtools::install_github("r-lib/devtools")
if(!require(pegas)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  # install from GitHub
  devtools::install_github("emmanuelparadis/pegas/pegas")
}
# read in package libraries
library(dplyr)
library(pegas)
library(ape)
#define overall working directory
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
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
# remove first column
df_dd02 <- df_dd02[,-1]
# then read in the fasta file which makes it a DNAbinxn object
dnb_pip4 <- read.dna(file=wd00_wd05_inf01,format="fasta")
#make the DNAbin object a genind object
geni_pip4 <- adegenet::DNAbin2genind(dnb_pip4)
#make the genind object a dataframe
df_pip4 <- adegenet::genind2df(geni_pip4)
#get the row names
orig_rwnm <- row.names(df_pip4)
# replace all NAs
df_pip4 <- df_pip4 %>% replace(is.na(.), "-")
#Grep among the list of files, to get the file that holds the different names
fnm_clo3 <- list.files(wd00_wd05)[grepl("clo03",list.files(wd00_wd05))]
fnm_df_lN02 <- list.files(wd00_wd05)[grepl("df_lN02",list.files(wd00_wd05))]
#read locations file
df_clo3<-read.csv(file=paste0(wd00_wd05,"/",fnm_clo3),header=TRUE,sep=";")
df_lN02<-read.csv(file=paste0(wd00_wd05,"/",fnm_df_lN02),header=TRUE,sep=",")
# split string and get lists nested in a list
lbls01 <- strsplit(as.character(row.names(df_pip4)), "_")
# # split to get a list of vectors
lrNm <- strsplit(as.character(rownames(dnb_pip4)), "_")
#get second object in nested list
#unique sample number
smplno1 <- sapply(lrNm, "[[", 1)
smplno2 <- sapply(lrNm, "[[", 2)
smplno3 <- sapply(lrNm, "[[", 3)
# #create haplotypes from dna.bin
# pip4Haps <- pegas::haplotype(pip4)
# #prepare hpt table
# ind.hap<-with(
#   stack(setNames(attr(pip4Haps, "index"), rownames(pip4Haps))),
#   table(hap=ind, pop=rownames(pip4)[values]))
# #make it a dataframe
# df_ihpt01 <- as.data.frame(ind.hap)
# #limit to include only 'Freq' that equals 1
# df_ihpt02 <- df_ihpt01[df_ihpt01$Freq == 1,]
# # get sample numbers from haplotype table
# liht2Nm <- strsplit(as.character(df_ihpt02$pop), "_")
# # get the first element of the split string
# liht2Nm.1 <- sapply(liht2Nm, "[[", 2)
# # copy the vector
# liht2Nm.1a <- liht2Nm.1
# # grep for long (5) numeric characters in the string, and substitute, to only get accession number 
# liht2Nm.1a[grepl("[0-9]{5}",liht2Nm.1a)] <- gsub("Mnelei","",liht2Nm.1[grepl("[0-9]{5}",liht2Nm.1)])
# # get the sampling year
# yearsmpl4 <- sapply(liht2Nm, "[[", 3)
# # bind them by column in a new data frame
# df_smpl.yer4 <- as.data.frame(cbind(liht2Nm.1a,yearsmpl4))
# # change the column names
# colnames(df_smpl.yer4) <- c("smplNo5","smplyer5")
# #copy the vector
# smplno1a <- smplno1
# # replace in the copied vector
# smplno1a[grepl("Mnemiopsis",smplno1)] <- smplno3[grepl("Mnemiopsis",smplno1)]
# # match back to get sampling year
# smplno5 <- df_smpl.yer4$smplyer5[match(smplno1a,df_smpl.yer4$smplNo5)]
# # combine to a data frame
# df_rNm_lblNm2 <- cbind(smplno1a,df_rNm_lblNm, smplno5)
# # change column name
# colnames(df_rNm_lblNm2) <- c("longNm1","longNm2","locat","smplyear")
# # change row names in data frame with sequences
# rownames(df_pip4) <- paste0(df_rNm_lblNm2$longNm1,"_",df_rNm_lblNm2$locat,"_",df_rNm_lblNm2$smplyear)
# # grep for "Bolinopsis" to exclude this sequence from the data frame
# df_pip4 <- df_pip4[!grepl("Bolinopsis",rownames(df_pip4)),]
# #make the date frame a matrix and a DNAbin object again
# dnb_pip4 <- as.DNAbin(as.matrix(df_pip4))
# make a distance matrix 
dst_pip4 <- ape::dist.dna(dnb_pip4, model= "raw")
# make it a matrix
mtx_pip4 <- as.matrix(dst_pip4)
# make it a data frame
df_pip5 <- as.data.frame(mtx_pip4)
# split the row name string by a character 
lpip5 <- strsplit(as.character(row.names(mtx_pip4)), "_")
# get the 2 nd element per vector in the list - this holds the abbreviated 
# location name
grp.loc5 <- sapply(lpip5, "[[", 2)
#make haplotype object
ht4 <- pegas::haplotype(dnb_pip4)
#which polyps belong to which haplotype?
ht4_indices<-attr(ht4, "index")
# get number of labels 
hlab4 <- length(labels(ht4))
#assign labels
names(ht4_indices)<-c(1:hlab4)
# make roman numerals arabic numerals instead
labs.ht4 <- as.numeric(as.roman(labels(ht4)))
# use arabian numerals instead of roman numerals for haplotypes
ht4 <- haplotype(dnb_pip4,labels=c(labs.ht4))
# make haplonet object
hN4 <- pegas::haploNet(ht4)
#prepare hpt table
ind.hap4<-with(
  stack(setNames(attr(ht4, "index"), rownames(ht4))),
  table(hap=ind, pop=rownames(dnb_pip4)[values]))
#make it a dataframe
df_ihpt04 <- as.data.frame(ind.hap4)
# make haplotype labels arabic numerals instead of roman numerals 
df_ihpt04$hap.ab <- as.numeric(as.roman(as.character(df_ihpt04$hap)))
# split the string with sequence name which holds the location name
lbf4  <- strsplit(as.character(df_ihpt04$pop), "_")
# get the second element from this string split
df_ihpt04$pop.loc <- sapply(lbf4, "[[", 2)
#remove any zero occurences
df_ihpt04 <- df_ihpt04[df_ihpt04$Freq >= 1,]	
# make a table of locations and haplotype numbers
hsl4 <- table(df_ihpt04$hap.ab, df_ihpt04$pop.loc)
#make the plot
plot(hN4, 
     size = sqrt(attr(hN4,"freq")/pi), 
     #size = log10(attr(hN4,"freq")), 
     scale.ratio = 0.6, 
     cex = 1.1, # set size of roman numerals on circles for haplotype ID
     #bg= colfh, 
     pie = hsl4, 
     show.mutation = 2, threshold = 0, labels(T))

# Use the replot function to be able to move pies around in the plot
# run the out-commented line below
#xy <- replot()
# and rearrange your pies as you want them to be placed.
# NOTICE !!! This replot function does not work in Rstudio. You must start up R in a terminal, and then 
# run all the code above, including the replot part
# Once you are happy with the rearranged circles, you can right click the plot, and exit the replot
# function.
# Because you write the new positions of the pies to the object 'xy' you might want to 
# use the 'dput' function when using R in a terminal , like this: dput(xy)
# This will give you a list of the coordinates, that you can use here in Rstudio
xy <- list(x = c(-3.78009492656885, 0.151374570335722, 0, -7.70394570412442, 
                 -0.413177807763105, -9.72915900866923, -8.31150969548786, -9.37474668037389, 
                 -2.81912360628134, 0.270331682520766, -2.07856204801223, 2.19338982144833, 
                 3.22502895264294, 0.847543264652827, 3.03507968042298, -6.33718905319336, 
                 -1.68543764284715, 0.104992462245932, -3.0626070633811, -1.35783397187625, 
                 0.502954879238477, 4.00856970055032, -0.714548112740307, 1.01401660595308, 
                 -0.951605419872329, -5.05917141159912, -1.84268740491318, 4.85353162973206, 
                 1.95751517834928, -4.28002434887989, -3.62846625836468, -7.2048744662272, 
                 -1.82084143217233, 1.94441103151045, 2.34215255505902, -3.02090226700123, 
                 2.78307642919596, -1.79027081755784, -1.4004692937287, 4.20652972471708, 
                 7.07114931100761, 1.31541198324631, -5.05917141159912, 3.09757595332802, 
                 5.77756760859653, 5.8596172932963, 2.10166079357648, -1.06604271516307, 
                 -5.81397012860837, -6.69133905185202, -8.49890774248904, -2.48928377455822, 
                 0.804350256531708, -1.37093811871508, -1.34472982503741, -4.93742968304924, 
                 -4.10958592891006), y = c(0.185371219435836, -1.97038269146326, 
                                           0, -0.269105781977764, -0.799328950293629, -0.243857059677008, 
                                           1.37206116757134, 0.9428328884585, 0.798043387686757, -3.8796730709858, 
                                           0.46938451101758, 2.220731379872, 2.46000064746, -2.25524327550549, 
                                           4.07033174790751, 1.37206116757134, 2.71738198447251, 4.44850022968935, 
                                           1.81798983789992, 1.39733695645537, 2.56054495144077, 3.80983701107041, 
                                           2.65266532385822, -1.46493889637387, -1.37345062710535, 1.50229212711966, 
                                           -1.7655432096847, -0.184103126614666, -1.68712469316883, -1.50897834493825, 
                                           -1.35480084091025, -1.62359484098738, -3.038898020258, -0.510846945430789, 
                                           0.409492359034124, -0.900323839496651, 1.00524437387602, -0.419358676162275, 
                                           -2.33950101063971, 2.23665152919641, 3.57760074848688, -0.706893236720463, 
                                           -1.46040946635667, -2.45824010557488, 5.01345158824423, 3.13420484655794, 
                                           2.92649802851483, -3.16032021671195, -1.33898726990273, 0.463107164744145, 
                                           -1.38458801442011, -0.344851948880031, 1.88091780830323, -0.837590764246912, 
                                           0.691570307812544, 2.25510974513413, 1.91512759506308))


# #where "............." are your options (colours, etc...) Once you found a nice layout you can save it:
# saveRDS(xy, "mynicelayout.rds")
# #and reuse it in a script:
# xy <- readRDS("mynicelayout.rds")
# plot(h, ............, xy = xy)
# #'xy' is a simple list with two vectors, so you can also edit it manually.
#plot(h, ............, xy = xy)
plot(hN4, 
     size = sqrt(attr(hN4,"freq")/pi), 
     #size = log10(attr(hN4,"freq")), 
     scale.ratio = 0.6, 
     cex = 1.1, # set size of roman numerals on circles for haplotype ID
     #bg= colfh, 
     pie = hsl4, 
     show.mutation = 2, threshold = 0, labels(T), xy=xy)

# split the string with sequence name which holds the location name
lbf4  <- strsplit(as.character(df_ihpt04$pop), "_")
# get the second element from this string split
df_ihpt04$pop.loc <- sapply(lbf4, "[[", 2)

smplyear3 <- sapply(lbf4, "[[", 3)
usmplyear3 <- unique(smplyear3)
df_ihpt05 <- df_ihpt04
df_ihpt05$pop.loc <- sapply(lbf4, "[[", 3)
#remove any zero occurences
df_ihpt04 <- df_ihpt04[df_ihpt04$Freq >= 1,]	
df_ihpt05 <- df_ihpt05[df_ihpt05$Freq >= 1,]	
# make a table of locations and haplotype numbers
hsl4 <- table(df_ihpt04$hap.ab, df_ihpt04$pop.loc)
hsl5 <- table(df_ihpt05$hap.ab, df_ihpt05$pop.loc)


# get colours to use for pies
colfh<- df_cll$colfcol_loc[match(colnames(hsl4),df_cll$coll_loc)]
colfh_y <- df_cfy3$colour[match(colnames(hsl5),df_cfy3$smplyr)]

flnm <- c(paste("Fig02_v03_haplotype_network_",inf01,"02.jpg",  sep = ""))
#paste output directory and filename together in a string
outflnm <- paste(wd00_wd05,"/",flnm,sep="")
# Exporting PFD files via postscript()           
# pdf(outflnm,
#     width=(1*1.0*8.2677),height=(4*1.0*2.9232))
jpeg(outflnm,
     width=(3200),height=(4800),res=300)
#define lpot arrangement
tbt.par <- par(mfrow=c(2, 1),
               oma=c(0,0,0,0), #define outer margins
               mai=c(0,0,0,0), #define inner margins
               mar=c(0,0,2,0))
# plot the network with colours
plot(hN4, 
     size = sqrt(attr(hN4,"freq")/pi), 
     #size = log10(attr(hN4,"freq")), 
     scale.ratio = 0.6, 
     cex = 1.1, # set size of roman numerals on circles for haplotype ID
     #bg= colfh,
     # use the 'alpha' function on the color to make the color fill more transparent
     #bg=alpha(c(colfh),c(0.7)),
     bg=alpha(c(colfh),c(1.0)),
     pie = hsl4, 
     show.mutation = 2, threshold = 0, labels(T), xy=xy)
#add a legend to the plot
legend("bottomleft",colnames(hsl4), 
       # use the 'alpha' function on the color to make the color fill more transparent
       #pt.bg=alpha(c(colfh),c(0.7))
       pt.bg=alpha(c(colfh),c(1.0))
       ,box.col=NA,
       #col=rainbow(ncol(new.hap.smplloc)), 
       pt.lwd=0.4,
       pch=21, ncol=2, cex=1.2)
# add subfigure letter
title(main = "a",
      cex.main = 2.2,   font.main= 2, col.main= "black",
      adj = 0.01, line = 0.1)

# plot the network with colours
plot(hN4, 
     size = sqrt(attr(hN4,"freq")/pi), 
     #size = log10(attr(hN4,"freq")), 
     scale.ratio = 0.6, 
     cex = 1.1, # set size of roman numerals on circles for haplotype ID
     #bg= alpha(c(colfh_y),c(0.7)),
     # use the 'alpha' function on the color to make the color fill more transparent
     bg= alpha(c(colfh_y),c(1.0)), 
     pie = hsl5, 
     show.mutation = 2, threshold = 0, labels(T), xy=xy)
#add a legend to the plot
legend("bottomleft",colnames(hsl5), 
       #col=rainbow(ncol(new.hap.smplye)), 
       # use the 'alpha' function on the color to make the color fill more transparent
       #pt.bg=alpha(c(colfh_y),c(0.7))
       pt.bg=alpha(c(colfh_y),c(1.0))
       ,box.col=NA,pt.lwd=0.4,
       pch=21, ncol=2, cex=1.2)
# add subfigure letter
title(main = "b",
      cex.main = 2.2,   font.main= 2, col.main= "black",
      adj = 0.01, line = 0.1)
# use the par settings defined above
par(tbt.par)
# end file to save as
dev.off()  
#reset this par parameter
par(mfrow = c(1, 1)) 



# #_______________________________________________________________________________
# Another solution which works only with pegas: set the lengths of the links between haplotypes equal to the sum of the sizes of their respective symbols:
#   ## extract the size of the haplotype symbols:
#   size <- attr(net3, "freq")
# ## make a copy of the network:
# mynet <- net3
# ## change the lengths of the links by changing the 3th column:
# mynet[, 3] <- size[net3[, 1]] + size[net3[, 2]]
# We now plot the modified network:
#   o <- plot(mynet, size = size, cex = 0.8, pie = ind.hap3, threshold = 0, show.mutation = 0, scale.ratio = 2)
# 'o' is a list that includes (among other things) the coordinates; we extract them and adjust the names:
#   xy <- setNames(o[c("xx", "yy")], c("x", "y"))
# We can plot the original network using these coordinates so that the mutations will be correctly displayed
# plot(net3, xy = xy, size = size, cex = 0.8, pie = ind.hap3, threshold = 0, scale.ratio = 2, legend = TRUE)
# The haplotype labels overlap quite a bit, so adding 'labels = FALSE' gives a nice result too (unless you want print these labels).
# Best,
# Emmanuel