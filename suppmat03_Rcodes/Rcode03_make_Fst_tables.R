#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#remove everything in the working environment, without a warning!!
rm(list=ls())

if(!require(dartR)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  library(devtools)
  install_github("green-striped-gecko/dartR")
  gl.install.vanilla.dartR()
}

if(!require(StAMPP)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  library(devtools)
  install_github("lpembleton/StAMPP")
}
library(StAMPP)
#load required package libraries
library(dartR)
library(adegenet)
library(StAMPP)
library(ape)
#https://search.r-project.org/CRAN/refmans/StAMPP/html/stamppFst.html
#define path for input dir
wd01 <- "/suppmat01_inp_files"
wd05 <- "/suppmat05_out_files"
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
wd00_wd01 <- paste0(wd00,wd01)
wd00_wd05 <- paste0(wd00,wd05)
#Define input file name
inpfnm2 <- "algn_Mnelei_18s_17.fas.aligned.fasta"
inpfnm3 <- "df_clo03.csv"
# paste path and file name together
wd00_wd05_inpfnm2 <- paste(wd00_wd05,"/",inpfnm2,sep="")
wd00_wd05_inpfnm3 <- paste(wd00_wd05,"/",inpfnm3,sep="")
# read in the fasta file as a DNAbin object
dnb_02 <- read.FASTA(wd00_wd05_inpfnm2)
# read in the csv file
df_clo03 <- read.csv(wd00_wd05_inpfnm3,sep=";")
#copy the DNAbin object again
dnb_pip02 <- dnb_02
# check the label names
#labels(dnb_02)
pip <- dnb_pip02
#Make the dnabin a genind object
gei_pip <- adegenet::DNAbin2genind(pip)
# convert genind object to genlight object with the dartR package
gl_pip <- dartR::gi2gl(gei_pip)
# assign individual names
gl_pip@ind.names <- labels(dnb_02)
# and assign population names
locglp2<- strsplit(as.character(gl_pip$ind.names), "_")
gl_pip@pop <- as.factor(sapply(locglp2, "[[", 2))
# load the "StAMPP" to be able to use functions from the "StAMPP" package
library("StAMPP")
# convert genlight object in to a stamp object
sta_mnelei <- StAMPP::stamppConvert(gl_pip,type = "genlight")
# split string and get lists nested in a list
loca2 <- strsplit(as.character(sta_mnelei$sample), "_")
#get second object in nested list
smplno2 <- sapply(loca2, "[[", 2)
# get the unique population names
pNms <- unique(sta_mnelei$pop.names)
npNms<- length(pNms)
#bind together as columns and make it a data frame
df_popno <- as.data.frame(cbind(pNms,seq(1:npNms)))
colnames(df_popno) <- c("popAbb","popNo")
# add back a population number
sta_mnelei$pop.num <- df_popno$popNo[match(sta_mnelei$pop.names,df_popno$popAbb)]
# reorder the data frame based on sampling year
sta_mnelei <- sta_mnelei[order(sta_mnelei$pop.names),]
library(StAMPP)
## replace all NAs with zero
sta_mnelei[is.na(sta_mnelei)] <- 0
# import genotype data and convert to allele frequecies

# data(potato.mini, package="StAMPP")
# potato.freq <- stamppConvert(potato.mini, "r")
# # Calculate pairwise Fst values between each population
# potato.fst <- stamppFst(potato.freq, 100, 95, 1)

# use 'StAMPP' package to calculate Fst and p-values on populations
Mnl.fst <- StAMPP::stamppFst(sta_mnelei, 100, 95, 1)


# # import genotype data and convert to allele frequecies
# data(potato.mini, package="StAMPP")
# potato.freq <- stamppConvert(potato.mini, "r")
# # Calculate pairwise Fst values between each population
# potato.fst <- stamppFst(potato.freq, 100, 95, 1)

# pause before continuing
#Sys.sleep(12)
# make the nested tables data frames
df_Mnelei.fst <- as.data.frame(Mnl.fst$Fsts)
df_Mnelei.pval <- as.data.frame(Mnl.fst$Pvalues)
df_Mnelei.pval <- t(df_Mnelei.pval)
# make it a matrix
# m3 holds the p-values
m3 <- as.matrix(df_Mnelei.pval)
# and m2 holds the Fst values
m2 <- as.matrix(df_Mnelei.fst)
# get the row names and columns names
rwnmMfst <- rownames(df_Mnelei.fst)
clnmMfst <- colnames(df_Mnelei.fst)

m2a <- ifelse(lower.tri(m2),m2,NA)
m3a <- ifelse(upper.tri(m3),m3,NA)
# merge the lower and upper triangles-
# here m2a with the Fstvalues will appear in the lower triangle
# and m3a with the p-values  will appear in the upper triangle
m4 <- ifelse(is.na(m2a),m3a,m2a)
rownames(m4) <- rwnmMfst
colnames(m4) <- clnmMfst
df_Mnelei.fst <- as.data.frame(m4)
# copy the stampp object to use for G statistics
sta_mnelei.forGst <- sta_mnelei
# replace individual names with population names
sta_mnelei.forGst$sample <- sta_mnelei.forGst$pop.names 
# use the other 'StAMPP' function 
Mnl.Gst <- StAMPP::stamppGmatrix(sta_mnelei.forGst)
Mnl.NesD <- StAMPP::stamppNeisD(sta_mnelei)
# turn these tables in to data frames too
df_Mnelei.Gst <- as.data.frame(Mnl.Gst)
df_Mnelei.NesD <- as.data.frame(Mnl.NesD)
# change column names
colnames(df_Mnelei.NesD) <- row.names(df_Mnelei.NesD)
#_______________________________________________________________________________
library(htmlTable)
d1 <- df_Mnelei.fst
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
# but do not consider the values that are NA
minval <- min(mdlo[!is.na(mdlo)])
maxval <- max(mdlo[!is.na(mdlo)])
#minval <- 0
#maxval <- 1
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# make a table caption
capt_tbl02 <-        "Table 3. Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling locations. The lower triangle shows the Fst values, the corresponding probability values are in the upper triangle. Sequences of nDNA ITS1-2 obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A low FST index indicates no variation among the individuals sampled, a high FST index indicates there is high variation among sequences compared. The color gradient goes from yellow for low FST values, and dark blue for high FST values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.Fst <- mdlo2
# show the table
t.HTML03 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML03
#_______________________________________________________________________________
# split string and get lists nested in a list
loca3 <- strsplit(as.character(sta_mnelei$sample), "_")
#get second object in nested list
smplno3 <- sapply(loca3, "[[", 3)


# get the unique population names
pNms <- unique(smplno3)
npNms<- length(pNms)
#bind together as columns and make it a data frame
df_popno <- as.data.frame(cbind(pNms,seq(1:npNms)))
colnames(df_popno) <- c("popAbb","popNo")
# add back a population number
sta_mnelei$pop.num <- df_popno$popNo[match(smplno3,df_popno$popAbb)]
# rename locations with years
sta_mnelei$pop.names <- smplno3
# reorder the data frame based on sampling year
sta_mnelei <- sta_mnelei[order(sta_mnelei$pop.names),]
# use 'StAMPP' package to calculate Fst and p-values on populations
Mnl.fst <- StAMPP::stamppFst(sta_mnelei, 100, 95, 1)
# pause before continuing
#Sys.sleep(12)
# make the nested tables data frames
df_Mnelei.fst <- as.data.frame(Mnl.fst$Fsts)

df_Mnelei.pval <- as.data.frame(Mnl.fst$Pvalues)
df_Mnelei.pval <- t(df_Mnelei.pval)
# make it a matrix
# m3 holds the p-values
m3 <- as.matrix(df_Mnelei.pval)
# and m2 holds the Fst values
m2 <- as.matrix(df_Mnelei.fst)
# get the row names and columns names
rwnmMfst <- rownames(df_Mnelei.fst)
clnmMfst <- colnames(df_Mnelei.fst)

m2a <- ifelse(lower.tri(m2),m2,NA)
m3a <- ifelse(upper.tri(m3),m3,NA)
# merge the lower and upper triangles-
# here m2a with the Fstvalues will appear in the lower triangle
# and m3a with the p-values  will appear in the upper triangle
m4 <- ifelse(is.na(m2a),m3a,m2a)
rownames(m4) <- rwnmMfst
colnames(m4) <- clnmMfst
df_Mnelei.fst <- as.data.frame(m4)
#_______________________________________________________________________________
library(htmlTable)
d1 <- df_Mnelei.fst
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
# but do not consider the values that are NA
minval <- min(mdlo[!is.na(mdlo)])
maxval <- max(mdlo[!is.na(mdlo)])
#minval <- 0
#maxval <- 1
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# make a table caption
capt_tbl02 <-        "Table 6. Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. The the lower triangle shows the Fst values, the corresponding probability values are in the upper triangle. Sequences of nDNA ITS1-2 obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A high FST index indicates no variation among the individuals sampled, a high FST index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low FST values, and dark blue for high FST values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.Fst <- mdlo2
# show the table
t.HTML06 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML06




#_______________________________________________________________________________
library(htmlTable)
d1 <- df_Mnelei.NesD
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# change column names - in the NeisD table the column names are replaced
colnames(mdlo) <- rownames(mdlo)
#_______________________________________________________________________________
# Get lower triangles from tables
#_______________________________________________________________________________

#http://www.sthda.com/french/wiki/ggplot2-heatmap-d-une-matrice-de-corr-lation-logiciel-r-et-visualisation-de-donn-es
#https://stat.ethz.ch/R-manual/R-devel/library/base/html/lower.tri.html
# Obtenir le triangle inférieur
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Obtenir le triangle supérieur
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
# use the funciton to get the low triangle of the table
lwtr_mdlo <- get_lower_tri(mdlo)
# lower.tri(mdlo)
# upper.tri(mdlo)
#replace the table to use later on
mdlo <- lwtr_mdlo
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
# but do not consider the values that are NA
minval <- min(mdlo[!is.na(mdlo)])
maxval <- max(mdlo[!is.na(mdlo)])
#minval <- 0
#maxval <- 1
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# add white color for all NAs
colourhtmlvalues[is.na(colourhtmlvalues)] <- "white"
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# paste together background color og font color
colourhtml_cell <- paste0(colourhtmlvalues,coltxtval)
colourhtml_cell <- gsub("NA","white",colourhtml_cell)
# convert back to matrix
# colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtml_cell <- t(matrix(data=colourhtml_cell,nrow=nr,ncol=nc))
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
# make a table caption
capt_tbl02 <-        "Table 4. Comparison of NeisD for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. Sequences of nDNA ITS1-2 obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A high NeisD index indicates no variation among the individuals sampled, a high NeisD index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low NeisD values, and dark blue for high NeisD values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.NeD <- mdlo2
# replace NAs with blank
mdlo2 <- gsub("NA","",mdlo2)
# show the table
t.HTML04 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML04
#_______________________________________________________________________________
#_______________________________________________________________________________
#copy the data frame
# make a hierfstat object from the genind object
# use the previously made genlight object population names
hfst_pip_lo <- hierfstat::genind2hierfstat(gei_pip,pop=gl_pip@pop)
hfst_pip_lo <- hfst_pip_lo %>% replace(is.na(.), 0)
hfst_pip_lo <- hfst_pip_lo[hfst_pip_lo$pop %in% hfst_pip_lo$pop[duplicated(hfst_pip_lo$pop)],]
hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
pw_Nei_pip_lo <- hierfstat::pairwise.neifst(hfst_pip_lo, diploid = FALSE)
pw_Nei_pip_lo <- pw_Nei_pip_lo %>% replace(is.na(.), 0)
df_pw_pip_lo <- as.data.frame(pw_Nei_pip_lo)
df_piplo2 <- df_pw_pip_lo
df_piplo2 <- pw_Nei_pip_lo
mdlo <- df_piplo2
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# use the funciton to get the low triangle of the table
lwtr_mdlo <- get_lower_tri(mdlo)
#replace the table to use later on
mdlo <- lwtr_mdlo
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
# but do not consider the values that are NA
minval <- min(mdlo[!is.na(mdlo)])
maxval <- max(mdlo[!is.na(mdlo)])
#minval <- 0
#maxval <- 1
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
# make a table caption
capt_tbl02 <-        "Table 5. Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. Sequences of nDNA ITS1-2 obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A high FST index indicates no variation among the individuals sampled, a high FST index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low FST values, and dark blue for high FST values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.NeD <- mdlo2
# replace NAs with blank
mdlo2 <- gsub("NA","",mdlo2)
# show the table
t.HTML05 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML05


#_______________________________________________________________________________
# start - make Fst tables with only Danish and German samples
#_______________________________________________________________________________

# make the dnabin object a phydata object
phyd_pip6 <- phangorn::as.phyDat(dnb_pip02)
# strsplit names and turn into characters
# and rowbind the nested lists to a dataframe
df_pp6Nm <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(names(phyd_pip6)),
                                  "_")))
# modify columns names
colnames(df_pp6Nm) <- c("smplNm","loc","yer")
# get unique location names
uloc <- unique(df_pp6Nm$loc)
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
Nmsloc <- names(phyd_pip6)
# use grepl instead of dplyr, to get all NEA full sample names
# but first collapse and paste with '|' sign to be able to grep for
# multiple chracteristics
# https://stackoverflow.com/questions/67377806/dplyr-select-and-starts-with-on-multiple-values-in-a-variable-list?noredirect=1&lq=1
NmslocNEA <- Nmsloc[grepl(paste(NEAloc, collapse="|"),Nmsloc)]
# now subset the phydat object by the full names that match bein NEA samples
phyd_pip7 <- subset(phyd_pip6, NmslocNEA)
# limit to only include samples collected in 2017 and 2018
Nmswy17_18 <- NmslocNEA[grepl(paste(c("2017","2018"), collapse="|"),NmslocNEA)]
# subset the phy dat object - I am replacing the previous object
phyd_pip7 <- subset(phyd_pip6, Nmswy17_18)
# make the phy-dat object a DNAbin object instead
dnb_pip7 <- as.DNAbin(phyd_pip7)
# make the DNAbin object a distance matrix 
dst_pip7 <- ape::dist.dna(dnb_pip7, model= "raw")
# make it a matrix
mtx_pip7 <- as.matrix(dst_pip7)
# make it a data frame
df_pip7 <- as.data.frame(mtx_pip7)
# split the row name string by a character 
lpip7 <- strsplit(as.character(row.names(mtx_pip7)), "_")
# get the 2 nd element per vector in the list - this holds the abbreviated 
# location name
grp.loc7 <- sapply(lpip7, "[[", 2)

#Make the dnabin a genind object
gei_pip7 <- adegenet::DNAbin2genind(dnb_pip7)
# convert genind object to genlight object with the dartR package
gl_pip7 <- dartR::gi2gl(gei_pip7)
# assign individual names
gl_pip7@ind.names <- labels(dnb_pip7)
# and assign population names
locglp7.2<- strsplit(as.character(gl_pip7$ind.names), "_")
gl_pip7@pop <- as.factor(sapply(locglp7.2, "[[", 2))
# load the "StAMPP" to be able to use functions from the "StAMPP" package
library("StAMPP")
# convert genlight object in to a stamp object
sta_mnelei.7 <- StAMPP::stamppConvert(gl_pip7,type = "genlight")


# split string and get lists nested in a list
loca7.2 <- strsplit(as.character(sta_mnelei.7$sample), "_")
#get second object in nested list
smplno7.2 <- sapply(loca7.2, "[[", 2)
# get the unique population names
pNms7 <- unique(sta_mnelei.7$pop.names)
npNms7<- length(pNms7)
#bind together as columns and make it a data frame
df_popno7 <- as.data.frame(cbind(pNms7,seq(1:npNms7)))
colnames(df_popno7) <- c("popAbb","popNo")
# add back a population number
sta_mnelei.7$pop.num <- df_popno7$popNo[match(sta_mnelei.7$pop.names,df_popno7$popAbb)]
# reorder the data frame based on sampling year
sta_mnelei.7 <- sta_mnelei.7[order(sta_mnelei.7$pop.names),]
# use 'StAMPP' package to calculate Fst and p-values on populations
Mnl.fst7 <- StAMPP::stamppFst(sta_mnelei.7, 100, 95, 1)
# pause before continuing
#Sys.sleep(12)
# make the nested tables data frames
df_Mnelei.fst7 <- as.data.frame(Mnl.fst7$Fsts)
df_Mnelei.pval7 <- as.data.frame(Mnl.fst7$Pvalues)
df_Mnelei.pval7 <- t(df_Mnelei.pval7)
# make it a matrix
# m3 holds the p-values
m3 <- as.matrix(df_Mnelei.pval7)
# and m2 holds the Fst values
m2 <- as.matrix(df_Mnelei.fst7)
# get the row names and columns names
rwnmMfst <- rownames(df_Mnelei.fst7)
clnmMfst <- colnames(df_Mnelei.fst7)

m2a <- ifelse(lower.tri(m2),m2,NA)
m3a <- ifelse(upper.tri(m3),m3,NA)
# merge the lower and upper triangles-
# here m2a with the Fstvalues will appear in the lower triangle
# and m3a with the p-values  will appear in the upper triangle
m4 <- ifelse(is.na(m2a),m3a,m2a)
rownames(m4) <- rwnmMfst
colnames(m4) <- clnmMfst
df_Mnelei.fst7 <- as.data.frame(m4)
# copy the stampp object to use for G statistics
sta_mnelei.forGst <- sta_mnelei.7
# replace individual names with population names
sta_mnelei.forGst$sample <- sta_mnelei.forGst$pop.names 
# use the other 'StAMPP' function 
Mnl.Gst <- StAMPP::stamppGmatrix(sta_mnelei.forGst)
Mnl.NesD <- StAMPP::stamppNeisD(sta_mnelei.7)
# turn these tables in to data frames too
df_Mnelei.Gst <- as.data.frame(Mnl.Gst)
df_Mnelei.NesD <- as.data.frame(Mnl.NesD)
# change column names
colnames(df_Mnelei.NesD) <- row.names(df_Mnelei.NesD)

library(htmlTable)
d1 <- df_Mnelei.fst7
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
# but do not consider the values that are NA
minval <- min(mdlo[!is.na(mdlo)])
maxval <- max(mdlo[!is.na(mdlo)])
#minval <- 0
#maxval <- 1
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# make a table caption
capt_tbl02 <-        "Table 3. vers 02 Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the Danish and German sampling locations collected in 2017-2018. The lower triangle shows the Fst values, the corresponding probability values are in the upper triangle. A low FST index indicates no variation among the individuals sampled, a high FST index that apporaches 1 indicates there is high variation among sequences compared. The color gradient goes from yellow for low FST values, and dark blue for high FST values. NA = data not available."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.Fst <- mdlo2
# show the table
t.HTML03 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML03
#____________________




# split string and get lists nested in a list
loca3 <- strsplit(as.character(sta_mnelei.7$sample), "_")
#get second object in nested list
smplno3 <- sapply(loca3, "[[", 3)
# get the unique population names
pNms <- unique(smplno3)
npNms<- length(pNms)
#bind together as columns and make it a data frame
df_popno <- as.data.frame(cbind(pNms,seq(1:npNms)))
colnames(df_popno) <- c("popAbb","popNo")
# add back a population number
sta_mnelei.7$pop.num <- df_popno$popNo[match(smplno3,df_popno$popAbb)]
# rename locations with years
sta_mnelei.7$pop.names <- smplno3
# reorder the data frame based on sampling year
sta_mnelei.7 <- sta_mnelei.7[order(sta_mnelei.7$pop.names),]
# use 'StAMPP' package to calculate Fst and p-values on populations
Mnl.fst <- StAMPP::stamppFst(sta_mnelei.7, 100, 95, 1)
# pause before continuing
#Sys.sleep(12)
# make the nested tables data frames
df_Mnelei.fst <- as.data.frame(Mnl.fst$Fsts)

df_Mnelei.pval <- as.data.frame(Mnl.fst$Pvalues)
df_Mnelei.pval <- t(df_Mnelei.pval)
# make it a matrix
# m3 holds the p-values
m3 <- as.matrix(df_Mnelei.pval)
# and m2 holds the Fst values
m2 <- as.matrix(df_Mnelei.fst)
# get the row names and columns names
rwnmMfst <- rownames(df_Mnelei.fst)
clnmMfst <- colnames(df_Mnelei.fst)

m2a <- ifelse(lower.tri(m2),m2,NA)
m3a <- ifelse(upper.tri(m3),m3,NA)
# merge the lower and upper triangles-
# here m2a with the Fstvalues will appear in the lower triangle
# and m3a with the p-values  will appear in the upper triangle
m4 <- ifelse(is.na(m2a),m3a,m2a)
rownames(m4) <- rwnmMfst
colnames(m4) <- clnmMfst
df_Mnelei.fst <- as.data.frame(m4)
#_______________________________________________________________________________
library(htmlTable)
d1 <- df_Mnelei.fst
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
# but do not consider the values that are NA
minval <- min(mdlo[!is.na(mdlo)])
maxval <- max(mdlo[!is.na(mdlo)])
#minval <- 0
#maxval <- 1
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# make a table caption
capt_tbl02 <-        "Table 6. vers02  Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. The the lower triangle shows the Fst values, the corresponding probability values are in the upper triangle. A high FST index indicates no variation among the individuals sampled, a high FST index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low FST values, and dark blue for high FST values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.Fst <- mdlo2
# show the table
t.HTML06 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML06








#_______________________________________________________________________________
library(htmlTable)
d1 <- df_Mnelei.NesD
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# change column names - in the NeisD table the column names are replaced
colnames(mdlo) <- rownames(mdlo)
#_______________________________________________________________________________
# Get lower triangles from tables
#_______________________________________________________________________________

#http://www.sthda.com/french/wiki/ggplot2-heatmap-d-une-matrice-de-corr-lation-logiciel-r-et-visualisation-de-donn-es
#https://stat.ethz.ch/R-manual/R-devel/library/base/html/lower.tri.html
# Obtenir le triangle inférieur
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Obtenir le triangle supérieur
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
# use the funciton to get the low triangle of the table
lwtr_mdlo <- get_lower_tri(mdlo)
# lower.tri(mdlo)
# upper.tri(mdlo)
#replace the table to use later on
mdlo <- lwtr_mdlo
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
# but do not consider the values that are NA
minval <- min(mdlo[!is.na(mdlo)])
maxval <- max(mdlo[!is.na(mdlo)])
#minval <- 0
#maxval <- 1
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# add white color for all NAs
colourhtmlvalues[is.na(colourhtmlvalues)] <- "white"
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# paste together background color og font color
colourhtml_cell <- paste0(colourhtmlvalues,coltxtval)
colourhtml_cell <- gsub("NA","white",colourhtml_cell)
# convert back to matrix
# colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtml_cell <- t(matrix(data=colourhtml_cell,nrow=nr,ncol=nc))
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
# make a table caption
capt_tbl02 <-        "Table 4. vers02 Comparison of NeisD for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. A high NeisD index indicates no variation among the individuals sampled, a high NeisD index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low NeisD values, and dark blue for high NeisD values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.NeD <- mdlo2
# replace NAs with blank
mdlo2 <- gsub("NA","",mdlo2)
# show the table
t.HTML04 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML04
#_______________________________________________________________________________
#_______________________________________________________________________________
#copy the data frame
# make a hierfstat object from the genind object
# use the previously made genlight object population names
hfst_pip_lo <- hierfstat::genind2hierfstat(gei_pip7,pop=gl_pip7@pop)
hfst_pip_lo <- hfst_pip_lo %>% replace(is.na(.), 0)
hfst_pip_lo <- hfst_pip_lo[hfst_pip_lo$pop %in% hfst_pip_lo$pop[duplicated(hfst_pip_lo$pop)],]
hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
pw_Nei_pip_lo <- hierfstat::pairwise.neifst(hfst_pip_lo, diploid = FALSE)
pw_Nei_pip_lo <- pw_Nei_pip_lo %>% replace(is.na(.), 0)
df_pw_pip_lo <- as.data.frame(pw_Nei_pip_lo)
df_piplo2 <- df_pw_pip_lo
df_piplo2 <- pw_Nei_pip_lo
mdlo <- df_piplo2
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# use the funciton to get the low triangle of the table
lwtr_mdlo <- get_lower_tri(mdlo)
#replace the table to use later on
mdlo <- lwtr_mdlo
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
# but do not consider the values that are NA
minval <- min(mdlo[!is.na(mdlo)])
maxval <- max(mdlo[!is.na(mdlo)])
#minval <- 0
#maxval <- 1
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
# make a table caption
capt_tbl02 <-        "Table 5. vers02  Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years.  A high FST index indicates no variation among the individuals sampled, a high FST index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low FST values, and dark blue for high FST values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.NeD <- mdlo2
# replace NAs with blank
mdlo2 <- gsub("NA","",mdlo2)
# show the table
t.HTML05 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML05


#_______________________________________________________________________________
# end - make Fst tables with only Danish and German samples
#_______________________________________________________________________________

