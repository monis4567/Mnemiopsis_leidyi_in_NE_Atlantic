#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
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
inpfnm2 <- "dnb_pip02_df.csv"
inpfnm3 <- "df_clo03.csv"
# paste path and file name together
wd00_wd05_inpfnm2 <- paste(wd00_wd05,"/",inpfnm2,sep="")
wd00_wd05_inpfnm3 <- paste(wd00_wd05,"/",inpfnm3,sep="")
# read in the csv file
df_pip02 <- read.csv(wd00_wd05_inpfnm2)
df_clo03 <- read.csv(wd00_wd05_inpfnm3,sep=";")
#make the date frame a matrix and a DNAbin object again
dnb_pip02 <- ape::as.DNAbin(as.matrix(df_pip02))
rownames(df_pip02) <- df_pip02$X
dnb_pip02 <- df_pip02[,-1]
pip <- dnb_pip02
#Make the dnabin a genind object
gei_pip <- adegenet::DNAbin2genind(pip)
# convert genind object to genlight object with the dartR package
gl_pip <- dartR::gi2gl(gei_pip)
gl_pip@ind.names <- rownames(df_pip02)
#gl_pip@loc.names <- rownames(df_pip02)
# load the "StAMPP" to be able to use functions from the "StAMPP" package
library("StAMPP")
# convert genlight object in to a stamp object
sta_mnelei <- StAMPP::stamppConvert(gl_pip,type = "genlight")
# split string and get lists nested in a list
loca2 <- strsplit(as.character(sta_mnelei$sample), "_")
#get second object in nested list
smplno2 <- sapply(loca2, "[[", 2)
# replace other location names
smplno2[grep("WaddenSeaBussumHaupstr",smplno2)] <-  "GermanyBusum" 
smplno2[grep("Germany:Maasholm",smplno2)] <- "GermanyKielFjord" 
# replace pop names in stamp object - you will need the identifiers for the 
# populations to be able 
sta_mnelei$pop.names <- df_clo03$locality8[match(smplno2,df_clo03$locality6)]
# use 'StAMPP' package to calculate Fst and p-values on populations
Mnl.fst <- StAMPP::stamppFst(sta_mnelei, 100, 95, 1)
# make the nested tables data frames
df_Mnelei.fst <- as.data.frame(Mnl.fst$Fsts)
df_Mnelei.pval <- as.data.frame(Mnl.fst$Pvalues)
df_Mnelei.pval <- t(df_Mnelei.pval)


m3 <- as.matrix(df_Mnelei.pval)
m2 <- as.matrix(df_Mnelei.fst)
rwnmMfst <- rownames(df_Mnelei.fst)
clnmMfst <- colnames(df_Mnelei.fst)
# m2 <- matrix(1:20, 4, 5)
# m3 <- matrix(rep(letters[1:5],4),4,5)
m2a <- ifelse(lower.tri(m2),m2,NA)
m3a <- ifelse(upper.tri(m3),m3,NA)
m4 <- ifelse(is.na(m2a),m3a,m2a)
rownames(m4) <- rwnmMfst
colnames(m4) <- clnmMfst
df_Mnelei.fst <- as.data.frame(m4)

# use the other 'StAMPP' function 
Mnl.Gst <- StAMPP::stamppGmatrix(sta_mnelei)
Mnl.NesD <- StAMPP::stamppNeisD(sta_mnelei)
# turn these tables in to data frames too
df_Mnelei.Gst <- as.data.frame(Mnl.Gst)
df_Mnelei.NesD <- as.data.frame(Mnl.NesD)
#_______________________________________________________________________________
library(htmlTable)
d1 <- df_Mnelei.fst
#d1 <- df_Mnelei.NesD
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
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
#colourtxthtml <- rep("color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals
# make a table caption
capt_tbl02 <-        "Table 3. Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. Sequences of nDNA ITS1-2 obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A high FST index indicates no variation among the individuals sampled, a high FST index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low FST values, and dark blue for high FST values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.3f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.Fst <- mdlo2
#mdlo2
# show the table
t.HTML03 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML03
#_______________________________________________________________________________
library(htmlTable)
d1 <- df_Mnelei.fst
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
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm<0.4,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
#colourtxthtml <- rep("color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals
# make a table caption
capt_tbl02 <-        "Table 4. Comparison of NeisD for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. Sequences of nDNA ITS1-2 obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A high NeisD index indicates no variation among the individuals sampled, a high NeisD index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low NeisD values, and dark blue for high NeisD values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.2f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.NeD <- mdlo2
#mdlo2
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
df_piplo2 <- pw_Nei_pip_lo
df_piplo2 <- as.matrix(pw_wc_pip_lo)
df_piplo2 <- as.matrix(pw_Gst_pip_lo)
# replace other location names
colnames(df_piplo2)[grep("WaddenSeaBussumHaupstr",colnames(df_piplo2))] <-  "GermanyBusum" 
rownames(df_piplo2)[grep("WaddenSeaBussumHaupstr",rownames(df_piplo2))] <-  "GermanyBusum" 
colnames(df_piplo2)[grep("Germany:Maasholm",colnames(df_piplo2))] <-  "GermanyKielFjord" 
rownames(df_piplo2)[grep("Germany:Maasholm",rownames(df_piplo2))] <-  "GermanyKielFjord" 
# replace long names with abbreviated names
rownames(df_piplo2) <- df_clo03$locality8[match(rownames(df_piplo2),df_clo03$locality6)]
colnames(df_piplo2) <- df_clo03$locality8[match(colnames(df_piplo2),df_clo03$locality6)]
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
# get lower triangle
df_df_piplo2 <- get_lower_tri(df_piplo2)
# get library for making 
library(htmlTable)
d1 <- df_df_piplo2 
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
# replace NA text colors with black
#coltxtval[grepl("NA",coltxtval)] <- gsub("NA","black",coltxtval[grepl("NA",coltxtval)])
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
# replace NA cell colors with white
#colourhtml_cell[grepl("NA",colourhtml_cell)] <- gsub("NA","white",colourhtml_cell[grepl("NA",colourhtml_cell)])
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
#colourhtml[grepl("NA",colourhtml)] <-  gsub("NA","transparent",colourhtml[grepl("NA",colourhtml)])
#colourtxthtml <- rep("color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
#colourtxthtml[grepl("NA",colourtxthtml)] <- gsub("NA","color:#000000;",colourtxthtml[grepl("NA",colourtxthtml)]) 
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals

colourhtml <- gsub("background-color:NA;","background-color:transparent;",colourhtml)
colourtxthtml <- gsub("color:NA;","color:#000000;",colourtxthtml)

#grepl("NA",colourhtml_cell)
# make a table caption
capt_tbl02 <-        "Table 5. Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. Sequences of nDNA ITS1-2 obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A high FST index indicates no variation among the individuals sampled, a high FST index indicates there is high variation among sequences compared. The color gradient reflects this as a yellow for low FST values, and dark blue for high FST values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.2f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)

mdlo2.NeiFst <- mdlo2
#mdlo2
# show the table
t.HTML05 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML05
