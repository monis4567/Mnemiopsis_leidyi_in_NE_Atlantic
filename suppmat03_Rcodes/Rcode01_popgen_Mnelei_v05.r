#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# This code is able  to run in:

#remove everything in the working environment, without a warning!!
rm(list=ls())

# inspired by these webpages:
#https://johnbhorne.wordpress.com/2016/09/15/still-making-haplotype-networks-the-old-way-how-to-do-it-in-r/
#https://arundurvasula.wordpress.com/2016/02/24/haplotype-networks-in-r/
#https://github.com/emmanuelparadis/pegas/issues/4
#https://cran.r-project.org/web/packages/haplotyper/vignettes/haplotyper_vignette.html
#https://popgen.nescent.org/PopDiffSequenceData.html#pairwise-population-differentiation
#https://cran.r-project.org/web/packages/poppr/vignettes/mlg.html
#https://rdrr.io/cran/poppr/f/vignettes/poppr_manual.Rmd
#https://grunwaldlab.github.io/poppr/articles/poppr_manual.html
#http://darwin.eeb.uconn.edu/uncommon-ground/eeb348/

# See how to get all possible sequences from IUPAC notation using R :
# see:https://stackoverflow.com/questions/66272593/how-to-get-all-possible-sequences-from-iupac-notation-using-r
##############################################################################

#read in packages and install libraries
#install.packages("ape")
if(!require(ape)){
  install.packages("ape")
}
#read in the 'ape' library
library(ape)
#read in the 'pegas' library

if(!require(pegas)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  # install from GitHub
  devtools::install_github("emmanuelparadis/pegas/pegas")
}
library(pegas)

library(devtools)
if(!require(bioseq)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  # install from GitHub
  remotes::install_github("fkeck/bioseq")
}
library(bioseq)
# 
if(!require(ggtree)){
  BiocManager::install("ggtree")
 
}
library(ggtree)
library(ggplot2)
# load additional packages
library(tidyr)
library(dplyr)


 # I installed the Bioconductor package  DECIPHER
# some of the dependencies that came with the DECIPHER
# package appear to block out some of the base functions
# even though I tried to call the namespace for each for the 
# functions I could not call the base functions, and I 
# could not uninstall single Bioconductor packages.
# Instead I found this website
#https://support.bioconductor.org/p/7071/
# That suggested that I can remove everything
# so in my R library I decided to run 'rm -rf' in a terminal
# to get rid of everything
# Unfortunately, I need to re-install everything. 

#define path for input dir
wd01 <- "/suppmat01_inp_files"
wd05 <- "/suppmat05_out_files"
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
wd00_wd01 <- paste0(wd00,wd01)
wd00_wd05 <- paste0(wd00,wd05)
# delete a directory -- must add recursive = TRUE
unlink(wd00_wd05, recursive = TRUE)
#create anew directory
dir.create(wd00_wd05)
#set the working dir
setwd(wd00_wd01)
#check the working dir
getwd()


#define input file as variable
inpf01 <- "Mnelei_seq_2024_aug_06_v01.fasta"
inpf02 <- "individual_seq_regions_for_ITS_and_SrRNA_in_Mnemiopsis_leydi_v01.fasta"
inpf03 <- "output01_accn_publication_Mnemiopsis_its1its2.csv"
inpf04 <- "collect_loc_Mnelei_smpls01.csv"
# paste path and input file together
pth_inpf01 <- paste(wd00_wd01,"/",inpf01,sep="")
pth_inpf02 <- paste(wd00_wd01,"/",inpf02,sep="")
pth_inpf03 <- paste(wd00_wd01,"/",inpf03,sep="")
pth_inpf04 <- paste(wd00_wd01,"/",inpf04,sep="")
#read in csv file
df_lN02 <- read.csv(pth_inpf03,sep=";")
df_MneCl01 <- read.csv(pth_inpf04,sep=";")

#_______________________________________________________________________________
# section 01 -  start - substitute in names and input files
#_______________________________________________________________________________

# change incorrect location names
df_lN02$location[grepl("taxon",df_lN02$location)] <- NA
df_lN02$location[grepl("[0-9]+.*",df_lN02$location)] <- NA
# replace incorrect positions      
df_lN02$latlonpos[grepl("fwd_",df_lN02$latlonpos)] <- NA
df_lN02$latlonpos <- gsub("near ([A-Z]) ([0-9]+)(\\.)([0-9]+) ([A-Z]) ([0-9]+)(\\.)([0-9]+)","\\2\\3\\4 \\1 \\6\\7\\8 \\5",df_lN02$latlonpos)
df_lN02$location[grepl("Iran",df_lN02$latlonpos)] <- df_lN02$latlonpos[grepl("Iran",df_lN02$latlonpos)]
df_lN02$latlonpos[grepl("Iran",df_lN02$latlonpos)] <- NA
# get Iranian lat lon pos
latlonIran <- unique(df_lN02$latlonpos[grepl("Iran: Caspian",df_lN02$location)])
# assign lat lon pos for Iran
df_lN02$latlonpos[grepl("Iran$",df_lN02$location)] <-  latlonIran
#unique(df_lN02$pubjournal)
df_lN02$location[df_lN02$pubjournal=="Mol. Phylogenet. Evol. 21 (2), 218-230 (2001)"] <- "France Villefranche_sur_mer"
# I looked up papers published location names and obtained lat lon from google earth
# 43 39 N 17 20 E #"France Villefranche_sur_mer"
# make an empty NA column to add non decimal degrees to
df_lN02$latlonpos_nd <- NA
# add a non decimal degree to the French point 
df_lN02$latlonpos_nd[df_lN02$location=="France Villefranche_sur_mer"] <- "43 39 N 17 20 E"
df_lN02$location[df_lN02$location=="France Villefranche_sur_mer"] <- "Mediterranean"
# see https://www.sciencedirect.com/science/article/pii/S2352485519303366
# add non decimal coordinates for Colombian samples
df_lN02$latlonpos_nd[grepl("Multiple introductions",df_lN02$pubtitle)] <- "11 02 N 74 51 W"
df_lN02$location[grepl("Multiple introductions",df_lN02$pubtitle)] <- "Central W Atlantic"

#get more lat lon positions from the literature
#https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2010.04701.x
#Table_01 from REUSCH, T.B.H., BOLTE, S., SPARWEL, M., MOSS, A.G. and JAVIDPOUR, J. (2010), Microsatellites reveal origin and genetic diversity of Eurasian invasions by one of the world’s most notorious marine invader, Mnemiopsis leidyi (Ctenophora). Molecular Ecology, 19: 2690-2699. https://doi.org/10.1111/j.1365-294X.2010.04701.x
Rta01 <- "Sampling site 	Invasion/sampling year 	Abbreviation 	rDNA sequencing 	N 	Co-ordinates 	Collector (Institution)
Black Sea Bulgaria 	∼1985/2009 	BSB 	Yes 	15 	43°11′N, 27°57′ E 	Kremena Stefanova (Bulgarian Institute of Oceanology)
Black Sea Turkey 	∼1985/2009 	BST 		23 	42° 01′ N, 35°08′ E 	Levent Bat (Sinop University)
Black Sea Ukraine 	∼1985/2009 	BSU 		48 	44°37′N 33°31′E 	Alexandra Gubanova (Institute of Biology of the Southern Seas)
Caspian Sea 	1999/2008 	CAS 	Yes 	22 	36°48′N, 53° 07′E 	Jamileh Javidpour (IFM-GEOMAR)
Fehmarn Belt 	2006/2009 	FEM 		48 	54°30′N, 11°20′E 	Sören Bolte (IFM- GEOMAR)
Kristineberg 	2006/2009 	KRS 		37 	58°24′N, 11°24′E 	Lene Fries Möller (University of Gothenburg)
Helgoland 	2006/2009 	HEL 	Yes 	46 	54°11′N, 07°53′E 	Philipp Schubert (IFM- GEOMAR)
Kiel Fjord 	2006/2009 	KIF 	Yes 	21 	54°25′N, 10°12′E 	Sören Bolte (IFM- GEOMAR)
Bornholm Basin 	2006/2009 	BBA 	Yes 	48 	55°11′N, 15°32′E 	Bastian Huver (DTU Aqua)
Maasholm 	2006/2008 	MAH 	Yes 	42 	54˚41′N, 10˚00′E 	Maxi Sparwel (Uni Münster/IFM- GEOMAR)
USA Woods Hole 	Native/2009 	USA-WH 	Yes 	45 	41° 31′N, 70°40′W 	Woods Hole specimen service
USA Panacea 	Native/2009 	USA-PC 	Yes 	33 	30°00′N, 84°20′W 	Jack Rudloe (Gulf Specimen Marine Lab)
USA Port Aransas 	Native/2009 	USA-PA 		19 	27°50′N, 97°03′W 	Anthony G. Moss (Auburn University)
USA Galveston Bay 	Native/2009 	USA-GB 	Yes 	20 	29°17′N, 94°52′W 	Anthony G. Moss (Auburn University)"

Rta01 <- gsub("BST 		","BST 		No 		",Rta01)
Rta01 <- gsub("BSU 		","BSU 		No 		",Rta01)
Rta01 <- gsub("FEM 		","FEM 		No 		",Rta01)
Rta01 <- gsub("KRS 		","KRS 		No 		",Rta01)
Rta01 <- gsub("USA-PA 		","USA-PA 		No 		",Rta01)
# replace space in text string
Rta01 <- gsub(" ","_", Rta01)
#Table_02 from REUSCH, T.B.H., BOLTE, S., SPARWEL, M., MOSS, A.G. and JAVIDPOUR, J. (2010), Microsatellites reveal origin and genetic diversity of Eurasian invasions by one of the world’s most notorious marine invader, Mnemiopsis leidyi (Ctenophora). Molecular Ecology, 19: 2690-2699. https://doi.org/10.1111/j.1365-294X.2010.04701.x
Rta02 <-  "Location, sample no. 	Accession number 	ITS1 	ITS2 	ITS2
285 	489 	508
USA—Woods Hole† 	AF293700 	C 	A 	T
The Netherlands‡ 	EF175463 	C 	– 	–
USA—Woods Hole_01 	HM147257 	C 	T 	C
USA—Woods Hole_16 	HM147272 	C 	– 	–
USA—Panacea_22 	HM147273 	T 	– 	–
USA—Panacea_03 	HM147274 	C 	– 	–
USA—Galveston Bay02 	HM147275 	– 	T 	–
Bornhom Bassin_01 	HM147258 	C 	T 	C
Bornholm Bassin_11 	HM147259 	C 	T 	C
Bornholm Bassin_12 	HM147260 	C 	T 	C
Bornholm Bassin_13 	HM147261 	C 	T 	C
Helgoland_052 	HM147262 	C 	T 	C
Helgoland_122 	HM147263 	– 	T 	T
Helgoland _82 	HM147264"

#  split string by end-of-line character
r <- unlist(strsplit(Rta01, "\n"))
# read in text string as a data frame 
#https://stackoverflow.com/questions/41076240/r-convert-large-character-string-to-dataframe
Rta01.1 <- read.table(text=gsub("(?<=[a-z])\\s+", "\t", r, perl=TRUE), 
           header=F, fill=T)
# replace column names
colnames(Rta01.1) <- Rta01.1[1,]
# remove line 1 and overwrite data frame
Rta01.1 <- Rta01.1[-1,]
# replace in column
Rta01.1$Sampling_site_ <- gsub("_$","",Rta01.1$Sampling_site_)
Rta01.1$`Co-ordinates_` <- gsub("°","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub("˚","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub("′","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub(",_","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub("__","_",Rta01.1$`Co-ordinates_`)
#https://stackoverflow.com/questions/41076240/r-convert-large-character-string-to-dataframe
Rta01.1 <- read.table(text=gsub("(?<=[a-z])\\s+", "\t", r, perl=TRUE), 
           header=F, fill=T)
# replace column names
colnames(Rta01.1) <- Rta01.1[1,]
# remove line 1 and overwrite data frame
Rta01.1 <- Rta01.1[-1,]
# replace in column
Rta01.1$Sampling_site_ <- gsub("_$","",Rta01.1$Sampling_site_)
Rta01.1$`Co-ordinates_` <- gsub("°","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub("˚","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub("′","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub(",_","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub("__","_",Rta01.1$`Co-ordinates_`)
Rta01.1$`Co-ordinates_` <- gsub("_"," ",Rta01.1$`Co-ordinates_`)
#replace space with underscore
Rta02 <- gsub(" ","_", Rta02)
#  split string by end-of-line character
r <- unlist(strsplit(Rta02, "\n"))
# read in text string as a data frame 
#https://stackoverflow.com/questions/41076240/r-convert-large-character-string-to-dataframe
Rta02.1 <- read.table(text=gsub("(?<=[a-z])\\s+", "\t", r, perl=TRUE), 
                      header=F, fill=T)
# remove line 2 and overwrite data frame
Rta02.1 <- Rta02.1[-2,]
# replace column names
colnames(Rta02.1) <- Rta02.1[1,]
# remove line 1 and overwrite data frame
Rta02.1 <- Rta02.1[-1,]
# replace in columns
Rta02.1$`Location,_sample_no._` <- gsub("‡_","",Rta02.1$`Location,_sample_no._`)
Rta02.1$`Location,_sample_no._` <- gsub("†_","",Rta02.1$`Location,_sample_no._`)
Rta02.1$`Location,_sample_no._` <- gsub("_[0-9]+_$","",Rta02.1$`Location,_sample_no._`)
Rta02.1$`Location,_sample_no._` <- gsub("Bornhom_","Bornholm_",Rta02.1$`Location,_sample_no._`)
Rta02.1$`Location,_sample_no._` <- gsub("[0-9]+_$","",Rta02.1$`Location,_sample_no._`)
Rta02.1$`Location,_sample_no._` <- gsub("_$","",Rta02.1$`Location,_sample_no._`)
Rta02.1$`Location,_sample_no._` <- gsub("—","_",Rta02.1$`Location,_sample_no._`)
Rta02.1$Accession_number_ <- gsub("_$","",Rta02.1$Accession_number_)
# match coordinates back to data frame
Rta02.1$`Location,_sample_no._` <- gsub("Bornholm_Bassin","Bornholm_Basin",Rta02.1$`Location,_sample_no._`)
Rta02.1$coordinates <- Rta01.1$`Co-ordinates_`[match(Rta02.1$`Location,_sample_no._`,Rta01.1$Sampling_site_)]

# Netherland sample is from https://www.vliz.be/imisdocs/publications/ocrd/132874.pdf
#The identical M. leidyi sequences were deposited
#in the NCBI Genbank database under the
#accession number EF175463,
Rta02.1$coordinates[Rta02.1$`Location,_sample_no._`== "The_Netherlands"] <- "52 56 N 4 36 E"
#change column name of first column
colnames(df_lN02)[1] <- c("accession_nmb")
#  match back to get coordinates
df_lN02$latlonpos_nd2 <- Rta02.1$coordinates[match(df_lN02$accession_nmb,Rta02.1$Accession_number_)]
# get position from USA Woods Hole
llposUSAwh<- unique(df_lN02$latlonpos[grepl("USA: Woods Hole",df_lN02$location)])
#
# The paper by Podar et al., (2001): Molecular Phylogenetics and Evolution , Vol. 21, No. 2, November, pp. 218 –230, 2001, doi:10.1006/mpev.2001.1036, available online at http://www.idealibrary.com on
# https://reader.elsevier.com/reader/sd/pii/S105579030191036X?token=D40BF30E74231EEA2E9D6971E571596BE20B428EB9392477B765DEF660253109213556AACBD2382676EC59112C53FAD3&originRegion=eu-west-1&originCreation=20220504094633
# Lists sequence AF293700 as coming from "USA: Woods Hole" 
df_lN02$location[grepl("AF293700",df_lN02$accession_nmb)] <- "USA: Woods Hole" 
# assign the Woods Hole position for the USA Woods Hole sample 'AF293700'
df_lN02$latlonpos[grepl("AF293700",df_lN02$accession_nmb)] <- llposUSAwh
# overwrite the incorrect position for the USA Woods Hole sample 'AF293700'
df_lN02$latlonpos_nd2[grepl("AF293700",df_lN02$accession_nmb)] <- NA
df_lN02$latlonpos_nd[grepl("AF293700",df_lN02$accession_nmb)] <- NA
# use grep to get KF435 accession numbers
#see https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081067#pone-0081067-t001
df_lN02$location[grepl("^KF435",df_lN02$accession_nmb)] <- "Mediterranean Sea"
df_lN02$latlonpos_nd[grepl("^KF435",df_lN02$accession_nmb)] <- "43 25 N 9 17 E" #"43 39 N 17 20 E" 
# add location for publication title
df_lN02$location[grepl("'Majaneva,S.'",df_lN02$pubauthor)] <- "Norway"
df_lN02$latlonpos_nd[grepl("Norway",df_lN02$location)] <- "57 57 N 7 57 E"
# add location for accession number
df_lN02$location[grepl("^KJ75",df_lN02$accession_nmb)] <- "Mediterranean Sea"
df_lN02$latlonpos_nd[grepl("^KJ75",df_lN02$accession_nmb)] <- "43 25 N 9 17 E" #"43 39 N 17 20 E" 
# add location for publication title
df_lN02$location[grepl("Blooms of the invasive ctenophore",df_lN02$pubtitle)] <- "Mediterranean Sea"
df_lN02$latlonpos_nd[grepl("Blooms of the invasive ctenophore",df_lN02$pubtitle)] <-"43 25 N 9 17 E" #"43 39 N 17 20 E" 
# add location for publication title
df_lN02$location[grepl("First molecular species identification of Mnemiopsis",df_lN02$pubtitle)] <- "Mediterranean Sea"
df_lN02$latlonpos_nd[grepl("First molecular species identification of Mnemiopsis",df_lN02$pubtitle)] <-"43 25 N 9 17 E" #"43 39 N 17 20 E" 
# add location for accession number
df_lN02$location[grepl("^JQ742005",df_lN02$accession_nmb)] <- "Caspian Sea" #"43 25 N 9 17 E"
# get unique location for Caspian sea
Caspsealoc <- unique(df_lN02$latlonpos[grepl("Iran: Caspian Sea",df_lN02$location)])
df_lN02$latlonpos[grepl("^JQ742005",df_lN02$accession_nmb)] <- Caspsealoc #"43 25 N 9 17 E"
# Assign location for the Netherlands sample
df_lN02$location[grepl("^EF175463",df_lN02$accession_nmb)] <- "Netherlands"

#replace in location names
df_lN02$location[grepl("Bornholm",df_lN02$location)] <- "Baltic Sea"
df_lN02$location[grepl("Denmark",df_lN02$location)] <- "NE Atlantic"
df_lN02$location[grepl("Iran",df_lN02$location)] <- "Caspian Sea"
df_lN02$location[grepl("Caspian",df_lN02$location)] <- "Caspian Sea"
df_lN02$location[grepl("Norway",df_lN02$location)] <- "NE Atlantic"
df_lN02$location[grepl("Mediterranean",df_lN02$location)] <- "Mediterranean"
# replace lat and lon for the NA entries
df_lN02$latlonpos_nd[is.na(df_lN02$latlonpos_nd)] <- df_lN02$latlonpos_nd2[is.na(df_lN02$latlonpos_nd)]
# identify the coordinates with decimal places and replace with spaces
pLaLo <- df_lN02$latlonpos[is.na(df_lN02$latlonpos_nd)]
pLaLo <- gsub('\\.| ',' ',pLaLo)
df_lN02$latlonpos[is.na(df_lN02$latlonpos_nd)] <- pLaLo
# get the mssing coordinates for the coordinates missing
df_lN02$latlonpos_nd[is.na(df_lN02$latlonpos_nd)] <- df_lN02$latlonpos[is.na(df_lN02$latlonpos_nd)]
# substitute extra spaces at the end
df_lN02$latlonpos_nd <- gsub(" $","",df_lN02$latlonpos_nd)
# copy the column 
df_lN02$latlonpos_nd2 <- df_lN02$latlonpos_nd
# split a column by a delimiter
df_lN02 <- 
  tidyr::separate(df_lN02, latlonpos_nd, 
                           sep = " ", into = paste0("ll", 1:6), 
                           fill = "right")
# paste together to get DMS coordinates
df_lN02$lat01<- paste0( df_lN02$ll1,"º",
                        df_lN02$ll2,"'",
                        df_lN02$ll3)
df_lN02$lon01<- paste0( df_lN02$ll4,"º",
                        df_lN02$ll5,"'",
                        df_lN02$ll6)

#_______________________________________________________________________________
# function : 'dms2dec' - start
#_______________________________________________________________________________
# get function from this website
# https://www.r-bloggers.com/2022/02/degree-minute-second-to-decimal-coordinates/
dms2dec <- function(dms, separators = c("º", "°", "\'", "’", "’’", "\"", "\'\'", "\\?")) {
  
  # version 1.4 (2 Feb 2022)
  # dms: a vector of latitude or longitude in degrees-minutes-seconds-hemisfere, e.g. 41° 34' 10.956" N (with or without spaces)
  # separators: the characters that are separating degrees, minutes and seconds in 'dms'; mind these are taken in the order in which they appear and not interpreted individually, i.e. 7'3º will be taken as 7 degrees, 3 minutes! input data are assumed to be properly formatted
  
  dms <- as.character(dms)
  dms <- gsub(pattern = " ", replacement = "", x = dms)
  for (s in separators) dms <- gsub(pattern = s, replacement = "_splitHere_", x = dms)
  
  splits <- strsplit(dms, split = "_splitHere_")
  n <- length(dms)
  deg <- min <- sec <- hem <- vector("character", n)
  
  for (i in 1:n) {
    deg[i] <- splits[[i]][1]
    min[i] <- splits[[i]][2]
    
    if (length(splits[[i]]) < 4) {
      hem[i] <- splits[[i]][3]
    } else {
      sec[i] <- splits[[i]][3]
      hem[i] <- splits[[i]][4]
    }
  }
  
  dec <- colSums(rbind(as.numeric(deg), (as.numeric(min) / 60), (as.numeric(sec) / 3600)), na.rm = TRUE)
  sign <- ifelse (hem %in% c("N", "E"), 1, -1)
  hem_miss <- which(is.na(hem))
  if (length(hem_miss) > 0) {
    warning("Hemisphere not specified at position(s) ", hem_miss, ", so the sign of the resulting coordinates may be wrong.")
  }
  dec <- sign * dec
  return(dec)
}  
#_______________________________________________________________________________
# function : 'dms2dec' - end
#_______________________________________________________________________________
# Use the 'dms2dec' function on the DMS pasted coordinates
df_lN02$lat02 <- dms2dec(df_lN02$lat01)
df_lN02$lon02 <- dms2dec(df_lN02$lon01)
# make new columns with decimal degress
df_lN02$declat <- df_lN02$lat02
df_lN02$declon <- df_lN02$lon02
# replace any spaces in location names
df_lN02$location <- gsub(" ","",df_lN02$location)
# add a column for sample year
df_lN02$smplyear <- gsub(".*(20[0-9]{2}).*","\\1",df_lN02$pubjournal)
df_lN02$smplyear[grepl("Unpublished",df_lN02$smplyear)] <- "unknown"
# define a filename to write the data frame to as a csv file
wd00_wd05_flnm3 <- paste(wd00_wd05,"/df_lN02.csv",sep="")
# write the data frame as a csv file
write.csv(df_lN02,file=wd00_wd05_flnm3)

#define input file
fl1 <- "collect_loc_Mnelei_smpls01.csv"
#paste path and filename together
pth_fl01 <- paste(wd00_wd01,"/",fl1,sep="")
#read in a csv file with positions for sampling locations
df_clo <- as.data.frame(read.csv(pth_fl01,sep = ";"))
df_clo$locality4 <-  df_clo$locality
#Modify location names
df_clo$locality4 <- gsub("FynBogense","Funen, Bogense", df_clo$locality4 )
df_clo$locality4 <- gsub("FynKerteminde","Funen, Kerteminde", df_clo$locality4 )
df_clo$locality4 <- gsub("JyllandMariagerfjord","Jutland, Mariager Fjord", df_clo$locality4 )
df_clo$locality4 <- gsub("NGermanyKielFjord","Germany, Kiel Fjord", df_clo$locality4 )
df_clo$locality4 <- gsub("NGermanyMecklenburgerBuchtWismarBucht","Germany, Wismar Bight", df_clo$locality4 )
df_clo$locality4 <- gsub("NJyllandLimfjord","Jutland, Limfjord", df_clo$locality4 )
df_clo$locality4 <- gsub("NJyllandLimfjordLogstoer","Jutland, Limfjord" , df_clo$locality4 )
df_clo$locality4 <- gsub("NWGermanyNSeaHelgolandRds", "North Sea, Helgoland Roads", df_clo$locality4 )
df_clo$locality4 <- gsub("NWGermanyWaddenSeaBussumHaupstr","Germany, Büsum", df_clo$locality4 )
df_clo$locality4 <- gsub("SamsoeBallen","Samsøe, Ballen", df_clo$locality4 )
df_clo$locality4 <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", df_clo$locality4 )

# paste together to get DMS coordinates
df_clo$DMSlat<- paste0( df_clo$lat_deg,"º",
                       df_clo$lat_min,"'",
                       df_clo$lat_sec,"''",
                       df_clo$lat_sph)
df_clo$DMSlon<- paste0( df_clo$lon_deg,"º",
                         df_clo$lon_min,"'",
                         df_clo$lon_sec,"''",
                         df_clo$lon_sph)
# Use the 'dms2dec' function on the DMS pasted coordinates
df_clo$dec_lat <- dms2dec(df_clo$DMSlat)
df_clo$dec_lon <- dms2dec(df_clo$DMSlon)

#exclude the location: 'SEDenmarkMecklenburgerBucht'
# as it is the same as "NGermanyMecklenburgerBuchtWismarBucht", and only
#one of these are needed
# same for the location 'NJyllandLimfjordLogstoer' which is identical to NJyllandLimfjord
df_clo <- df_clo[!df_clo$locality=="SEDenmarkMecklenburgerBucht",]
df_clo <- df_clo[!df_clo$locality=="NJyllandLimfjordLogstoer",]


lN02loc <- unique(df_lN02$location)
lN02lat <- df_lN02$declat[match(lN02loc,df_lN02$location)]
lN02lon <- df_lN02$declon[match(lN02loc,df_lN02$location)]
df_locN02 <- as.data.frame(cbind(lN02loc,lN02lat,lN02lon))
# get number of locations for Genbank records
nlocgb <- nrow(df_locN02)
df_clo03 <- df_clo
nclloc <- ncol(df_clo)
df_clo03 [ nrow(df_clo03) + nlocgb , ] <- NA
df_clo03$locality4[is.na(df_clo03$locality4)] <- df_locN02$lN02loc
df_clo03$dec_lat[is.na(df_clo03$locality2)] <- df_locN02$lN02lat
df_clo03$dec_lon[is.na(df_clo03$locality2)] <- df_locN02$lN02lon


#replace non ascii characteres
df_clo03$locality5 <- textclean::replace_non_ascii(df_clo03$locality4)
# Replace "Germany:Maasholm" with "Germany, Kiel Fjord"
df_clo03$locality5[df_clo03$locality5=="Germany:Maasholm" ] <- "Germany, Kiel Fjord" 
df_clo03$locality4[df_clo03$locality4=="Germany:Maasholm" ] <- "Germany, Kiel Fjord" 
df_clo03$locality6 <- df_clo03$locality5
df_clo03$locality6[df_clo03$locality6=="Germany:Maasholm" ] <- "Germany, Kiel Fjord" 
#add a column with sampling locations to be able to 
# match the localities
df_clo03$locality2 <- abbreviate(df_clo03$locality5, 1, named = FALSE)
# gep for NCBI occurence in locality names
df_clo03 <- df_clo03[!grepl("NCBI",df_clo03$locality4),]
# substitute i locality names
df_clo03$locality6 <- gsub(" ","",df_clo03$locality5)
df_clo03$locality6 <- gsub(",","",df_clo03$locality6)
df_clo03$locality6 <- gsub("Funen","Fyn",df_clo03$locality6)
df_clo03$locality6 <- gsub("Jutland","Jylland",df_clo03$locality6)
df_clo03$locality6 <- gsub("Sealand","Sjaelland",df_clo03$locality6)
df_clo03$locality6 <- gsub("JyllandMariagerFjord","JyllandMariagerfjord",df_clo03$locality6)
df_clo03$locality6 <- gsub("Germany:KielFjord", "NGermanyKielFjord" ,df_clo03$locality6)
df_clo03$locality6 <- gsub("JyllandLimfjord","NJyllandLimfjord"  ,df_clo03$locality6)
df_clo03$locality6 <- gsub("GermanyWismarBight" ,"NGermanyMecklenburgerBuchtWismarBucht" ,df_clo03$locality6)
df_clo03$locality6 <- gsub("NorthSeaHelgolandRoads","NWGermanyNSeaHelgolandRds"  ,df_clo03$locality6)
#substitue abbreviations
df_clo03$locality2[df_clo03$locality2=="GWB"] <- "GWi"
df_clo03$locality2[df_clo03$locality2=="G:K"] <- "GKi"
df_clo03$locality2[df_clo03$locality2=="G:H"] <- "GHe"
df_clo03$locality2[df_clo03$locality2=="G:M"] <- "GKi"
df_clo03$locality2[df_clo03$locality2=="G:M"] <- "GKi"
df_clo03$locality2[df_clo03$locality2=="SS"] <- "SSk"
df_clo03$locality2[df_clo03$locality2=="SB"] <- "SBa"
df_clo03$locality2[df_clo03$locality2=="GB"] <- "GBu"
df_clo03$locality2[df_clo03$locality2=="GKF"] <- "GKi"
df_clo03$locality2[df_clo03$locality2=="JMF"] <- "JMa"
df_clo03$locality2[df_clo03$locality2=="FK"] <- "FKe"
df_clo03$locality2[df_clo03$locality2=="FB"] <- "FBo"
df_clo03$locality2[df_clo03$locality2=="JL"] <- "JLi"
df_clo03$locality2[df_clo03$locality2=="NSHR"] <- "GHe"
df_clo03$locality2[df_clo03$locality2=="CW"] <- "CWA"
df_clo03$locality2[df_clo03$locality2=="AO"] <- "NWA"
df_clo03$locality2[df_clo03$locality2=="NE"] <- "NEA"
df_clo03$locality2[df_clo03$locality2=="A"] <- "NWA"

#order the data fram with location abbreviations to get 
df_clo03 <- df_clo03[order(df_clo03$locality2),]
# make vector of long location names
llocnms <- c("Fyn Bogense",
             "Fyn Kerteminde",
             "Jylland Mariager fjord",
             "N Germany Kiel Fjord",
             "N Germany Mecklenburger Bucht Wismar Bucht",
             "N Jylland Limfjord",
             "N Jylland Limfjord Logstoer",
             "NW Germany N Sea Helgoland Rds",
             "NW Germany Wadden Sea Bussum Haupstr",
             "Samsoe Ballen",
             "SE Denmark Mecklenburger Bucht",
             "Sjaelland Skovshoved"
)
# make a vector of abbreviated location names
abbloccnms <- 
  c("FBog",
    "FKer",
    "JMar",
    "GKie",
    "GWis",
    "JLim",
    "JLog",
    "GHel",
    "GBus",
    "SBal",
    "GWie",
    "Ssko")
# if the vectors are of equal length
if (length(llocnms) == length(abbloccnms)){
#combine to a data frame
df_abblocnms <- as.data.frame(cbind(llocnms,abbloccnms))}

#substitue in data frame
df_abblocnms$locnm2 <- gsub(" ","",df_abblocnms$llocnms)
#copy data frame
df_clo03$locality7 <- df_clo03$locality2

# # round decimals in lat lon positions to and paste together
df_lN02$declat2 <- round(df_lN02$declat,3)
df_lN02$declon2 <- round(df_lN02$declon,3)
df_lN02$latlonpos3 <- paste(df_lN02$declat2,"; ",df_lN02$declon2,sep="")
# # round decimals in lat lon positions to and paste together
df_clo03$dec_lat2 <-  round(as.numeric(df_clo03$dec_lat),3)
df_clo03$dec_lon2 <-  round(as.numeric(df_clo03$dec_lon),3)
df_clo03$latlonpos3 <- paste(df_clo03$dec_lat2,"; ",df_clo03$dec_lon2,sep="")
# replace sampling year for publication year, or use sampling year
# as listed in NCBI GenBank records
df_lN02$smplyear2 <- df_lN02$smplyear
df_lN02$smplyear2[df_lN02$smplyear2=="unknown"]
df_lN02$smplyear2[grepl("JX26",df_lN02$accession_nmb)] <- "2009"
df_lN02$smplyear2[grepl("GU",df_lN02$accession_nmb)] <- "2009"
df_lN02$smplyear <- df_lN02$smplyear2
# reduce to only having 4 characters
df_lN02$smplyear2 <- gsub("unknown","unkn",df_lN02$smplyear2)
# write the df_lN02 data frame to a csv file
write.csv(df_lN02,paste0(wd00_wd05,"/df_lN02.csv"))
# substitute in locality names
df_clo03$locality4 <- gsub("Bogense,","Bogense",df_clo03$locality4)
# copy column in data frame
df_clo03$locality8 <- df_clo03$locality7
# grep for long location names and replace with 4 letter abbreviations
df_clo03$locality8[grepl("Bogense",df_clo03$locality4)] <- "FBog"
df_clo03$locality8[grepl("Kerteminde",df_clo03$locality4)] <- "FKer"
df_clo03$locality8[grepl("North Sea, Helgoland",df_clo03$locality4)] <- "GHel"
df_clo03$locality8[grepl("Kiel Fjord",df_clo03$locality4)] <- "GKie"
df_clo03$locality8[grepl("Limfjord",df_clo03$locality4)] <- "JLim"
df_clo03$locality8[grepl("Mariager",df_clo03$locality4)] <- "JMar"
df_clo03$locality8[grepl("Samsøe",df_clo03$locality4)] <- "SBal"
df_clo03$locality8[grepl("Skovshoved",df_clo03$locality4)] <- "SSko"

# subset the data frame by two criteria that both must be TRUE, To remove
# the redundant Kiel locations
df_clo04 <- df_clo03[ which( !grepl("Kiel",df_clo03$locality4) | !is.na(df_clo03$locality)) , ]
df_clo04 <- df_clo04[ which( !grepl("Helgoland",df_clo04$locality4) | !is.na(df_clo04$locality)) , ]
#make a vector that will match the location names used in the sequences
llongNm <- c("KielFjord", "Ballen", "NSeaHelgolandRds", 
          "MecklenburgerBuchtWismarBucht", 
          "Loegstoer", "Bogense", "MecklenburgerBucht", "Kerteminde", "Mariagerfjord", 
          "Skovshoved", "WaddenSeaBussumHaupstr",
          "BalticSea", "CaspianSea", "CentralWAtlantic", "Funen, Bogense", 
          "Funen, Kerteminde", "Germany, Büsum", "North Sea, Helgoland Roads", 
          "Germany, Kiel Fjord", "Germany, Wismar Bight", "Jutland, Limfjord", 
          "Jutland, Mariager Fjord", "Mediterranean", "NEAtlantic", "Netherlands", 
          "AtlanticOcean:NWAtlantic", "Samsøe, Ballen", "Sealand, Skovshoved", 
          "USA:GalvestonBay", "USA:Panacea", "USA:WoodsHole",
          
          "Germany:KielFjord","Germany:Maasholm","Germany:Helgoland")  

llabbr <- c("GKie", "SBal", "GHel", 
            "GWis", 
            "JLog", "FBog", "GWis", 
            "FKer", "JMar", 
            "SSko", "GBus", 
            "BalS", "CasS", "CWAt", "Fbog", 
            "FKer", "GBus", "GHel", 
            "GKie", "GWis", "JLim", 
            "JMar", "Medi", "NEAt", "Neth", 
            "NWAt", "SBal", "SSko", 
            "USAG", "USAP", "USAW",
            
            "GKie","GKie","GHel")
# combine to a dataframe
df_4LtabbrNm <- as.data.frame(cbind(llongNm,llabbr))
# match to get the 4 letter abbreviation added to the data frame
df_clo04$loc4Lett<- df_4LtabbrNm$llabbr[match(df_clo04$locality4,df_4LtabbrNm$llongNm)]
# use the four letter abbreviation data frame to get the 
# short names for the NCBI collected sequences
df_lN02$llabbr <- df_4LtabbrNm$llabbr[match(df_lN02$location,df_4LtabbrNm$llongNm)]
# make long sequence names to use for the sequences
# the alignment
df_lN02$lngsqNm <- paste(   df_lN02$accession_nmb,
                    df_lN02$llabbr,
                    df_lN02$smplyear2,
                    "Nmn",sep="_")
#_______________________________________________________________________________
# section 01 -  end - substitute in names and input files
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 02 -  start - read in fasta file with alignment, identify coding
# regions, trim alignment, and save only coding regions alignments
#_______________________________________________________________________________
# read in the alignment with sequences
al01 <- bioseq::read_fasta(pth_inpf01)
# remove any sequences that have been given the term 'copy' in the 
# sequence name
al01 <- al01[(!grepl("Copy|copy",names(al01)))]
# remove any sequence-names that appear in duplicates
al01 <- al01[!duplicated(names(al01))]
# Read in a fasta file with gene regions
Fas_Ref_seq <- bioseq::read_fasta(pth_inpf02)
# substitute the sequence names
lb01 <- labels(al01)
lb01 <- gsub("Mnemiopsis_leidyi_","",lb01)
lb01 <- gsub("_$","",lb01)
#lb01 <- gsub("__Copy$","",lb01)
lb01 <- gsub("_consensus_sequence","",lb01)
# identify the Mnelei sequences collected for this study
DGsmp <- lb01[grepl("Mnelei[0-9]{+}",lb01)]
# arrange the sequence names as a data frame
dfDGsmp <- data.frame(do.call('rbind', strsplit(as.character(DGsmp),'_',fixed=TRUE)))
colnames(dfDGsmp) <- c("smplNm","Lc1","Lc2","clcY","clcm")
dfDGsmp$clcm[grepl("Mnel",dfDGsmp$clcm)] <- "Jun"
# match to get short names
dfDGsmp$llabbr <- df_4LtabbrNm$llabbr[match(dfDGsmp$Lc2,df_4LtabbrNm$llongNm)]
# make replace names
rplcNms_DGsmp <- paste(dfDGsmp$smplNm,
      dfDGsmp$llabbr,
      dfDGsmp$clcY,
      dfDGsmp$clcm,sep="_")
# and use them to replace
lb01[grepl("Mnelei[0-9]{+}",lb01)] <- rplcNms_DGsmp

# identify the Mnelei sequences collected for this study
NCBIsmp <- lb01[!grepl("Mnelei[0-9]{+}",lb01)]
# match the accession number to get the long sequence names
NCBIsmp <- df_lN02$lngsqNm[match(NCBIsmp,df_lN02$accession_nmb)]
# and replace the NCBI labels
lb01[!grepl("Mnelei[0-9]{+}",lb01)] <- NCBIsmp
# use all the modified sequence names to get new 
# sequence header in the bioseq fasta file alignment read in
names(al01) <- lb01
# subset to exclude the samples from 'FKer_2017_May'
# as these stem from a different species
al01 <- al01[!grepl("FKer_2017_May",names(al01))]
# make a function to extract start and end positions
# of a gene region using a fasta input file (fasta_ref_seq) where
#  the headers have a genename (GnNm) to search for
# included, and specify the number of bases (noofbases)
# to search for in both ends of the gene
Get.st_en.of_geneR <- function(GnNm,fasta_ref_seq,noofbases,alignM){
  # make the fasta reference sequences  a tibble
  REF.seq <- dplyr::as_tibble(fasta_ref_seq)
  # get index numbers for GnNm
  idxN.GnNm <- which(grepl(GnNm,names(fasta_ref_seq)) )
  # get the reference seq for GnNm
  REF.seq.GnNm <- as.character(REF.seq[idxN.GnNm,])
  # get the first no of characters in the sequence
  REF.seq.GnNm_fnoof <- substr(REF.seq.GnNm,1,noofbases)
  # get the last no of characters in the sequence
  REF.seq.GnNm_lnoof <- substr(REF.seq.GnNm,nchar(REF.seq.GnNm)-noofbases,
                             nchar(REF.seq.GnNm))
  # use the longest sequence as reference
  # this will only work if the longest sequence covers the entire 
  # alignment - I know it is a rubbish function
  idxLng <- which(max(nchar(gsub("-","",alignM)))==nchar(gsub("-","",alignM)))
  # find the start position of the last and the first characters
  # of the longest sequence
  st.pos.end <- gregexpr(REF.seq.GnNm_lnoof, alignM[idxLng])[[1]][1]
  st.pos.sta <- gregexpr(REF.seq.GnNm_fnoof, alignM[idxLng])[[1]][1]
  st.pos.end <- st.pos.end+nchar(REF.seq.GnNm_lnoof)-1
  # return the start and end position of the gene name
  # in a list
  return(list(st.pos.sta,st.pos.end,REF.seq.GnNm_fnoof,REF.seq.GnNm_lnoof))
  
  }
# the first element in the list is the start position of the gene
# the second element is the end position of the gene
ITS1_s.p <- Get.st_en.of_geneR("ITS1",Fas_Ref_seq,12,al01)[[1]]
ITS1_e.p <- Get.st_en.of_geneR("ITS1",Fas_Ref_seq,12,al01)[[2]]
# use Bioseq 'seq_crop_pattern' function to cut out the block 
# with the coding ITS1 region
al01_ITS1 <- bioseq::seq_crop_position(al01,
                                       ITS1_s.p,
                                       ITS1_e.p)
# the function also returns the sequence in the first and last part
# of the region sought
ITS1_s.sq <- Get.st_en.of_geneR("ITS1",Fas_Ref_seq,12,al01)[[3]]
ITS1_e.sq <- Get.st_en.of_geneR("ITS1",Fas_Ref_seq,12,al01)[[4]]
# this can be pasted together and used to extract a pattern
gnpat <- paste0(ITS1_s.sq,"*.*",ITS1_e.sq)
# which pretty much return the same sequences, although not for 
# the sequences that are too disimilar
al01_ITS1p <- bioseq::seq_extract_pattern(al01,gnpat)

# the first element in the list is the start position of the gene
# the second element is the end position of the gene
ITS2_s.p <- Get.st_en.of_geneR("ITS2",Fas_Ref_seq,12,al01)[[1]]
ITS2_e.p <- Get.st_en.of_geneR("ITS2",Fas_Ref_seq,12,al01)[[2]]

# use Bioseq 'seq_crop_pattern' function to cut out the block 
# with the coding ITS2 region
al01_ITS2 <- bioseq::seq_crop_position(al01,
                                       ITS2_s.p,
                                       ITS2_e.p)
# use a function to find unique characters in a string
#https://stackoverflow.com/questions/31814548/function-that-extracts-each-unique-character-in-a-string
uniqchars <- function(x) unique(strsplit(x, "")[[1]]) 
# use this function inside another function
# to remove sequences that are missing all nucleotides
rmv.bnk_sq_from_algnm <- function(al01){
  # count the number of sequnces
  sqq<- al01 %>% as.data.frame(as.matrix()) %>% pull(1)
  nsq <- length(sqq)
  # make a number series to itersate over
  noFsq <- seq(1,nsq,1)
  # make an empty list to add to
  lst.unqchr <- list()
  # iterate over the series of numbers
  for (n in noFsq)
    { sqline <- sqq[n]
    #use the function to get the unique characters per sequence
    unqch <- uniqchars(sqline)
    #paste the characters found together into a single string
    chrFound <- paste(unqch, collapse = "")
    # then append this string to the list
    lst.unqchr[n] <- chrFound
  }
  #make the list a data frame
  df_usqqCh <- do.call(rbind.data.frame, lst.unqchr)
  # identify the row index numbers that do not have only '-'
  idx.mss <- which(df_usqqCh!="-")
  # then subset the data frame so it only comprises the sequences that 
  # at least have nucleotides
  al02 <- al01[idx.mss]
  return(al02)}

#then use the function on the two cut out gene regions.
al02_ITS1 <- rmv.bnk_sq_from_algnm(al01_ITS1)
al02_ITS2 <- rmv.bnk_sq_from_algnm(al01_ITS2)
# identify sequences that have 3 dashes, as these
# will be difficult to handle when making IUPAC code variants
# later on . - so it is better to remove them already now
idxNo.indelITS1 <- which(!grepl("---",as.matrix(al02_ITS1)[,1]))
idxNo.indelITS2 <- which(!grepl("---",as.matrix(al02_ITS2)[,1]))
al02_ITS1 <- al02_ITS1[c(idxNo.indelITS1)]
al02_ITS2 <- al02_ITS1[c(idxNo.indelITS2)]

# define file name to write to
outF.ITS1 <- paste0(wd00_wd05,"/algn_Mnelei_ITS1_v01.fasta")
outF.ITS2 <- paste0(wd00_wd05,"/algn_Mnelei_ITS2_v01.fasta")
#write the  alignment as a fasta file
bioseq::write_fasta(al02_ITS1,file=outF.ITS1,line_length = Inf,block_length = Inf)
bioseq::write_fasta(al02_ITS2,file=outF.ITS2,line_length = Inf,block_length = Inf)

#_______________________________________________________________________________
# Trim the alignment -  start
#_______________________________________________________________________________
# ------------- Options -------------
#set fraction OK for all other files 
FracOK <- 0.6
FracOK <- 0.99

# fraction of samples required to have data (NOT equal to "-") at position for it to be considered OK
fraction_ok <- FracOK
max_gap_width <- 100 # 
min_grp_width <- 100 #  
# values not in the list will be marked "others"" in the plots
valuelist <- c("a","g","c","t","n","-")
valuelist <- c("A","G","C","T","N","-")
othervalues <- "Others" 
# colours for plotting letters
plotcolors <- c("red",     # A
                "yellow",  # G
                "green",   # C
                "blue",    # T
                "grey",    # N   
                "black",   # -
                "white")   # others
# crop the sequence identifiers when plotting (limit to ncrop characters)
ncrop <- 20 
# ensure packages are loaded
library(tidyverse)
library(bioseq)
library(ape)
# get the working dir
prwd<- getwd()
# set working directory to be able to read in functions
setwd(paste0(wd00,"/suppmat03_Rcodes"))
source("Rfunction_DNA_sequence_subset.R")
source("Rfunction_ReadFasta.R")
setwd(prwd)
# ------------- Load data using ReadFasta function -------------
df <- ReadFasta(outF.ITS2)
# df <- seqinr::read.fasta(file = pth_inpf04)
# ------------- Calculate positions to be included (OK)  -------------
res <- DNA_sequence_subset(df,
                           required_fraction=fraction_ok,
                           max_gap_width=max_gap_width,
                           min_grp_width=min_grp_width)

df_sum <- res[[1]]
sListOK <- res[[2]]

# join df_sum to the input data marking which postions should be included
df <- df %>%
  left_join(select(df_sum,Posn,OK),by="Posn")

df_samples <- df %>%
  distinct(Sample)

nsamples <- nrow(df_samples)
nrequired <- fraction_ok*nsamples

# ------------- Output dataframe -------------
# output dataframe with sequences where OK==TRUE
df_out <- df %>% 
  filter(OK==TRUE) %>%
  select(-OK)
df_out <- df_out %>%
  pivot_wider(values_from=Value,names_from=Posn)

# 
# filename<-paste0(substr(inpf01,1,nchar(inpf01)-4),"_cropped.csv")
# folder_out <- wd00_wd05
# write.table(df_out,file=paste0(folder_out,filename),row.names=F,col.names=T,sep=",",quote=F)


# ------------- Plot count of OK data at each position -------------
xlabel <- paste0("Position (OK: ",sListOK,")")

# plot the number of OKs at each position
p1 <- ggplot(data=df_sum,aes(x=Posn,y=Count)) +
  geom_point()+ 
  geom_hline(yintercept=nrequired,color="red",linetype=2) +
  ylab("Count")+
  xlab(label="") + 
  geom_ribbon(aes(ymin=0, ymax=(1-OK)*nsamples),alpha=0.1) +
  theme_minimal()

p1

df_plot <- df %>%
  mutate(Value=ifelse(Value %in% valuelist,Value,othervalues))

df_plot$Sample <- factor(df_plot$Sample,levels=rev(df_samples$Sample),labels=rev(substr(df_samples$Sample,1,ncrop)))
df_plot <- df_plot %>%
  mutate(a=ifelse(OK==T,1,0.7))

df_plot$Value <- factor(df_plot$Value,levels=c(valuelist,othervalues))
plotlabels <- c(toupper(paste0(valuelist,"     ")),othervalues)

p2 <- ggplot() +
  geom_raster(data=df_plot,aes(x=Posn,y=Sample,fill=Value)) + #
  ylab("Sample")+
  xlab(label="") + 
  scale_fill_manual(values=plotcolors,name="",labels=plotlabels) +
  scale_alpha_continuous(guide="none") +
  theme_minimal() +
  theme(axis.text.y=element_text(size=5),
        legend.text = element_text(size=8),
        legend.position="bottom",
        legend.key=element_rect(color="black"),
        legend.key.size=unit(0.6,"lines")) +
  guides(fill=guide_legend(nrow=1))
p2

#____________________________________________________________
# END trim the alignment
#____________________________________________________________
#_______________________________________________________________________________
# section 02 -  end - read in fasta file with alignment, identify coding
# regions, trim alignment, and save only coding regions alignments
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 03 -  start - make all disambiguity variants, and rename
#_______________________________________________________________________________


#https://www.biostars.org/p/158250/
# Create all possible variants from ambiguities
aITS1_vrs <- bioseq::seq_disambiguate_IUPAC(al02_ITS1)
aITS2_vrs <- bioseq::seq_disambiguate_IUPAC(al02_ITS2)
#_______________________________________________________________________________
# function : 'Rnm_seq_var' - start
#_______________________________________________________________________________
# make a function that can rename the sequence variants that 
# was produced from the 'bioseq::seq_disambiguate_IUPAC' function
Rnm_seq_var <- function(M.algn)
{
  # count the number of entries in the list of sequences
  nsq <- length(M.algn)
  # make series of number that covers this
  # number of entries in the list of sequences
  sqNoAlgn <- seq(1,nsq,1)
  # make an empty list to add variants of sequences to
  lst_vseq <- list()
  # iterate over sequences
  for (n in sqNoAlgn)
  {
    disam.s <- M.algn[n]
    unldisam.s <- unlist(disam.s)
    unldisam.s <- as.character(unldisam.s)
    sqNm <- names(disam.s)
    disam.s1 <- unlist(disam.s)
    n.disam <- length(disam.s1)
    sqn.disam <- seq(1,n.disam,1)
    vrNms <- paste0(sqNm,"_var",LETTERS[sqn.disam])
    df_vNmsq <- as.data.frame(cbind(vrNms,unldisam.s))
    lst_vseq[[n]] <- df_vNmsq
  }
  #bind the rows in each list in to one data frame
  df_disam <- data.table::rbindlist(lst_vseq, fill=T)
  df_disam <- as.data.frame(df_disam)
  # the data frame cannot be made a bioseq tibble directly
  # use only the columns with the sequences to make a bioseq
  # object with disambiguity sequences 
  bsq_disam <- bioseq::as_dna(df_disam[,2])
  # then make this a DNAbin object - however this DNAbin object need names
  dnb_disam <- bioseq::as_DNAbin(bsq_disam)
  # get the names from column 1 in the dataframe
  names(dnb_disam) <- df_disam[,1]
  # make the DNAbin object a bioseq-DNA-object
  bsq_disam <- bioseq::as_dna(dnb_disam)
  return(bsq_disam)
}
#_______________________________________________________________________________
# function : 'Rnm_seq_var' - end
#_______________________________________________________________________________

# use this function to rename the disambiguity sequences prepared
alvrs.M.ITS1 <- Rnm_seq_var(aITS1_vrs)
alvrs.M.ITS2 <- Rnm_seq_var(aITS2_vrs)


# define file name to write to
outF.ITS1 <- paste0(wd00_wd05,"/algn_Mnelei_ITS1_v02_all_IUPACode_vars.fasta")
outF.ITS2 <- paste0(wd00_wd05,"/algn_Mnelei_ITS2_v02_all_IUPACode_vars.fasta")
#write the  alignment as a fasta file
bioseq::write_fasta(alvrs.M.ITS1,file=outF.ITS1,line_length = Inf,block_length = Inf)
bioseq::write_fasta(alvrs.M.ITS2,file=outF.ITS2,line_length = Inf,block_length = Inf)

#_______________________________________________________________________________
# section 03 -  end - make all disambiguity variants, and rename
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 04 -  start - make trees from disambiguity variants
#_______________________________________________________________________________
#
#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)
#iterate over the  sequence of numbers that  reflect the genes 
for (ng in n.f.alvrs.M)
{
  print(ng)
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm<- gene.nms[ng]
  # Now the alignment has been picked from the list,
  # Then prepare the NJ trees, and the label categories
  # and a color range to assign the labels in the trees
  # also get the labels on the sequences
  LalvrsM <- names(alvrsM)
  alvrsM <- bioseq::as_DNAbin(alvrsM)
  names(alvrsM) <- LalvrsM
  # Try and plot a neighbour join tree
  mtr_dd_pip <- ape::dist.dna(alvrsM)
  dmtr_dd_pip <- as.data.frame(as.matrix(mtr_dd_pip))
  # replace NAn with zero
  mtr_dd_pip[is.na(mtr_dd_pip)] <- 0
  tre_pip <- ape::njs(mtr_dd_pip)
  plot(tre_pip)
  #sort the branches in the tree
  tre_pipr <- ape::ladderize(tre_pip, right = TRUE)
  #https://joey711.github.io/phyloseq/plot_tree-examples.html
  plot(tre_pipr, cex=0.4)
  
  library(ggtree)
  p <- ggtree(tre_pipr, right = TRUE,
              options(ignore.negative.edge=TRUE)) + 
    xlim(0, 0.025) + # to allow more space for labels
    geom_treescale() # adds the scale
  
  df_tiplb01 <- as.data.frame(cbind(c(tre_pipr$tip.label)))
  df_tiplb01$cat <- NA
  colnames(df_tiplb01) <- c("seqNm", "cat")
  
  # arrange the sequence names as a data frame
  df_tiplb02 <- data.frame(do.call('rbind', strsplit(as.character(df_tiplb01$seqNm),'_',fixed=TRUE)))
  colnames(df_tiplb02) <- c("smplNm","LocNm","smplYer","smplMn","varNm")
  seqNm <- df_tiplb01$seqNm
  df_tiplb01 <- cbind(seqNm,df_tiplb02)
  
  #________________
  # make a vector with the column names to color lables by as categories in the
  # NJ tree
  lblvars <- c("LocNm","varNm")
  # also make a vector with longer and more descriptive names for the variables
  # that are to be used in the legend title
  lblvars.LngNm <- c("Locaction","Gene variant")
  # make an empty list to put generated plots into
  lst_plts <- list()
  # iterate over the columns that are to be used as tip label categories
  for (l in lblvars)
  {
    
    #subset data frame
    df_tiplb02 <- df_tiplb01[,c("seqNm",l)]
    #rename columns
    colnames(df_tiplb02) <- c("seqNm", "cat")
    plt.idx.n <- which(grepl(l,lblvars))
    # also get e legend title replacemnet
    lg.title.rpl<- lblvars.LngNm[plt.idx.n]
    # rename the data frame
    tipcategories <- df_tiplb02
    # make colors for the text inside the labels, this requires a small
    # matrix with headers for text labels and values as colors
    tpcats <- unique(tipcategories[,2])
    tpcats <- tpcats[order(tpcats)]
    n.tpcats <- length(tpcats)
    # it is about 60% of the colurs in the color range that are
    # dark colours and require a white text on the label
    # by rounding the proportion, an amount of catagories can be
    # assigned white text
    no.white <- round(n.tpcats*0.6,digits = 0)
    clfltxtlb <- rev(c(rep("black",(n.tpcats-no.white)),rep("white",no.white)))
    # make the color categories a 
    tplbcls <- rbind(tpcats,clfltxtlb)
    colnames(tplbcls) <- tplbcls[1,]
    txtlblcols <-tplbcls[-1,]
    # make the tipcategories a dataframe
    dd <- as.data.frame(tipcategories)
    plt_titl_LET <- LETTERS[plt.idx.n]
    # check out color scales here: https://sjspielman.github.io/introverse/articles/color_fill_scales.html
    plt_n <- p %<+% dd + 
      ggtree::geom_tiplab(aes(fill = factor(cat),
                              color= factor(cat)),
                          #color = "black", # color for label font
                          geom = "label",  # labels not text
                          label.padding = unit(0.09, "lines"), # amount of padding around the labels
                          label.size = 0, # size of label border
                          size=2.4) + 
      #scale_fill_viridis_d(option = "cividis") +
      scale_fill_viridis_d(option = "viridis") +
      #scale_color_viridis_d(option = "G", direction = -1) +
      scale_color_manual(values = txtlblcols) +
      
      labs(fill=lg.title.rpl) +
      labs(color=lg.title.rpl) +
      # theme(#legend.position = c(0.5,0.2),
      #      legend.title = element_text("lege.Titl") ) +
      #   #, #element_blank(), # no title
      #      #legend.key = element_blank()) # no keys +
      #use a title
      ggtitle(plt_titl_LET)
    # add the plot to the list that is for gathering plots
    lst_plts[[plt.idx.n]] <- plt_n
    # end iteration over columns to make variables from that are used
    # to make plots
}

library(patchwork)
# assemble plots with patchwork
plt_assmbl <- lst_plts[[1]] + lst_plts[[2]]
#make filename to save plot to
flNm.Fig <- paste0(wd00_wd05,"/Fig01_NJ_tree_disambig_seq_",GnNm,"_v01.png")
# make an if test to check whether the plot should be saved
# i.e. set to FALSE if you do not want the plot to be saved
bSaveFigures=T
# if the 'bSaveFigures' is TRUE the plot will be saved
if(bSaveFigures==T){
  ggsave(plt_assmbl,file=flNm.Fig,
         width=210,height=297,
         #width=210,height=(297*0.5),
         units="mm",dpi=150)
}
#________________
# end iteration over alignments
}
#_______________________________________________________________________________
# section 04 -  end - make trees from disambiguity variants
#_______________________________________________________________________________