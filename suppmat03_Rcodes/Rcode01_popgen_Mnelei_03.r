#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# This code is able  to run in:
# version
# platform       x86_64-pc-linux-gnu         
# arch           x86_64                      
# os             linux-gnu                   
# system         x86_64, linux-gnu           
# status                                     
# major          3                           
# minor          6.2                         
# year           2019                        
# month          12                          
# day            12                          
# svn rev        77560                       
# language       R                           
# version.string R version 3.6.2 (2019-12-12)
# nickname       Dark and Stormy Night 

#NOTICE that the 'haplotypes' package requires 
# R v3.6.2

# To get 'rstudio' to run R v3.6.2.
# First install R v3.6.2. 
#The in a unix terminal
# Execute:
#$ export RSTUDIO_WHICH_R=/usr/bin/R-3.6.2/bin/R
#The start 'rstudio' from the terminal by typing:
#$ rstudio

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
#read in packages and install libraries
#read in the 'ape' library
library(ape)
#read in the 'pegas' library
library(pegas)

#install packages if required
# gaston had issues installing. Returning the error: “failed to create lock directory”
# I looked here: https://stackoverflow.com/questions/14382209/r-install-packages-returns-failed-to-create-lock-directory
if(!require("gaston")){
  install.packages("gaston", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("gaston")
}
# The 'hierfstat' package has many functions for various pop gen 
# calculations
if(!require("hierfstat")){
  install.packages("hierfstat", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("hierfstat")
}
# the 'haplotypes' package is used for calculating Phist values
# unfortunately 'haplotypes will not work in R v3.3, and will
# not work in R v4.0.2, which is why this code is prepared for
# R v3.6.2
# Here the installation of the 'haplotypes' package
# is commented out, as the 'hierfstat' package
# appears to be capable of doing the same thing
# if(!require("haplotypes")){
#   install.packages("haplotypes")
#   library("haplotypes")
# }

# get ips packge to write dnabin to fasta files
if(!require("ips")){
  install.packages("ips", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("ips")
}
library("ips")
#The 'tidyverse' package is required for rearranging the data frames before
# making the colored tables
if(!require("tidyverse")){
  install.packages("tidyverse", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("tidyverse")
}
# the 'pals' enables you to make color gradients - here used for colored tables
if(!require("pals")){
  install.packages("pals", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("pals")
}
# strataG package can calculate Tajimas D,

#But the 'strata' package is not available for R v3.6.2
# The strataG package works for R v4.0.2
# if(!require("strataG")){
#   install.packages("strataG", dependencies = TRUE, INSTALL_opts = '--no-lock')
#   library("strataG")
# }

#read in libraries
library("apex")
library("adegenet")
library("pegas")
library("mmod")
#library("poppr")
library("hierfstat")
library("tidyverse")
library("pals")

library(poppr)
library(ape)
library(RColorBrewer)
library(vegan)
# get the poppr package
if(!require(poppr)){
  install.packages("poppr", repos='http://cran.us.r-project.org')
}
library(poppr)
# get the vegan package
if(!require(vegan)){
  install.packages("vegan", repos='http://cran.us.r-project.org')
}
library(vegan)
# get the stringi package which is required for the adegenet package
if(!require(stringi)){
  install.packages("stringi", repos='http://cran.us.r-project.org')
}
library(stringi)
# get the adegenet package -  will allow you to work with genlight objects
if(!require(adegenet)){
  install.packages("adegenet", repos='http://cran.us.r-project.org')
}
library(adegenet)

# if not installed , then install package: PopGenReport -  you will need this package
# for the SNPRelate package and the dartR package
# the cran deposited version fails to install, so instead get it from github
if(!require(PopGenReport)){
  #install.packages("devtools")
  devtools::install_github("green-striped-gecko/PopGenReport")
}
library(PopGenReport)
# get the SNPRelate package - which is needed for the dartR package
# the cran deposited version fails to install, so instead get it from github
if(!require(SNPRelate)){
  #install.packages("devtools")
  devtools::install_github("zhengxwen/SNPRelate")
}
library(SNPRelate)

# get the dartR package -  can convert genind objects to genlight objects
if(!require(dartR)){
  install.packages("dartR")
  gl.install.vanilla.dartR()
  }
library(dartR)
library(biogeo)


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
inpf01 <- "algn_Mnelei_18s_09.fas"
inpf01 <- "algn_Mnelei_18s_10.aligned.fasta.fas"
inpf02 <- "output01_accn_publication_Mnemiopsis_its1its2.csv"
pth_inpf01 <- paste(wd00_wd01,"/",inpf01,sep="")
pth_inpf02 <- paste(wd00_wd01,"/",inpf02,sep="")
#read in csv file
df_lN02 <- read.csv(pth_inpf02,sep=";")
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
unique(df_lN02$pubjournal)
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

#Rta01.1

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
# ovewrite the incorrect position for the USA Woods Hole sample 'AF293700'
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
# Assign location for the Holland sample
df_lN02$location[grepl("^EF175463",df_lN02$accession_nmb)] <- "Holland"

#replace in location names
df_lN02$location[grepl("Bornholm",df_lN02$location)] <- "Baltic Sea"
df_lN02$location[grepl("Denmark",df_lN02$location)] <- "NE Atlantic"
df_lN02$location[grepl("Iran",df_lN02$location)] <- "Caspian Sea"
df_lN02$location[grepl("Caspian",df_lN02$location)] <- "Caspian Sea"
df_lN02$location[grepl("Norway",df_lN02$location)] <- "NE Atlantic"
df_lN02$location[grepl("Mediterranean",df_lN02$location)] <- "Mediterranean"
# replace lat and lon for the NA entries
df_lN02$latlonpos_nd[is.na(df_lN02$latlonpos_nd)] <- df_lN02$latlonpos_nd2[is.na(df_lN02$latlonpos_nd)]
# split a column
df_lN02 <- tidyr::separate(df_lN02, latlonpos_nd, sep = " ", into = paste0("ll", 1:6), fill = "right")
#get decimal degrees
dd <- as.numeric(df_lN02$ll1)
mm <- as.numeric(df_lN02$ll2)
ns <- as.character(df_lN02$ll3)
df_lN02$lat1 <- biogeo::dms2dd(dd,mm,0,ns)
#get decimal degrees
dd <- as.numeric(df_lN02$ll4)
mm <- as.numeric(df_lN02$ll5)
ns <- as.character(df_lN02$ll6)
df_lN02$lon1 <- biogeo::dms2dd(dd,mm,0,ns)
# split a column
df_lN02 <- tidyr::separate(df_lN02, latlonpos, sep = " ", into = paste0("lld", 1:4), fill = "right")
#evaluate hemisphere of coordinate
lld1 <- as.numeric(df_lN02$lld1)
df_lN02$lat2 <- ifelse(df_lN02$lld2=="N" | df_lN02$lld2=="S",(ifelse(df_lN02$lld2=="N",lld1,-lld1)),0)
#evaluate hemisphere of coordinate
lld3 <- as.numeric(df_lN02$lld3)
df_lN02$lon2 <- ifelse(df_lN02$lld4=="E" | df_lN02$lld4=="W",(ifelse(df_lN02$lld4=="E",lld3,-lld3)),0)
# check if column has NA and if if does, then add from the other column
df_lN02$lat1[is.na(df_lN02$lat1)] <- df_lN02$lat2[is.na(df_lN02$lat1)]
df_lN02$lon1[is.na(df_lN02$lon1)] <- df_lN02$lon2[is.na(df_lN02$lon1)]
# make it a new column
df_lN02$declat <- df_lN02$lat1
df_lN02$declon <- df_lN02$lon1
# replace any spaces in location names
df_lN02$location <- gsub(" ","",df_lN02$location)
# add a column for sample year
df_lN02$smplyear <- gsub(".*(20[0-9]{2}).*","\\1",df_lN02$pubjournal)
df_lN02$smplyear[grepl("Unpublished",df_lN02$smplyear)] <- "unknown"




#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#define input file
fl1 <- "collect_loc_Mnelei_smpls01.csv"
#paste path and filename together
pth_fl01 <- paste(wd00_wd01,"/",fl1,sep="")
#read in a csv file with positions for sampling locations
df_clo <- as.data.frame(read_delim(pth_fl01,delim = ";"))

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
#df_clo$locality4 <- gsub("SEDenmarkMecklenburgerBucht","", df_clo$locality4 )
df_clo$locality4 <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", df_clo$locality4 )

# df_clo
#use biogeo dms2dd funciton to convert to decimal degrees
df_clo$dec_lat <- biogeo::dms2dd(df_clo$lat_deg,
                                 df_clo$lat_min,
                                 df_clo$lat_sec,
                                 df_clo$lat_sph)
#use biogeo dms2dd function to convert to decimal degrees
df_clo$dec_lon <- biogeo::dms2dd(df_clo$lon_deg,
                                 df_clo$lon_min,
                                 df_clo$lon_sec,
                                 df_clo$lon_sph)
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
  c("FBo",
    "FKe",
    "JMa",
    "GKi",
    "GMe",
    "JLi",
    "JLo",
    "GHe",
    "GWa",
    "SBa",
    "DMe",
    "Ssk")
#combine to a data frame
df_abblocnms <- as.data.frame(cbind(llocnms,abbloccnms))
#substitue in data frame
df_abblocnms$locnm2 <- gsub(" ","",df_abblocnms$llocnms)
#copy data frame
df_clo03$locality7 <- df_clo03$locality2

# # round decimals in lat lon positions to and paste together
df_lN02$declat2 <- round(df_lN02$declat,3)
df_lN02$declon2 <- round(df_lN02$declon,3)
df_lN02$latlonpos3 <- paste(df_lN02$declat2,"; ",df_lN02$declon2,sep="")
#match back lat lon positions to data frame
#df_ihpt03$latlon3 <- df_lN02$latlonpos3[match(df_ihpt03$spcAccn1,df_lN02$accession_nmb)]
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
# substitute in locality names
df_clo03$locality4 <- gsub("Bogense,","Bogense",df_clo03$locality4)
# copy column in data frame
df_clo03$locality8 <- df_clo03$locality7
# grep for long location names and replace with abbreviations
df_clo03$locality8[grepl("Bogense",df_clo03$locality4)] <- "FBo"
df_clo03$locality8[grepl("Kerteminde",df_clo03$locality4)] <- "FKe"
df_clo03$locality8[grepl("North Sea, Helgoland",df_clo03$locality4)] <- "GHe"
df_clo03$locality8[grepl("Kiel Fjord",df_clo03$locality4)] <- "GKi"
df_clo03$locality8[grepl("Limfjord",df_clo03$locality4)] <- "JLi"
df_clo03$locality8[grepl("Mariager",df_clo03$locality4)] <- "JMa"
df_clo03$locality8[grepl("Samsøe",df_clo03$locality4)] <- "SBa"
df_clo03$locality8[grepl("Skovshoved",df_clo03$locality4)] <- "SSk"
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
# Read in the fasta file
#_______________________________________________________________________________
#_______________________________________________________________________________
#read FASTA as dna.bin
pip <- ape::read.dna(pth_inpf01, format = "fasta")
#make the DNAbin object a genind object
geni_pip <- adegenet::DNAbin2genind(pip)
#make the genind object a dataframe
df_pip <- adegenet::genind2df(geni_pip)
#get the row names
orig_rwnm <- row.names(df_pip)
#substitute in the row names
mdf_rnm01 <- gsub("_consensus_sequence","",orig_rwnm)
mdf_rnm02 <- gsub("Mnelei","Mne_lei", mdf_rnm01)
mdf_rnm03 <- gsub("Mnemiopsis_leidyi","Mne_lei", mdf_rnm02)
#nms_index_1st2ndprt <- gsub("^([^_]*_[^_]*)_.*$", "\\1", nms_index)
#nms_index_2ndprt <- gsub("^([^_]*)(_[^_]*)_.*$", "\\2", nms_index)
#substitute to get year and month
df_pip$rwnm04 <- gsub("(.*)_(.*)_(.*)$","\\2_\\3",mdf_rnm03)
#replace 
df_pip$rwnm04 <- gsub("Loegstoer_","",df_pip$rwnm04)
#substitute to get year 
df_pip$rwnm06 <- gsub("(.*)_(.*)$","\\1",df_pip$rwnm04)
#substitute to get location
df_pip$rwnm05 <- gsub("(.*)_(.*)_(.*)_(.*)_(.*)$","\\2_\\3",mdf_rnm03)
#substitute to get specific number
mdf_rnm04 <- gsub("Mne_lei_","Mne_lei",mdf_rnm03)
mdf_rnm04  <- gsub("Mne_lei","Mnelei",mdf_rnm03)
df_pip$rwnm07 <- gsub("(.*)_(.*)_(.*)_(.*)_(.*)$","\\1\\2",mdf_rnm04)
df_pip$rwnm07 <- gsub("(Mnelei[0-9]{3})(.*)$","\\1",df_pip$rwnm07)
df_pip$rwnm07 <- gsub("_","",df_pip$rwnm07)
df_pip$rwnm05 <- gsub("Mne_lei_","",df_pip$rwnm05)
#substitute some of the odd names
df_pip$rwnm05[grepl("_Limfjord", df_pip$rwnm05)] <- gsub("Limfjord","NJylland_Limfjord",df_pip$rwnm05[grepl("_Limfjord", df_pip$rwnm05)])
df_pip$rwnm05[grepl("Limfjord", df_pip$rwnm05)] <- gsub("Limfjord_Loegstoer","NJylland_LimfjordLogstoer",df_pip$rwnm05[grepl("Limfjord", df_pip$rwnm05)])
df_pip$rwnm05 <- gsub("lei[0-9]{3}_","",df_pip$rwnm05)
# make the synonym sample locations identical in name
df_pip$rwnm05 <- gsub("SEDenmark_MecklenburgerBucht","NGermany_MecklenburgerBuchtWismarBucht",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("NJyllandLimfjordLogstoer","NJyllandLimfjord",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("Limfjord_Loegstoer","Limfjord",df_pip$rwnm05)
#unique(df_pip$rwnm05)
# #Modify location names
# df_pip$rwnm05 <- gsub("FynBogense","Funen, Bogense", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("FynKerteminde","Funen, Kerteminde", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("JyllandMariagerfjord","Jutland, Mariager Fjord", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("NGermanyKielFjord","Germany, Kiel Fjord", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("NGermanyMecklenburgerBuchtWismarBucht","Germany, Wismar Bight", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("NJyllandLimfjord","Jutland, Limfjord-Aalborg", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("NJyllandLimfjordLogstoer","Jutland, Limfjord-Løgstør" , df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("NWGermanyNSeaHelgolandRds", "North Sea, Helgoland Roads", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("NWGermanyWaddenSeaBussumHaupstr","Germany, Büsum", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("SamsoeBallen","Samsøe, Ballen", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("SEDenmarkMecklenburgerBucht","Germany, Wismar Bight", df_pip$rwnm05 )
# df_pip$rwnm05 <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", df_pip$rwnm05 )

df_pip$rwnm04 <- gsub("_","",df_pip$rwnm04)
df_pip$rwnm05 <- gsub("_","",df_pip$rwnm05)
df_pip$rwnm06 <- gsub("_","",df_pip$rwnm06)
df_pip$rwnm07 <- gsub("_","",df_pip$rwnm07)
df_pip$rwnm09 <- df_pip$rwnm05
df_pip$rwnm10 <- df_pip$rwnm06
# replace accession numbers with locations
df_pip$rwnm09[!is.na(match(df_pip$rwnm05,df_lN02$accession_nmb))] <- df_lN02$location[!is.na(match(df_pip$rwnm05,df_lN02$accession_nmb))]
# replace accession numbers with sample year
df_pip$rwnm10[!is.na(match(df_pip$rwnm05,df_lN02$accession_nmb))] <- df_lN02$smplyear[!is.na(match(df_pip$rwnm05,df_lN02$accession_nmb))]
# df_pip$rwnm07 <- gsub("(MneleiHM)(.*)$","\\1_HMNCBI",df_pip$rwnm07)
# df_pip$rwnm07 <- gsub("(MneleiJQ)(.*)$","\\1_JQNCBI",df_pip$rwnm07)
# df_pip$rwnm07 <- gsub("(MneleiGU)(.*)$","\\1_GUNCBI",df_pip$rwnm07)
# df_pip$rwnm07 <- gsub("(MneleiAF)(.*)$","\\1_AFNCBI",df_pip$rwnm07)
# unique(df_pip$rwnm07)
df_pip$rwnm09 <- gsub("Germany:KielFjord" ,"NGermanyKielFjord" ,df_pip$rwnm09)
df_pip$rwnm09 <- gsub("Germany:Helgoland" ,"NWGermanyNSeaHelgolandRds",df_pip$rwnm09)
df_pip$rwnm09 <- gsub("USA:GalvestonBay" ,"USA:WoodsHole",df_pip$rwnm09)
df_pip$rwnm09 <- gsub("USA:Panacea" ,"USA:WoodsHole",df_pip$rwnm09)

df_pip$rwnm05 <- df_pip$rwnm09 
df_pip$rwnm06 <- df_pip$rwnm10
df_pip$rwnm06 <- gsub("(HM)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(JQ)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(GU)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(AF)(.*)$","\\1NCBI",df_pip$rwnm06)

df_pip$rwnm06 <- gsub("(JX)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(EF)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(KF)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(KJ)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(KY)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(OM)(.*)$","\\1NCBI",df_pip$rwnm06)
df_pip$rwnm06 <- gsub("(OU)(.*)$","\\1NCBI",df_pip$rwnm06)
#unique(df_pip$rwnm06)

df_pip$rwnm05 <- gsub("(HM)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(JQ)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(GU)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(AF)(.*)$","\\1NCBI",df_pip$rwnm05)

df_pip$rwnm05 <- gsub("(EF)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(JX)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(KF)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(KJ)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(KY)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(OM)(.*)$","\\1NCBI",df_pip$rwnm05)
df_pip$rwnm05 <- gsub("(OU)(.*)$","\\1NCBI",df_pip$rwnm05)
#unique(df_pip$rwnm05)
df_pip$rwnm04 <- gsub("(HM)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(JQ)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(GU)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(AF)(.*)$","\\1NCBI",df_pip$rwnm04)

df_pip$rwnm04 <- gsub("(EF)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(JX)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(KF)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(KJ)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(KY)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(OM)(.*)$","\\1NCBI",df_pip$rwnm04)
df_pip$rwnm04 <- gsub("(OU)(.*)$","\\1NCBI",df_pip$rwnm04)
#unique(df_pip$rwnm04)

df_pip$rwnm08 <- paste(df_pip$rwnm07,df_pip$rwnm05,df_pip$rwnm06,df_pip$rwnm04,sep="_")
#unique(df_pip$rwnm08)
df_pip$rwnm08 <- gsub("NJyllandLimfjordLogstoer","NJyllandLimfjord",df_pip$rwnm08)
#df_pip$rwnm08
#replace the row names in the data frame
row.names(df_pip) <- mdf_rnm03

row.names(df_pip) <- df_pip$rwnm08
df_pip$rwnm02 <- row.names(df_pip)

# subset to only include the rows that match
df_pip02 <- subset(df_pip, grepl("Mnelei", rwnm02))
#row.names(df_pip02)
#delete the column no longer needed
df_pip02$rwnm02 <- NULL
df_pip02$rwnm03 <- NULL
df_pip02$rwnm04 <- NULL
df_pip02$rwnm05 <- NULL
df_pip02$rwnm06 <- NULL
df_pip02$rwnm07 <- NULL
df_pip02$rwnm08 <- NULL
# replace all NAs
df_pip02 <- df_pip02 %>% replace(is.na(.), "-")
#make the date frame a matrix and a DNAbin object again
dnb_pip02 <- as.DNAbin(as.matrix(df_pip02))

pip <- dnb_pip02
#________________________________________________________________
#start - Plot haplotype network:
#________________________________________________________________

#create haplotypes from dna.bin
pipHaps <- pegas::haplotype(pip)
#view haplotype 
as.matrix(pipHaps)
#prepare hpt table
ind.hap<-with(
	stack(setNames(attr(pipHaps, "index"), rownames(pipHaps))),
	table(hap=ind, pop=rownames(pip)[values]))
#make it a dataframe
df_ihpt01 <- as.data.frame(ind.hap)
#limit to include only 'Freq' that equals 1
df_ihpt02 <- df_ihpt01[df_ihpt01$Freq == 1,]	
#check object type
class(df_ihpt02)
#view(df_ihpt02)
# split string and get lists nested in a list
locations <- strsplit(as.character(df_ihpt02$pop), "_")
# AF293700.1 comes from the study:
#   Podar M, Haddock SH, Sogin ML, Harbison GR. A molecular phylogenetic framework for the phylum Ctenophora using 18S rRNA genes. Mol Phylogenet Evol. 2001 Nov;21(2):218-30. doi: 10.1006/mpev.2001.1036. PMID: 11697917.
# 
# GU062752.1 comes from the study:
#   Ghabooli,S., Shiganova,T.A., Zhan,A., Cristescu,M.E., Eghtesadi-Araghi,P. and MacIsaac,H.J. (unpublished). Multiple introductions and invasion pathways for the invasive ctenophore Mnemiopsis leidyi in Eurasia. Unpublished
# 
# HM007195.1 comes from the study:
#   Fuentes,V., Angel,D.L., Bayha,K.M., Atienza,D., Edelist,D., Bordehore,C. and Purcell,J.E. (2010). Blooms of the invasive ctenophore, Mnemiopsis leidyi, span the Mediterranean Sea in 2009.  Hydrobiologia (2010). In press.
#
# JQ071530.1 comes from the study:
#   Saponari,L. (Unpublished). First molecular species identification of Mnemiopsis leidyi in Italian Seas. Unpublished

#get second object in nested list
#unique sample number
smplno <- sapply(locations, "[[", 1)
#locality
smplloca <- sapply(locations, "[[", 2)
#unique(smplloca)
#Modify location names
smplloca <- gsub("FynBogense","Funen, Bogense", smplloca )
smplloca <- gsub("FynKerteminde","Funen, Kerteminde", smplloca )
smplloca <- gsub("JyllandMariagerfjord","Jutland, Mariager Fjord", smplloca )
smplloca <- gsub("NGermanyKielFjord","Germany, Kiel Fjord", smplloca )
smplloca <- gsub("NGermanyMecklenburgerBuchtWismarBucht","Germany, Wismar Bight", smplloca )
smplloca <- gsub("NJyllandLimfjord","Jutland, Limfjord", smplloca )
#smplloca <- gsub("NJyllandLimfjordLogstoer","Jutland, Limfjord-Løgstør" , smplloca )
smplloca <- gsub("NWGermanyNSeaHelgolandRds", "North Sea, Helgoland Roads", smplloca )
smplloca <- gsub("NWGermanyWaddenSeaBussumHaupstr","Germany, Büsum", smplloca )
smplloca <- gsub("SamsoeBallen","Samsøe, Ballen", smplloca )
#smplloca <- gsub("SEDenmarkMecklenburgerBucht","", smplloca )
smplloca <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", smplloca )

accnms <- gsub("Mnelei","",sapply(locations, "[[", 1))
#df_lN02$smplyear[match(accnms,df_lN02$accession_nmb)]

# only year
smplye <- sapply(locations, "[[", 3)
# year and month
smplym <- sapply(locations, "[[", 4)
#class(smplym)
# see the first header of this list of characters
#head(smplym)
#make it a table
new.hap.smplloc <- table(df_ihpt02$hap, smplloca)
new.hap.smplye <- table(df_ihpt02$hap, smplye)
# check what kind of object it is
#class(locations)
#see the object
#new.hap
# make it a data frame
df_nhpt01.loc <- as.data.frame(new.hap.smplloc)
df_nhpt01.ye  <- as.data.frame(new.hap.smplye)
#class(df_nhpt01.ye)
#use reshape2 package to recast df from long to wide: see this example : https://www.datacamp.com/community/tutorials/long-wide-data-R
df_nhpt02.loc <- reshape2::dcast(df_nhpt01.loc, Var1 ~ smplloca, value.var="Freq")
df_nhpt02.ye <- reshape2::dcast(df_nhpt01.ye, Var1 ~ smplye, value.var="Freq")

#dev.off()
#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 1.2, cex = 0.1, pie = new.hap.smplloc, 
     show.mutation = 0, threshold = 0, labels(FALSE))
#add a legend to the plot
legend("topright",colnames(new.hap.smplloc), 
       col=rainbow(ncol(new.hap.smplloc)), 
       pch=19, ncol=1)

#define variable
inp.f.fnm <- "Mnelei_per_location"
flnm <- c(paste("haplotype_network_",inp.f.fnm,".pdf",  sep = ""))
#paste output directory and filename together in a string
outflnm <- paste(wd00_wd05,"/",flnm,sep="")
# Exporting PFD files via postscript()           
pdf(outflnm,
    width=(1*1.0*8.2677),height=(2*1.0*2.9232))

#define the plotting area with par settings
op <- par(mfrow=c(1,1), # set number of panes inside the plot - i.e. c(2,2) would make four panes for plots
          oma=c(1,1,0,0), # set outer margin (the margin around the combined plot area) - higher numbers increase the number of lines
          mar=c(5,5,5,5) # set the margin around each individual plot 
)

#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 0.6, cex = 0.1, pie = new.hap.smplloc, 
     show.mutation = 0, threshold = 0, labels(FALSE))
#add a legend to the plot
legend("topleft",colnames(new.hap.smplloc), 
       col=rainbow(ncol(new.hap.smplloc)), 
       pch=19, ncol=1, cex=0.6)

#end plot area
par(op)
# end svg file to save as
dev.off()  


#define variable
inp.f.fnm <- "Mnelei_per_year"
flnm <- c(paste("haplotype_network_",inp.f.fnm,".pdf",  sep = ""))
#paste output directory and filename together in a string
outflnm <- paste(wd00_wd05,"/",flnm,sep="")
# Exporting PFD files via postscript()           
pdf(outflnm,
    width=(1*1.0*8.2677),height=(2*1.0*2.9232))

#define the plotting area with par settings
op <- par(mfrow=c(1,1), # set number of panes inside the plot - i.e. c(2,2) would make four panes for plots
          oma=c(1,1,0,0), # set outer margin (the margin around the combined plot area) - higher numbers increase the number of lines
          mar=c(5,5,5,5) # set the margin around each individual plot 
)

#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 0.6, cex = 0.1, pie = new.hap.smplye, 
     show.mutation = 0, threshold = 0, labels(FALSE))
#add a legend to the plot
legend("topleft",colnames(new.hap.smplye), 
       col=rainbow(ncol(new.hap.smplye)), 
       pch=19, ncol=1, cex=0.6)

#end plot area
par(op)
# end svg file to save as
dev.off()  



#________________________________________________________________
# end - Plot haplotype network:
#________________________________________________________________


#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
flnm <- c(paste("haplotype_network_",inp.f.fnm,"02.pdf",  sep = ""))
#paste output directory and filename together in a string
outflnm <- paste(wd00_wd05,"/",flnm,sep="")
# Exporting PFD files via postscript()           
pdf(outflnm,
    width=(1*1.0*8.2677),height=(4*1.0*2.9232))
#define lpot arrangement
tbt.par <- par(mfrow=c(2, 1))
#_______________________________________________________________________________
#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 0.6, cex = 0.1, 
     #bg= colfpl1, 
     pie = new.hap.smplloc, 
     show.mutation = 0, threshold = 0, labels(FALSE))
#add a legend to the plot
legend("bottomleft",colnames(new.hap.smplloc), 
       #col=colfpl1,
       col=rainbow(ncol(new.hap.smplloc)), 
       pch=19, ncol=1, cex=0.5)
title(main = "a",
      cex.main = 1.4,   font.main= 2, col.main= "black",
      adj = 0.01, line = 0.1)
#_______________________________________________________________________________
#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 0.6, cex = 0.1, pie = new.hap.smplye, 
     show.mutation = 0, threshold = 0, labels(FALSE))
#add a legend to the plot
legend("bottomleft",colnames(new.hap.smplye), 
       col=rainbow(ncol(new.hap.smplye)), 
       pch=19, ncol=1, cex=0.6)
title(main = "a",
      cex.main = 1.4,   font.main= 2, col.main= "black",
      adj = 0.01, line = 0.1)
#_______________________________________________________________________________
par(tbt.par)
# end svg file to save as
dev.off()  
#reset this parameter
par(mfrow = c(1, 1)) 
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________

#________________________________________________________________
# Try and plot a neighbour join tree
#
mtr_dd_pip <- ape::dist.dna(pip)
tre_pip <- ape::njs(mtr_dd_pip)
#sort the branches in the tree
tre_pipr <- ape::ladderize(tre_pip, right = TRUE)

#https://joey711.github.io/phyloseq/plot_tree-examples.html
plot(tre_pipr, cex=0.4)

# Try and plot a neighbour join tree
#_______________________________________________________________________________
# see this website for inspiration
# https://aschuerch.github.io/posts/2017-04-24-blog-post-1
#https://bioconductor.org/packages/release/bioc/html/ggtree.html
if (!require("BiocManager", quietly = TRUE))
  if(!require(BiocManager)){
    install.packages("BiocManager")
    library(BiocManager)
  }
library(BiocManager)
if(!require(ggtree)){
  BiocManager::install("ggtree")
  library(ggtree)
}
library(textclean)
if(!require(textclean)){
  install.packages("textclean")
  library(textclean)
}

library(ggtree)
library("ggplot2")
library("ggtree")

p <- ggtree(tre_pipr) + 
  xlim(0, 0.025) + # to allow more space for labels
  geom_treescale() # adds the scale

df_tiplb01 <- as.data.frame(cbind(c(tre_pipr$tip.label)))
df_tiplb01$cat <- NA
colnames(df_tiplb01) <- c("seqNm", "cat")
df_tiplb01$cat[grepl("NCBI",df_tiplb01$seqNm)] <- "NCBI"
df_tiplb01$cat[grepl("Germany",df_tiplb01$seqNm)] <- "Germany"
df_tiplb01$cat[grepl("Jylland",df_tiplb01$seqNm)] <- "Jylland"
df_tiplb01$cat[grepl("_Fyn",df_tiplb01$seqNm)] <- "Fyn"
df_tiplb01$cat[grepl("SamsoeBallen_",df_tiplb01$seqNm)] <- "Samsoe"

tipcategories <- df_tiplb01

dd = as.data.frame(tipcategories)

p %<+% dd + 
  geom_tiplab(aes(fill = factor(cat)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
  theme(legend.position = c(0.5,0.2), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys
#_______________________________________________________________________________
#define variable
inp.f.fnm <- "Mnelei"
flnm <- c(paste("Fig03_NJ_tree_",inp.f.fnm,".pdf",  sep = ""))
#paste output directory and filename together in a string
outflnm <- paste(wd00_wd05,"/",flnm,sep="")
# # Exporting PFD files via postscript()           
# pdf(outflnm,
#     width=(1*1.0*8.2677),height=(4*1.0*2.9232))
# #plot the tree
# plot(tre_pipr, cex=0.4)
# #end plot area
# par(op)
# # end svg file to save as
# dev.off()  

#Make the dnabin a genind object
gei_pip <- adegenet::DNAbin2genind(pip)
#Make the genind object a data frame
df_geipip <- genind2df(gei_pip)
# replace any NAs with dash
df_geipip[is.na(df_geipip)] <- "-"
#make the data frame a matrix and make this a dnabin object
dnb_pip02 <- as.DNAbin(as.matrix(df_geipip))

#substitute in filename to make new output file names
nexfnm01 <- gsub("\\.fas","\\.nex",inpf01)
#write nexus format
pth_outpnxf <- paste(wd00_wd05,"/",nexfnm01,sep="")
#write the DNAbin object as a nexus file
ape::write.nexus.data(dnb_pip02,pth_outpnxf,
                      interleaved = F)


#substitute in the input filename
inpf02 <- gsub("09\\.fas","10\\.fas",inpf01)
inpf04 <- gsub("10\\.aligned","11\\.aligned",inpf01)
# paste path to directory and output filename together
pth_inpf04 <- paste(wd00_wd05,"/",inpf04,sep="")
#write the output file
ips::write.fas(dnb_pip02,pth_inpf04)

#________________________________________________________________


#________________________________________________________________
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#define input file as variable
#inpf01 <- "algn_Mnelei_18s_09.fas"
pth_inpf02
pth_inpf01 <- paste(wd00_wd01,"/",inpf01,sep="")
#Test the 'apex' package :
#read FASTA as multi.dna
mltdn_pip <- apex::read.multiFASTA(pth_inpf04)
#check what kind of object it is
#class(mltdn_pip)
#view as plotted alignment
plot(mltdn_pip, cex = 0.2)

#string split the character labels by underscore
# in the multidna object
# rowbind the nested lists to a dataframe
df_pip02 <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(mltdn_pip@labels),
                                  "_")))
#check object
#class(df_pip02)
#rename columns names
colnames(df_pip02) <- c("gnnm1","spnm1","spnm2",
                        "gnnm2","spnm3","spnm3",
                        "gnnm3")
# see the header rows of df
#head(df_pip02,7)
#subset multidna object - to exclude species that are not Mnemiopsis leidyii
mltdn_pip02 <- mltdn_pip[grepl("Mnelei",mltdn_pip@labels)]
#view as plotted alignment
plot(mltdn_pip02, cex = 0.2)

#_______________________________________________________________________________
#_______________________________________________________________________________

# ------------- Options -------------
#set fraction OK for all other files 
FracOK <- 0.90
# if the input file is a sequence file with ets2 gene fragment, then set
# the fraction higher
if (grepl("ets2",pth_inpf04)==T)
{FracOK <- 0.99}
# fraction of samples required to have data (NOT equal to "-") at position for it to be considered OK
fraction_ok <- 0.80
fraction_ok <- FracOK
max_gap_width <- 100 # 
min_grp_width <- 100 #  

# values not in the list will be marked "others"" in the plots
valuelist <- c("a","g","c","t","n","-")
#valuelist <- c("A","G","C","T","N","-")
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
# save output figures 
bSaveFigures<-F

#_______________________________________________________________________________
library(tidyverse)
library(ape)
#install.packages("ape")
if(!require(ape)){
  install.packages("ape")
  library(ape)
}
# get the working dir
prwd<- getwd()
#source("ReadFasta.R")
#source("DNA_sequence_subset.R")
# set working directory to be able to read in functions
setwd(paste0(wd00,"/suppmat03_Rcodes"))
source("Rfunction_DNA_sequence_subset.R")
source("Rfunction_ReadFasta.R")
setwd(prwd)
# ------------- Load data using ReadFasta function -------------
df <- ReadFasta(paste0(pth_inpf04))
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


filename<-paste0(substr(inpf01,1,nchar(inpf01)-4),"_cropped.csv")
folder_out <- wd00_wd05
write.table(df_out,file=paste0(folder_out,filename),row.names=F,col.names=T,sep=",",quote=F)

#pth_inpf03 <- paste0(folder_out,filename)
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

#p1

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
#p2

#____________________________________________________________
# trim the alignment
#https://www.biostars.org/p/158250/
folder <- wd00_wd05
#pth_inpf01 
# read in the alignment
fnm_toread <- pth_inpf04
# read in the alignment
al_dt01 <- ape::read.dna(fnm_toread,format="fasta", as.matrix=TRUE)
# get the upper and lower limits to cut the alignment by
uppbp <- as.numeric(gsub(".*-","",sListOK))
lowbp <- as.numeric(gsub("-.*","",sListOK))
if (uppbp==Inf){
  uppbp <- ncol(al_dt01)
lowbp <- 1}
#trim the alignment
al_dt02 <- al_dt01[,lowbp:uppbp]

outdir <- wd00_wd05
# define a file to write to
fnm01_1 <- gsub("_fastafiles\\.","_",inpf01)
fnm01_2 <- gsub("fas\\.","trimmed\\.",fnm01_1)
fnm01_2 <- gsub(".fasta.fas","\\.trimmed\\.fas",fnm01_1)
#fnm01_3 <- gsub("triplefin_out04_","triplefin_out06_",fnm01_2)
fnm02 <- gsub("\\.txt","",fnm01_2)

fnm_towrite <- paste0(outdir,"/",fnm02)
#write the trimmed alignment as a fasta file
ape::write.dna(al_dt02, file=fnm_towrite, format="fasta", nbcol=-1, colsep="")
#_______________________________________________________________________________
#_______________________________________________________________________________
#define input file as variable
#inpf01 <- "algn_Mnelei_18s_10.fas"
pth_inpf03 <- fnm_towrite
#Test the 'apex' package :
#read FASTA as multi.dna
mltdn_pip03 <- apex::read.multiFASTA(pth_inpf04)
# see the names of the sequences
mltdn_pip03@labels
# rowbind the nested lists to a dataframe
df_pip03 <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(mltdn_pip03@labels),
                                  "_")))
#check object
#class(df_pip03)
#rename columns names
colnames(df_pip03) <- c("unqnm","locality","year",
                        "yearmonth")
#make the multidna object a genind object
gi_pip03 <- apex::multidna2genind(mltdn_pip03)
#Make the genind object a data frame
df_geipip03 <- genind2df(gi_pip03)

#replace all NAs with "-"
df_geipip04 <- df_geipip03 %>% replace(is.na(.), "-")

df_geipip03 <- df_geipip04
#make the df object a matrix and make this a dnabin object
dnb_pip03 <- as.DNAbin(as.matrix(df_geipip04))
#see column names of the data frame
#colnames(df_pip03)
#colnames(df_geipip03)
#make genind and hierfstat objects by region
gi_pip_unqnm <- adegenet::DNAbin2genind(dnb_pip03,pop=df_pip03$unqnm) 
gi_pip_local <- adegenet::DNAbin2genind(dnb_pip03,pop=df_pip03$locality) 
gi_pip_years <- adegenet::DNAbin2genind(dnb_pip03,pop=df_pip03$year) 
gi_pip_yearm <- adegenet::DNAbin2genind(dnb_pip03,pop=df_pip03$yearmonth) 
#make hierfstat objects - Note that it is required to set 'pop' to the 
# vector that matches the population to which the indivudal belongs
hfst_pip_un <- hierfstat::genind2hierfstat(gi_pip_unqnm,pop=df_pip03$unqnm)
hfst_pip_lo <- hierfstat::genind2hierfstat(gi_pip_local,pop=df_pip03$locality)
hfst_pip_ye <- hierfstat::genind2hierfstat(gi_pip_years,pop=df_pip03$year)
hfst_pip_ym <- hierfstat::genind2hierfstat(gi_pip_yearm,pop=df_pip03$yearmonth)
#replace all NAs with 0
hfst_pip_un2 <- hfst_pip_un %>% replace(is.na(.), 0)
hfst_pip_lo2 <- hfst_pip_lo %>% replace(is.na(.), 0)
hfst_pip_ye2 <- hfst_pip_ye %>% replace(is.na(.), 0)
hfst_pip_ym2 <- hfst_pip_ym %>% replace(is.na(.), 0)
#overwrite the previous dataframes
hfst_pip_un <- hfst_pip_un2
hfst_pip_lo <- hfst_pip_lo2
hfst_pip_ye <- hfst_pip_ye2
hfst_pip_ym <- hfst_pip_ym2
#pairwise fst test Nei - this takes an 'hierfstat' object as input
pw_Nei_pip_un <- hierfstat::pairwise.neifst(hfst_pip_un, diploid = FALSE)
pw_Nei_pip_lo <- hierfstat::pairwise.neifst(hfst_pip_lo, diploid = FALSE)
pw_Nei_pip_ye <- hierfstat::pairwise.neifst(hfst_pip_ye, diploid = FALSE)
pw_Nei_pip_ym <- hierfstat::pairwise.neifst(hfst_pip_ym, diploid = FALSE)
# # Replace NAs with zeros 
# pw_Nei_pip_un2 <- pw_Nei_pip_un %>% replace(is.na(.), 0)
# pw_Nei_pip_lo2 <- pw_Nei_pip_lo %>% replace(is.na(.), 0)
# pw_Nei_pip_ye2 <- pw_Nei_pip_ye %>% replace(is.na(.), 0)
# pw_Nei_pip_ym2 <- pw_Nei_pip_ym %>% replace(is.na(.), 0)
# #overwrite the previous dataframes
# pw_Nei_pip_un <- pw_Nei_pip_un2
# pw_Nei_pip_lo <- pw_Nei_pip_lo2
# pw_Nei_pip_ye <- pw_Nei_pip_ye2
# pw_Nei_pip_ym <- pw_Nei_pip_ym2

# Remove rows where the pop name is unique - 
# as the 'hierfstat::pairwise.WCfst'
# function cannot handle the unique pop rows only 
# represented by only a single individual
hfst_pip_un2 <- hfst_pip_lo[hfst_pip_un$pop %in% hfst_pip_un$pop[duplicated(hfst_pip_un$pop)],]
hfst_pip_lo2 <- hfst_pip_lo[hfst_pip_lo$pop %in% hfst_pip_lo$pop[duplicated(hfst_pip_lo$pop)],]
hfst_pip_ye2 <- hfst_pip_ye[hfst_pip_ye$pop %in% hfst_pip_ye$pop[duplicated(hfst_pip_ye$pop)],]
hfst_pip_ym2 <- hfst_pip_ye[hfst_pip_ym$pop %in% hfst_pip_ym$pop[duplicated(hfst_pip_ym$pop)],]
#overwrite previous data frames
hfst_pip_un <- hfst_pip_un2
hfst_pip_lo <- hfst_pip_lo2
hfst_pip_ye <- hfst_pip_ye2
hfst_pip_ym <- hfst_pip_ym2
# The Pairwise fst test used in the heatmap colored table
# below makes use of the 'hierfstat' prepared object
#Pairwise fst test WC - this takes an 'hierfstat' object as input
pw_wc_pip_un <- hierfstat::pairwise.WCfst(hfst_pip_un, diploid = FALSE)
pw_wc_pip_lo <- hierfstat::pairwise.WCfst(hfst_pip_lo, diploid = FALSE)
pw_wc_pip_ye <- hierfstat::pairwise.WCfst(hfst_pip_ye, diploid = FALSE)
pw_wc_pip_ym <- hierfstat::pairwise.WCfst(hfst_pip_ym, diploid = FALSE)
#color the df with a heat map
#https://stackoverflow.com/questions/47733031/how-to-plot-dataframe-in-r-as-a-heatmap-grid
library(tidyverse)
# try it out on the pw_wc_pip_rg object prepared with
# the 'hierfstat' package
# make the data frame as a heat map
df_pw_wc_pip <- as.data.frame(pw_wc_pip_ye)
df_pw_wc_pip2 <- df_pw_wc_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(df_pw_wc_pip2)
#plot data frame as heat map
ggplot(df_pw_wc_pip2, aes(x = rowname,
                          y = colname, 
                          fill = value)) +
  geom_tile()

#_________________________________________________________________
# make the data frame as a heat map
df_pw_wc_pip <- as.data.frame(pw_wc_pip_lo)
df_pw_wc_pip2 <- df_pw_wc_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#head(df_pw_wc_pip2)
#plot data frame as heat map
ggplot(df_pw_wc_pip2, aes(x = rowname,
                          y = colname, 
                          fill = value)) +
  geom_tile()
#_________________________________________________________________
# make the data frame as a heat map
df_pw_wc_pip <- as.data.frame(pw_wc_pip_ym)
df_pw_wc_pip2 <- df_pw_wc_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#head(df_pw_wc_pip2)
#plot data frame as heat map
ggplot(df_pw_wc_pip2, aes(x = rowname,
                          y = colname, 
                          fill = value)) +
  geom_tile()
#_________________________________________________________________
#_________________________________________________________________
# make the data frame as a heat map with 
df_Nei_pip <- as.data.frame(pw_Nei_pip_lo)
df_Nei_pip3 <- df_Nei_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#head(df_Nei_pip3)
#plot data frame as heat map
ggplot(df_Nei_pip3, aes(x = rowname,
                        y = colname, 
                        fill = value)) +
  geom_tile()
#_________________________________________________________________
# make the data frame as a heat map
df_Nei_pip <- as.data.frame(pw_Nei_pip_ye)
df_Nei_pip3 <- df_Nei_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#head(df_Nei_pip3)
#plot data frame as heat map
ggplot(df_Nei_pip3, aes(x = rowname,
                        y = colname, 
                        fill = value)) +
  geom_tile()
#_________________________________________________________________
# make the data frame as a heat map
df_Nei_pip <- as.data.frame(pw_Nei_pip_ym)
df_Nei_pip3 <- df_Nei_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#head(df_Nei_pip3)
#plot data frame as heat map
ggplot(df_Nei_pip3, aes(x = rowname,
                        y = colname, 
                        fill = value)) +
  geom_tile()
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#_________________________________________________________________
#
# start - Make Phist table with numbers in cells
#
#_________________________________________________________________


#Pairwise Gst test Nei  - this takes a 'genind' object as input
# the leafy seadragon study: Stiller, J. et al. (2017). The Leafy Seadragon, Phycodurus eques, a Flagship Species with Low But Structured Genetic Variability. Journal of Heredity, 2017, 152–162
# Makes use of Pairwise Gst test Nei for preparing heatmap
# tables
pw_Gst_pip_ym <- mmod::pairwise_Gst_Nei(gi_pip_yearm, linearized = FALSE)
pw_Gst_pip_lo <- mmod::pairwise_Gst_Nei(gi_pip_local, linearized = FALSE)

# The Pairwise fst test used in the heatmap colored table
# below makes use of the 'hierfstat' prepared object
#Pairwise fst test WC - this takes an 'hierfstat' object as input
pw_wc_pip_lo <- hierfstat::pairwise.WCfst(hfst_pip_lo, diploid = FALSE)
pw_wc_pip_ym <- hierfstat::pairwise.WCfst(hfst_pip_ym, diploid = FALSE)

# make the data frame as a heat map
df_pw_wc_pip_lo <- as.data.frame(pw_wc_pip_lo)
df_pw_wc_pip_ym <- as.data.frame(pw_wc_pip_ym)
#replace NAs with zeroes
df_pw_wc_pip_lo[is.na(df_pw_wc_pip_lo)] = 0
df_pw_wc_pip_ym[is.na(df_pw_wc_pip_ym)] = 0
#_____________________________________________________________
# Prepare min and max values for plot
#_____________________________________________________________

df_pw_wc_pip <- df_pw_wc_pip_ym
#get the dplyr library to summarize
library(dplyr)
#get the max value per col
mx_wc_pip <- df_pw_wc_pip %>% summarise_if(is.numeric, max)
#get the min value per col
mn_wc_pip <- df_pw_wc_pip %>% summarise_if(is.numeric, min)
# get the max value among the max values obtained
mxx_pip <- max(mx_wc_pip)
mnn_pip <- max(mn_wc_pip)
#define upper and lower limits
upppl1 <- round(mxx_pip, 0)+1
lowmi1 <- round(mnn_pip, 0)-1
#rearrange dataframe to long 
df_pw_wc_pip_2 <- df_pw_wc_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
# see first header lines

#make a data frame for abbreviated sampling locations
locnm <- c("FynBogense","FynKerteminde","JyllandMariagerfjord",
           "NGermanyKielFjord","NGermanyMecklenburgerBuchtWismarBucht",
           "NJyllandLimfjord","NJyllandLimfjordLogstoer",
           "NWGermanyNSeaHelgolandRds","NWGermanyWaddenSeaBussumHaupstr",
           "SamsoeBallen","SEDenmarkMecklenburgerBucht","SjaellandSkovshoved")
locabb <- c("FBo","FKe","JMa","GKi","GMe","JLi","JLo","GHe","GWa","SBa","DMe","SSk")
# bind it together to a data frame
df_abloc <- as.data.frame(cbind(locnm,locabb))

# get the location names and abbreviated location names not already 
# included in the data frame with abbreviated location names
ncbilocab <- df_clo03$locality8[!df_clo03$locality8 %in% df_abloc$locabb]
ncbiloclo <- df_clo03$locality6[!df_clo03$locality8 %in% df_abloc$locabb]
# bind these to vectors together in a data frame
df_ncbiloc <- as.data.frame(cbind(ncbiloclo,ncbilocab))
# change the column names so that they match
colnames(df_ncbiloc) <- colnames(df_abloc)
# bind them by rows
df_abloc <- rbind(df_abloc,df_ncbiloc)
#head(df_pw_wc_pip_2)
#get the library 'pals' to be able to make use of colour gradients
library("pals")
#replace column names
colnames(df_pw_wc_pip_2) <- c("rwnm_spcs","clnm_spcs","value")
#replace with abbreviated location names
locnm2 <- df_abloc$locabb[match(df_pw_wc_pip_2$clnm_spcs, df_abloc$locnm)]
# assign back NCBI codes
locnm2[is.na(locnm2)]  <- df_pw_wc_pip_2$clnm_spcs[is.na(match(df_pw_wc_pip_2$clnm_spcs, df_abloc$locnm))]
#replace 'clnm_spcs' column
df_pw_wc_pip_2$clnm_spcs <- locnm2

#prepare plot
h<-ggplot(df_pw_wc_pip_2, aes(x = rwnm_spcs, 
                              y= reorder(clnm_spcs,
                                         desc(clnm_spcs)))) +
  geom_tile(aes(fill = value),height=1) +
  geom_text(size=4.0, aes(label = sprintf("%0.2f", value),
                          fontface = "bold"),colour="white") +
  #
  #scale_colour_manual(values=c("white")) +
  #limits is to determine within what range the colours must change
  #Values informs ggplot when
  #when colors in the color gradient is to change
  #see the different color gradients  here:
  #https://www.rdocumentation.org/packages/pals/versions/1.6/topics/kovesi
  scale_fill_gradientn(colours=pals::kovesi.linear_bgyw_15_100_c68(600), 
                       guide = "colourbar",
                       limits = c(0,1), # change limit to get diff range
                       values=c(0,1),na.value = "grey50")+
  scale_y_discrete(position = "left")+
  ggtitle("Phist_year")+
  labs(x = expression("rwnm_year"), 
       y=expression("clnm_year"))+
  
  theme(axis.text.y = element_text( size = 18,hjust=0.9),
        axis.text.x = element_text( size = 18,hjust=0.5),
        aspect.ratio =0.5,
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text( size = 18),
    #   axis.text.x.bottom = element_text( size = 14),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18))

print(h)


#_________________________________________________________________

df_pw_wc_pip <- df_pw_wc_pip_lo
#get the dplyr library to summarize
library(dplyr)
#get the max value per col
mx_wc_pip <- df_pw_wc_pip %>% summarise_if(is.numeric, max)
#get the min value per col
mn_wc_pip <- df_pw_wc_pip %>% summarise_if(is.numeric, min)
# get the max value among the max values obtained
mxx_pip <- max(mx_wc_pip)
mnn_pip <- max(mn_wc_pip)
#define upper and lower limits
upppl1 <- round(mxx_pip, 0)+1
lowmi1 <- round(mnn_pip, 0)-1

#rearrange dataframe to long 
df_pw_wc_pip_2 <- df_pw_wc_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
# see first header lines
#head(df_pw_wc_pip_2)
#get the library 'pals' to be able to make use of colour gradients
library("pals")
#replace column names
colnames(df_pw_wc_pip_2) <- c("rwnm_spcs","clnm_spcs","value")

#replace with abbreviated location names
locnm2 <- df_abloc$locabb[match(df_pw_wc_pip_2$clnm_spcs, df_abloc$locnm)]
# assign back NCBI codes
locnm2[is.na(locnm2)]  <- df_pw_wc_pip_2$clnm_spcs[is.na(match(df_pw_wc_pip_2$clnm_spcs, df_abloc$locnm))]
#replace 'clnm_spcs' column
df_pw_wc_pip_2$clnm_spcs <- locnm2

#replace with abbreviated location names
locnm2 <- df_abloc$locabb[match(df_pw_wc_pip_2$rwnm_spcs, df_abloc$locnm)]
# assign back NCBI codes
locnm2[is.na(locnm2)]  <- df_pw_wc_pip_2$rwnm_spcs[is.na(match(df_pw_wc_pip_2$rwnm_spcs, df_abloc$locnm))]
#replace 'clnm_spcs' column
df_pw_wc_pip_2$rwnm_spcs <- locnm2


#prepare plot
h2<-ggplot(df_pw_wc_pip_2, aes(x = rwnm_spcs, 
                              y= reorder(clnm_spcs,
                                         desc(clnm_spcs)))) +
  geom_tile(aes(fill = value),height=1) +
  geom_text(size=4.0, aes(label = sprintf("%0.2f", value),
                          fontface = "bold"),colour="white")  +
  #limits is to determine within what range the colours must change
  #Values informs ggplot when
  #when colors in the color gradient is to change
  #see the different color gradients  here:
  #https://www.rdocumentation.org/packages/pals/versions/1.6/topics/kovesi
  scale_fill_gradientn(colours=pals::kovesi.linear_bgyw_15_100_c68(600), 
                       guide = "colourbar",
                       limits = c(0,1), # change limit to get diff range
                       values=c(0,1),na.value = "grey50")+
  scale_y_discrete(position = "left")+
  ggtitle("Phist_locality")+
  labs(x = expression("rwnm_loc"), 
       y=expression("clnm_loc"))+
  
  theme(axis.text.y = element_text( size = 8,hjust=0.9),
        axis.text.x = element_text( angle = 90, size = 8,vjust=1,hjust=0.5),
        aspect.ratio =0.5,
        axis.title.x = element_text(angle = 90, size = 8,vjust=1,hjust=0.5),
        axis.title.y = element_text( size = 8),
        #axis.text.x.bottom = element_text( size = 14),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8))

print(h2)


#__________________________________________________________________

# Use section below to print out svg or pdf plot in a file
# # set working directory
# setwd (wd00)
# getwd()
# # Exporting EPS files via postscript()
# postscript(c(paste("plot_optimal_primer_conc_",qPCRNospec,".eps", sep = "")),
#            width=(1.6*8.2677),height=(1.5*1.6*2.9232),
#            #family = "Arial", 
#            paper = "special", onefile = FALSE,
#            horizontal = FALSE)
# 
# print(h)
# # end pdf file to save as
# dev.off()
#wd05

#defien variable for speceis
spc <- "Mnelei"
#define output filename
outpfl01 <- paste("heatmap_Phist_",spc,"smplyear_01.pdf",  sep = "")
outpfl02 <- paste("heatmap_Phist_",spc,"smplloca_02.pdf",  sep = "")
# paste output filename and path to outputdirectory together
pth_outpfl01 <- paste(wd00_wd05,"/",outpfl01,sep="")
# print the plot as pdf
pdf(c(pth_outpfl01)
    ,width=(1.6*8.2677),height=(1.6*2.9232*2))
print(h)
dev.off()

# paste output filename and path to outputdirectory together
pth_outpfl02 <- paste(wd00_wd05,"/",outpfl02,sep="")
# print the plot as pdf
pdf(c(pth_outpfl02)
    ,width=(1.6*8.2677),height=(1.6*2.9232*2))
print(h2)
dev.off()

#h
#h2


#_______________________________________________________________________________
# Make table 02 and table 03
#_______________________________________________________________________________
#copy plots
h01 <- h
h02 <- h2
#change labels on axis
h01 <- h01 + ylab("sample year") + xlab("sample year") 
h02 <- h02 + ylab("sample location") + xlab("sample location") 
#change labels on legend
h01 <- h01 + labs(fill='PhiST')
h02 <- h02 + labs(fill='PhiST')
# Add titles
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#h01t <- h01 + labs(title = "A")#,
h01t <- h01 + labs(title = "")#,
# Add titles
#h02t <- h02 + labs(title = "B")#,
h02t <- h02 + labs(title = "")#,
# ------------- plot Combined figure -------------
library(patchwork)
# set a variable to TRUE to determine whether to save figures
bSaveFigures <- T
#see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
# on how to arrange plots in patchwork
ph <-   h01t +
  #h02t +
  
  plot_layout(nrow=1,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") +
  plot_annotation(caption=inpf01) #& theme(legend.position = "bottom")
#p
#make filename to save plot to
figname01 <- paste0("Table03_PhiST_year_02",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(ph,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}
#_______________________________________________________________________________
#_______________________________________________________________________________
ph <-   #h01t +
  h02t +
  
  plot_layout(nrow=1,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") +
  plot_annotation(caption=inpf01) #& theme(legend.position = "bottom")
#p
#make filename to save plot to
figname01 <- paste0("Table02_PhiST_locality_02",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(ph,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}
#_______________________________________________________________________________
#_________________________________________________________________
#
# end - Make Phist table with numbers in cells 
#
#_________________________________________________________________


#
#_________________________________________________________________
# start - plot pies on map with ggplot
#_________________________________________________________________
# #transpose the haplotype data frame prepared earlier
#this makes the data frame a table
tbl_hap_loc01 <- t(new.hap.smplloc)
# see the class of the object
#class(tbl_hap_loc01)
#make it a data frame instead
df_hap_loc01 <- as.data.frame(tbl_hap_loc01)
#check the column names
#colnames(df_hap_loc01)
#df_hap_loc01$smplloca
# reshape the data frame for long to wide
df_hap_loc02 <- reshape(data=df_hap_loc01,idvar="smplloca",
                          v.names = "Freq",
                          timevar = "Var2",
                          direction="wide")
# #define input file
# fl1 <- "collect_loc_Mnelei_smpls01.csv"
# #paste path and filename together
# pth_fl01 <- paste(wd00_wd01,"/",fl1,sep="")
# #read in a csv file with positions for sampling locations
# df_clo <- as.data.frame(read_delim(pth_fl01,delim = ";"))
# 
# df_clo$locality4 <-  df_clo$locality
# #Modify location names
# df_clo$locality4 <- gsub("FynBogense","Funen, Bogense", df_clo$locality4 )
# df_clo$locality4 <- gsub("FynKerteminde","Funen, Kerteminde", df_clo$locality4 )
# df_clo$locality4 <- gsub("JyllandMariagerfjord","Jutland, Mariager Fjord", df_clo$locality4 )
# df_clo$locality4 <- gsub("NGermanyKielFjord","Germany, Kiel Fjord", df_clo$locality4 )
# df_clo$locality4 <- gsub("NGermanyMecklenburgerBuchtWismarBucht","Germany, Wismar Bight", df_clo$locality4 )
# df_clo$locality4 <- gsub("NJyllandLimfjord","Jutland, Limfjord", df_clo$locality4 )
# df_clo$locality4 <- gsub("NJyllandLimfjordLogstoer","Jutland, Limfjord" , df_clo$locality4 )
# df_clo$locality4 <- gsub("NWGermanyNSeaHelgolandRds", "North Sea, Helgoland Roads", df_clo$locality4 )
# df_clo$locality4 <- gsub("NWGermanyWaddenSeaBussumHaupstr","Germany, Büsum", df_clo$locality4 )
# df_clo$locality4 <- gsub("SamsoeBallen","Samsøe, Ballen", df_clo$locality4 )
# #df_clo$locality4 <- gsub("SEDenmarkMecklenburgerBucht","", df_clo$locality4 )
# df_clo$locality4 <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", df_clo$locality4 )
# 
# # df_clo
# #use biogeo dms2dd funciton to convert to decimal degrees
# df_clo$dec_lat <- biogeo::dms2dd(df_clo$lat_deg,
#                    df_clo$lat_min,
#                    df_clo$lat_sec,
#                    df_clo$lat_sph)
# #use biogeo dms2dd function to convert to decimal degrees
# df_clo$dec_lon <- biogeo::dms2dd(df_clo$lon_deg,
#                                  df_clo$lon_min,
#                                  df_clo$lon_sec,
#                                  df_clo$lon_sph)
# #exclude the location: 'SEDenmarkMecklenburgerBucht'
# # as it is the same as "NGermanyMecklenburgerBuchtWismarBucht", and only
# #one of these are needed
# # same for the location 'NJyllandLimfjordLogstoer' which is identical to NJyllandLimfjord
# df_clo <- df_clo[!df_clo$locality=="SEDenmarkMecklenburgerBucht",]
# df_clo <- df_clo[!df_clo$locality=="NJyllandLimfjordLogstoer",]
# 
# lN02loc <- unique(df_lN02$location)
# lN02lat <- df_lN02$declat[match(lN02loc,df_lN02$location)]
# lN02lon <- df_lN02$declon[match(lN02loc,df_lN02$location)]
# df_locN02 <- as.data.frame(cbind(lN02loc,lN02lat,lN02lon))
# # get number of locations for Genbank records
# nlocgb <- nrow(df_locN02)
# df_clo03 <- df_clo
# nclloc <- ncol(df_clo)
# df_clo03 [ nrow(df_clo03) + nlocgb , ] <- NA
# df_clo03$locality4[is.na(df_clo03$locality4)] <- df_locN02$lN02loc
# df_clo03$dec_lat[is.na(df_clo03$locality2)] <- df_locN02$lN02lat
# df_clo03$dec_lon[is.na(df_clo03$locality2)] <- df_locN02$lN02lon
# #replace non ascii characteres
# df_clo03$locality5 <- textclean::replace_non_ascii(df_clo03$locality4)
# #add a column with sampling locations to be able to 
# # match the localities
# df_clo03$locality2 <- abbreviate(df_clo03$locality5, 1, named = FALSE)
# 
# df_clo03 <- df_clo03[!grepl("NCBI",df_clo03$locality4),]
# #match between data frames
 df_hap_loc02$dec_loc2 <- df_clo03$locality2[match(df_hap_loc02$smplloca,df_clo03$locality4)]

# sum for duplicated values in row
# this is to add up the multiple entries row for the same 
# localities
# https://stackoverflow.com/questions/10180132/consolidate-duplicate-rows
library(plyr)
df_hap_loc03 <- plyr::ddply(df_hap_loc02,"dec_loc2",numcolwise(sum))
df_hap_loc03 <- df_hap_loc03[!is.na(df_hap_loc03[,1]),]
df_hap_loc02 <-   df_hap_loc03

#match between data frames
df_hap_loc02$dec_lat <- df_clo03$dec_lat[match(df_hap_loc02$dec_loc2,df_clo03$locality2)]
df_hap_loc02$dec_lon <- df_clo03$dec_lon[match(df_hap_loc02$dec_loc2,df_clo03$locality2)]

#df_hap_loc02$dec_loc2
# limit the data frame to remove any rows that have NAs
df_hap_loc03 <-  df_hap_loc02[complete.cases(df_hap_loc02), ] 
#make a c
df_hap_loc03 <- df_hap_loc02
#modify the colnames
colnames(df_hap_loc03) <- c(gsub("Freq\\.","",colnames(df_hap_loc02)))



# https://towardsdatascience.com/using-ggplot-to-plot-pie-charts-on-a-geographical-map-bb54d22d6e13
#count the number of columns, and subtract 2
enc1 <- ncol(df_hap_loc03)-2 
#head(df_hap_loc03,3)
#sum up for each row
df_hap_loc03$rws <- rowSums(df_hap_loc03[,2:enc1])
rws <- as.numeric(df_hap_loc03$rws)
#count the columns and subtract 3
enc <- ncol(df_hap_loc03)-3

#install the package that allows for making pit charts in ggplot
if(!require("scatterpie")){
  install.packages("scatterpie", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("scatterpie")
}
library(scatterpie)

#make a viridis colour range
cl03 <- pals::viridis(length(unique(df_hap_loc03[,c(2:enc)])))
cl03 <- pals::inferno(length(unique(df_hap_loc03[,c(2:enc)])))
ncolHptloc03 <- length(unique(colnames(df_hap_loc03)[2:enc]))
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
colfunc <- colorRampPalette(cbbPalette2)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#
cl03 <- colfunc(ncolHptloc03)
# https://towardsdatascience.com/using-ggplot-to-plot-pie-charts-on-a-geographical-map-bb54d22d6e13
# Using map_data()
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
if (!exists("world"))
{  
  #world <- ne_countries(scale = 10, returnclass = "sf")
  world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
}

# #____
# world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")

df_hap_loc03$rws2 <- df_hap_loc03$rws
#replace NAs with zeros
df_hap_loc04 <- df_hap_loc03[!is.na(df_hap_loc03$rws),]
df_hap_loc05 <- df_hap_loc04
df_hap_loc04$dec_lat <- as.numeric(df_hap_loc04$dec_lat)
df_hap_loc04$dec_lon <- as.numeric(df_hap_loc04$dec_lon)
#https://guangchuangyu.github.io/2016/12/scatterpie-for-plotting-pies-on-ggplot/
world <- ggplot2::map_data('world')
jitlvl <- 0.017
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p01 <- ggplot(data = world) +
  geom_map(map=world, aes(map_id=region), fill="grey",
           color="black") +
  #geom_sf(color = "black", fill = "azure3") +
  #https://ggplot2.tidyverse.org/reference/position_jitter.html
  # use 'geom_jitter' instead of 'geom_point' 
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.60), 
                              data = df_hap_loc04, 
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  
  scale_color_manual(values=c(rep("black",
     length(unique(df_hap_loc04[,c(2:enc)]))))) +

  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
#scale_colour_viridis_d(breaks=11) +
#scale_color_brewer(palette="Dark2") +
#https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  
geom_scatterpie_legend(df_hap_loc04$rws*0.60, x=-40, y=30) +
# set alpha values for color intensity of fill color in point
#https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
#define limits of the plot
ggplot2::coord_sf(xlim = c(-90, 70),
                  ylim = c(0, 60.0),
                  expand = FALSE)
# change labels on axis
p01 <- p01 + xlab("Longitude") + ylab("Latitude")
#dev.off()
# change label for legend
p01 <- p01 + labs(fill='Haplotype')
#https://www.statology.org/ggplot-background-color/
#p01 <- p01 + theme_minimal() #no background annotations
# #https://www.statology.org/ggplot-background-color/
p01 <- p01 + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                   panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,
# see the plot
p01

df_hap_loc03$dec_lat<- as.numeric(df_hap_loc03$dec_lat)
df_hap_loc03$dec_lon<- as.numeric(df_hap_loc03$dec_lon)
#make a viridis colour range
#cl03 <- pals::viridis(length(unique(df_hap_loc03[,c(2:enc)])))
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p02 <- ggplot(data = world) +
  geom_map(map=world, aes(map_id=region), fill="grey",
           color="black") +
  
  #  geom_sf(color = "black", fill = "azure3") +
  #https://ggplot2.tidyverse.org/reference/position_jitter.html
  # use 'geom_jitter' instead of 'geom_point' 
  
  
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.04), 
                              data = df_hap_loc03, 
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  
  
  scale_color_manual(values=c(rep("black",
                length(unique(df_hap_loc04[,c(2:enc)]))))) +
  
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
#scale_colour_viridis_d(breaks=11) +
#scale_color_brewer(palette="Dark2") +
geom_scatterpie_legend(df_hap_loc04$rws*0.04, x=14, y=57) +
  
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
# set alpha values for color intensity of fill color in point
#https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
#define limits of the plot
ggplot2::coord_sf(xlim = c(6, 15.4),
                  ylim = c(53.8, 58.0),
                  expand = FALSE)
#https://www.statology.org/ggplot-background-color/
#p02 <- p02 + theme_minimal() #no background annotations
# #https://www.statology.org/ggplot-background-color/
p02 <- p02 + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                   panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,

# see the plot
#p02
#dev.off()
#change labels on axis
p02 <- p02 + xlab("Longitude") + ylab("Latitude")
#change labels on legend
p02 <- p02 + labs(fill='Haplotype')
p02
# Add titles
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p01t <- p01 + labs(title = "a")#,
# Add titles
# p02t <- p02 + labs(title = "eDNA samples attempted",
#                    subtitle = "at least approv controls and 1 or 2 pos repl")#,
p02t <- p02 + labs(title = "b")#,


# ------------- plot Combined figure -------------
library(patchwork)
# set a variable to TRUE to determine whether to save figures
bSaveFigures <- T
#see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
# on how to arrange plots in patchwork
p <-  p01t +
  p02t +
  
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") +
  plot_annotation(caption=inpf01) #& theme(legend.position = "bottom")
#p
#make filename to save plot to
figname01 <- paste0("map_haplotype_pie_diagr",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}

#make filename to save plot to
figname01 <- paste0("Fig04_map_haplotype_pie_diagr",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}

#_________________________________________________________________
# end - plot pies on map with ggplot
#_________________________________________________________________

#_________________________________________________________________
#_________________________________________________________________
# plot collection points
#_________________________________________________________________

if(!require(sf)){
  install.packages("sf")

}
library(sf)
#install package if needed
if(!require(rnaturalearth)){
  install.packages("rnaturalearth")
  #install.packages("rnaturalearth", repos = "http://packages.ropensci.org", type = "source")
  
}
library(rnaturalearth)
#install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
# # 
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
 if (!exists("denm_map"))
 {
  denm_map <- rnaturalearth::ne_countries(country=c("denmark","sweden","germany",
                                     "norway", "poland"),scale = 10, returnclass = "sf")
  #wor_map <- ne_countries(country="world",scale = 10, returnclass = "sf")
}


#df_clo$locality
jitlvl <- 0.24
jitlvl <- 0.024
dev.off()
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
colfunc <- colorRampPalette(cbbPalette2)

cl <- cbbPalette
# subset to exclude NCBI sequences
df_clo2 <- df_clo[!grepl("NCBI",df_clo$locality),]
df_clo2$locality3 <- df_clo2$locality
#Modify location names
df_clo2$locality <- gsub("FynBogense","Funen, Bogense", df_clo2$locality )
df_clo2$locality <- gsub("FynKerteminde","Funen, Kerteminde", df_clo2$locality )
df_clo2$locality <- gsub("JyllandMariagerfjord","Jutland, Mariager Fjord", df_clo2$locality )
df_clo2$locality <- gsub("NGermanyKielFjord","Germany, Kiel Fjord", df_clo2$locality )
df_clo2$locality <- gsub("NGermanyMecklenburgerBuchtWismarBucht","Germany, Wismar Bight", df_clo2$locality )
df_clo2$locality <- gsub("NJyllandLimfjord","Jutland, Limfjord", df_clo2$locality )
#df_clo2$locality <- gsub("NJyllandLimfjordLogstoer","Jutland, Limfjord-Løgstør" , df_clo2$locality )
df_clo2$locality <- gsub("NWGermanyNSeaHelgolandRds", "North Sea, Helgoland Roads", df_clo2$locality )
df_clo2$locality <- gsub("NWGermanyWaddenSeaBussumHaupstr","Germany, Büsum", df_clo2$locality )
df_clo2$locality <- gsub("SamsoeBallen","Samsøe, Ballen", df_clo2$locality )
#df_clo2$locality <- gsub("SEDenmarkMecklenburgerBucht","", df_clo2$locality )
df_clo2$locality <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", df_clo2$locality )

#identify unique localities
nloc <- length(unique(df_clo2$locality))
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#
cl <- colfunc(nloc)
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p04 <- 
  ggplot(data = denm_map) +
  #ggplot(st_transform(norw_map, 9122)) +
  geom_sf(color = "black", fill = "azure3") +
  geom_jitter(data = df_clo2, 
              aes(x = dec_lon, y = dec_lat, #, 
                  color=locality,
                  fill=locality,
                  shape=locality),
              width = jitlvl, #0.07, jitter width 
              height = jitlvl, #0.07, # jitter height
              size = 6) +
  #manually set the pch shape of the points
  scale_shape_manual(values=c(rep(21,nloc))) +
  #set the color of the points
  #here it is black, and repeated the number of times
  #matching the number of species listed
  scale_color_manual(values=c(rep("black",nloc))) +
  #set the color of the points
  #use alpha to scale the intensity of the color
  scale_fill_manual(values=alpha(
    c(cl),
    c(0.7)
  ))+
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(4, 16.4),
                    ylim = c(53.4, 58.4),
                    expand = FALSE)
# change label for legend - Notice that you need to change for all 3 variables
# you called 'aes' in 'geom_jitter'
p04 <- p04 + labs(fill='Location')
p04 <- p04 + labs(color='Location')
p04 <- p04 + labs(shape='Location')
p04 <- p04 + xlab("Longitude") + ylab("Latitude")
# #https://www.statology.org/ggplot-background-color/
p04 <- p04 + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                     panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,
# add border around plot
p04 <- p04 + theme(panel.border = element_rect(color = "white",
                                               fill = NA,
                                               size = 1.0))
# change background of legend
p04 <- p04 + theme(legend.key = element_rect(fill = "white"))
#p04
#make filename to save plot to
figname06 <- paste0("map_samples_",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname06,sep="")
if(bSaveFigures==T){
  ggsave(p04,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}


#make filename to save plot to
figname06 <- paste0("Fig01_map_samples_",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname06,sep="")
if(bSaveFigures==T){
  ggsave(p04,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}

# # 
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
world2 <- ne_countries(scale = 10, returnclass = "sf")
#df_hap_loc04[(df_hap_loc04$dec_loc2=="M"),]
#cl03
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p05 <- ggplot(data = world2) +
  geom_sf(color = "black", fill = "azure3") +
  #geom_sf(color = "black", fill = "azure3") +
  #https://ggplot2.tidyverse.org/reference/position_jitter.html
  # use 'geom_jitter' instead of 'geom_point' 
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.20), 
                              data = df_hap_loc04, 
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  
  scale_color_manual(values=c(rep("black",
                                  length(unique(df_hap_loc04[,c(2:enc)]))))) +
  
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
  #scale_colour_viridis_d(breaks=11) +
  #scale_color_brewer(palette="Dark2") +
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  
  geom_scatterpie_legend(df_hap_loc04$rws*0.20, x=-10, y=47) +
  # set alpha values for color intensity of fill color in point
  #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(-12, 16),
                    ylim = c(34.8, 58.0),
                    expand = FALSE) +
  theme(aspect.ratio=7/8)
# copy teh data frame to use later on
df_hap_loc07 <- df_hap_loc04
# change labels on axis
p05 <- p05 + xlab("Longitude") + ylab("Latitude")
#dev.off()
# change label for legend
p05 <- p05 + labs(fill='haplotype')
# see the plot
p05
#make a viridis colour range
#cl03 <- pals::viridis(length(unique(df_hap_loc03[,c(2:enc)])))
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p06 <-  ggplot(data = world2) +
  geom_sf(color = "black", fill = "azure3") +
  
  #  geom_sf(color = "black", fill = "azure3") +
  #https://ggplot2.tidyverse.org/reference/position_jitter.html
  # use 'geom_jitter' instead of 'geom_point' 
  
  
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.04), 
                              data = df_hap_loc03, 
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  
  
  scale_color_manual(values=c(rep("black",
                                  length(unique(df_hap_loc04[,c(2:enc)]))))) +
  
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
  #scale_colour_viridis_d(breaks=11) +
  #scale_color_brewer(palette="Dark2") +
  geom_scatterpie_legend(df_hap_loc04$rws*0.04, x=14, y=57) +
  
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  # set alpha values for color intensity of fill color in point
  #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(4, 17.4),
                    ylim = c(53.8, 58.0),
                    expand = FALSE) +
  theme(aspect.ratio=2/7)
# see the plot
p06
#dev.off()
#change labels on axis
p06 <- p06 + xlab("Longitude") + ylab("Latitude")
#change labels on legend
p06 <- p06 + labs(fill='haplotype')

#https://www.statology.org/ggplot-background-color/
p06 <- p06 + theme(panel.background = element_rect(fill = 'white', color = 'white')) #,
# panel.grid.major = element_line(color = 'red' , linetype = 'dotted'),
# panel.grid.minor = element_line(color = 'green', size = 2))
p06
# Add titles
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p05t <- p05 + labs(title = "a")#,
# Add titles
# p06t <- p06 + labs(title = "eDNA samples attempted",
#                    subtitle = "at least approv controls and 1 or 2 pos repl")#,
p06t <- p06 + labs(title = "b")#,
# ------------- plot Combined figure -------------
library(patchwork)
# set a variable to TRUE to determine whether to save figures
bSaveFigures <- T
#see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
# on how to arrange plots in patchwork
p <-  p05t +
  p06t +
  
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") +
  plot_annotation(caption=inpf01) #& theme(legend.position = "bottom")
#p
#make filename to save plot to
figname01 <- paste0("Fig04_v03_map_haplotype_pie_diagr_02",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}
#_________________________________________________________________
# end - plot pies on map with ggplot
#_________________________________________________________________
csvf <- paste0(wd00_wd05,"/df_haplo03.csv")
write.csv(df_hap_loc03, file = csvf)
#



# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
cbbPalette1 <- c("firebrick4","firebrick2",#"orange",
                 "gold3")
colfunc2 <- colorRampPalette(cbbPalette2)
colfunc1 <- colorRampPalette(cbbPalette1)
cl <- cbbPalette
# subset to include NCBI sequences
cloc1 <-  unique(colnames(new.hap.smplloc)[!grepl(", ",colnames(new.hap.smplloc))])
lcloc1 <- length(cloc1)
# subset to exclude NCBI sequences
cloc2 <-  unique(colnames(new.hap.smplloc)[!grepl("NCBI",colnames(new.hap.smplloc))])
lcloc2 <- length(cloc2)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#
#for NCBI categories
colfpl1 <- colfunc1(lcloc1)
#for non-NCBI categories
colfpl2 <- colfunc2(lcloc2)
# bind to dataframes
df_Cl1 <- as.data.frame(cbind(cloc1,colfpl1))
df_Cl2 <- as.data.frame(cbind(cloc2,colfpl2))
colnames(df_Cl1) <- c("coll_loc","colfcol_loc")
colnames(df_Cl2) <- c("coll_loc","colfcol_loc")
df_cll <- rbind(df_Cl1,df_Cl2)
#match to get same order of colors as for collection localities
colfh <- df_cll$colfcol_loc[match((colnames(new.hap.smplloc)),df_cll$coll_loc)]
# uncomment to see colours
#plot(rep(1,20),col=colfunc(20),pch=19,cex=3) #check colors for species

NCBIsmpl <- colnames(new.hap.smplye)[grepl("NCBI",colnames(new.hap.smplye))]
yearsmpl <- colnames(new.hap.smplye)[!grepl("NCBI",colnames(new.hap.smplye))]
clfy1 <- df_cll$colfcol_loc[match(NCBIsmpl,df_cll$coll_loc)]

cbbPalette3 <- c("black","cyan","white")
cbbPalette3 <- c("black","firebrick2","brown","blue","cyan","white")
colfunc3 <- colorRampPalette(cbbPalette3)
clfy2 <- colfunc3(length(yearsmpl))

df_cfy3 <- as.data.frame(cbind(yearsmpl,clfy2))
df_cfy4 <- as.data.frame(cbind(NCBIsmpl,clfy1))
colnames(df_cfy3) <- c("smplnm","colour")
colnames(df_cfy4) <- c("smplnm","colour")
df_cfy5 <- as.data.frame(rbind(df_cfy3,df_cfy4))

clrf6 <- df_cfy5$colour[match(colnames(new.hap.smplye),df_cfy5$smplnm)]
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
flnm <- c(paste("Fig02_haplotype_network_",inp.f.fnm,"02.pdf",  sep = ""))
flnm <- c(paste("Fig02_haplotype_network_",inp.f.fnm,"02.jpg",  sep = ""))
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
#_______________________________________________________________________________
#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 0.6, 
     cex = 0.1, # set size of roman numerals on circles for haplotype ID
     bg= colfh, 
     pie = new.hap.smplloc, 
     show.mutation = 2, threshold = 0, labels(T))
#add a legend to the plot
legend("bottomleft",colnames(new.hap.smplloc), 
       pt.bg=colfh,box.col=NA,
       #col=rainbow(ncol(new.hap.smplloc)), 
       pt.lwd=0.4,
       pch=21, ncol=1, cex=0.8)
title(main = "a",
      cex.main = 1.8,   font.main= 2, col.main= "black",
      adj = 0.01, line = 0.1)
#_______________________________________________________________________________
#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 0.6, 
     cex = 0.1, # set size of roman numerals on circles for haplotype ID
     bg=clrf6,
     pie = new.hap.smplye, 
     show.mutation = 2, threshold = 0, labels(TRUE))
#add a legend to the plot
legend("bottomleft",colnames(new.hap.smplye), 
       #col=rainbow(ncol(new.hap.smplye)), 
       pt.bg=clrf6,box.col=NA,pt.lwd=0.4,
       pch=21, ncol=1, cex=0.8)
title(main = "b",
      cex.main = 1.8,   font.main= 2, col.main= "black",
      adj = 0.01, line = 0.1)
#_______________________________________________________________________________
par(tbt.par)
# end svg file to save as
dev.off()  
#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________

#______________________________________________________________________________
# Make table for samples collected - start
#______________________________________________________________________________
# Get rownames from dataframe
rwnm01 <- rownames(df_pip)
# split by delimeter and make a new data frame
df_r02 <- data.frame(do.call('rbind', strsplit(as.character(rwnm01),'_',fixed=TRUE)))
#change column names on data frame
colnames(df_r02) <- c("smplnm","collloc","collyear","collmnth")
#replace in location names
df_clo2$locality <- gsub(",","",df_clo2$locality)
#df_clo03$locality
# get sampling position
df_r02$lat_deg <- df_clo2$lat_deg[match(df_r02$collloc,df_clo2$locality)]
df_r02$lat_min <- df_clo2$lat_min[match(df_r02$collloc,df_clo2$locality)]
df_r02$lat_sec <- df_clo2$lat_sec[match(df_r02$collloc,df_clo2$locality)]
df_r02$lat_sph <- df_clo2$lat_sph[match(df_r02$collloc,df_clo2$locality)]
df_r02$lon_deg <- df_clo2$lon_deg[match(df_r02$collloc,df_clo2$locality)]
df_r02$lon_min <- df_clo2$lon_min[match(df_r02$collloc,df_clo2$locality)]
df_r02$lon_sec <- df_clo2$lon_sec[match(df_r02$collloc,df_clo2$locality)]
df_r02$lon_sph <- df_clo2$lon_sph[match(df_r02$collloc,df_clo2$locality)]
# make an empty column for accession numbers
df_r02$NCBIaccsNo <- NA
#replace in this column with accession numbers
df_r02$NCBIaccsNo[grepl("Mnelei[A-Z]{2}",df_r02$smplnm)] <- gsub("Mnelei","",df_r02$smplnm[grepl("Mnelei[A-Z]{2}",df_r02$smplnm)])
# make an empty column for species names
df_r02$spcNm <- NA
df_r02$spcNm[grepl("Mnelei",df_r02$smplnm)] <- "Mnemiopsis leidyi"
#replace the non Mnemiopsis leidyi samples
notMnelei <- df_r02$smplnm[!grepl("Mnelei",df_r02$smplnm)]
notMnelei[grepl("^.*[A-Z]{2}.*$",notMnelei)] <- gsub("^(.*)([A-Z]{2}.*)$","\\1_\\2",notMnelei[grepl("^.*[A-Z]{2}.*$",notMnelei)])
notMnelei[grepl("^.*\\.[A-Z]{1}.*$",notMnelei)] <- gsub("\\.","_",notMnelei[grepl("^.*\\.[A-Z]{1}.*$",notMnelei)])
# split by delimeter and make a new data frame
df_nMl <- data.frame(do.call('rbind', strsplit(as.character(notMnelei),'_',fixed=TRUE)))
colnames(df_nMl) <- c("spcNm","NCBIaccsNos")
#get non Mnemiopsis samples and modify
df_r02$spcNm[!grepl("Mnelei",df_r02$smplnm)] <- df_nMl$spcNm
df_r02$NCBIaccsNo[!grepl("Mnelei",df_r02$smplnm)] <- df_nMl$NCBIaccsNos
#replace accession numbers
df_r02$NCBIaccsNosmplnm <- df_r02$NCBIaccsNo
df_r02$NCBIaccsNosmplnm[is.na(df_r02$NCBIaccsNosmplnm)] <- df_r02$smplnm[is.na(df_r02$NCBIaccsNosmplnm)]
#replace species names
df_r02$spcNm <- gsub("Bolinopsissp","Bolinopsis sp", df_r02$spcNm)
df_r02$spcNm <- gsub("maculata"," maculata",df_r02$spcNm)
df_r02$spcNm <- gsub("crystallina"," crystallina",df_r02$spcNm)
#paste collection position in to one string
df_r02$posloc <- paste(df_r02$lat_deg,"º",df_r02$lat_min,"'",df_r02$lat_sec,"''",df_r02$lat_sph,
                       ";",
                       df_r02$lon_deg,"º",df_r02$lon_min,"'",df_r02$lon_sec,"''",df_r02$lon_sph,
                       sep="")
#replace NA containing positions
df_r02$posloc[grepl("NA",df_r02$posloc)] <- "NA"
#copy columns
df_r02$collmnth2 <- df_r02$collmnth
df_r02$collloc2 <- df_r02$collloc
# replace with NAs
df_r02$collloc2[grepl("NCBI",df_r02$collloc2)] <- "NA"
df_r02$collmnth2[grepl("NCBI",df_r02$collmnth2)] <- "NA"
df_r02$collmnth2[grepl("^[A-Z]{1}[0-9]{+}",df_r02$collmnth2)] <- "NA"
df_r02$collloc2[grepl("crystal",df_r02$collloc2)] <- "NA"
df_r02$NCBIaccsNosmplnm[grepl("crystal",df_r02$NCBIaccsNosmplnm)] <- "NA"
df_r02$collloc2[grepl("Bolinop",df_r02$collloc2)] <- "NA"

#match to get abbreviated location name
df_r02$collloc3 <- df_abloc$locabb[match(df_r02$collloc2,df_abloc$locnm)]

df_r02$collloc3[df_r02$collloc3=="GMe"] <- "GWB"
df_r02$collloc3[df_r02$collloc3=="GWa"] <- "GWB"
# round decimal degrees and paste together in a string
df_clo03$dec_lat2 <- round(as.numeric(df_clo03$dec_lat),3)
df_clo03$dec_lon2 <-  round(as.numeric(df_clo03$dec_lon),3)
df_clo03$dec_latlon2 <- paste0(df_clo03$dec_lat2,"; ",df_clo03$dec_lon2)
#match to get sampling location 
df_r02$posloc <- df_clo03$dec_latlon2[match(df_r02$collloc3,df_clo03$locality8)]
#substitute to get accession number
df_r02$accnmb <- gsub("Mnelei","",df_r02$smplnm)
df_r02$collyear2 <- df_lN02$smplyear2[match(df_r02$accnmb,df_lN02$accession_nmb)]
# add years and month for local DK Mnelei samples 
df_r02$collyear2[!grepl("NCBI",df_r02$collmnth)] <- df_r02$collmnth[!grepl("NCBI",df_r02$collmnth)]

# define columns to keep
keeps <- c("spcNm",
           "NCBIaccsNosmplnm",
           "posloc",
           "collyear2",
           "collloc3" )
#keep defined columns
df_r03 <-  df_r02[keeps]
#rename columns
colnames(df_r03) <- c("Genus species",
                      "sample Number or NCBI Accession number",
                      "Latitude Longitude for sample collection",
                      "Year and month collected",
                      "Locality abbr.")

#https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
if(!require(tableHTML)){
  install.packages("tableHTML")
  library(tableHTML)
}
library(tableHTML)
require(tableHTML)
#try the tableHTML with no border
tableHTMLt03 <- df_r03 %>% 
  tableHTML(border = 0) 
#paste path and file name together
pth_fl02 <- paste(wd00_wd05,"/","Table01a_samples_collected.html",sep="")
pth_fl03 <- paste(wd00_wd05,"/","Table01b_samples_collected.csv",sep="")
#and to export in a file a html file
write_tableHTML(tableHTMLt03, file = pth_fl02)
# and to a csv file
write.csv(df_r03,file=pth_fl03)


# Sort by vector name [z] then [x]
df_r03 <- df_r03[
  with(df_r03, order(df_r03$`Genus species`, df_r03$`sample Number or NCBI Accession number`)),]

row.names(df_r03) <- NULL
# get package
if(!require("kableExtra")){
  # see https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html#Table_Styles
  #For dev version
  install.packages("devtools")
  devtools::install_github("haozhu233/kableExtra")
  
  library("kableExtra")
}
# get unique location abbrevaitions, and match to get location names
abbr8 <- unique(df_clo03$locality8)
lcnm8 <- df_clo03$locality4[match(abbr8, df_clo03$locality8)]
# substitute in names of locations to remove colons and missing spaces
lcnm8 <- gsub(":",", ",lcnm8)
lcnm8 <- gsub("CaspianSea","Caspian Sea",lcnm8)
lcnm8 <- gsub("BalticSea","Baltic Sea",lcnm8)
lcnm8 <- gsub("CentralWAtlantic","Central W Atlantic",lcnm8)
lcnm8 <- gsub("NEAtlantic","NE Atlantic",lcnm8)
lcnm8 <- gsub("NWAtlantic","NW Atlantic",lcnm8)
lcnm8 <- gsub("AtlanticOcean","Atlantic Ocean",lcnm8)
#paste long nmae and abbreviatoin together to use for table header
ablo8 <- paste0(lcnm8," (",abbr8,")")
# paste the vector together to be one single string
ablo8 <- paste(ablo8,collapse="; ")
#paste together a text to use as header in the table
Tbltxt01 <- paste0("Table 1. List samples collected from warty comb jelly (Mnemiopsis leidyi) in the seas around Denmark together with sequences obtained from genetic databases. Sampled locations are abbreviated: ",
                   ablo8,                   
". Additional sequences was obtained from Mnemiopsis and other ctenophore species from NBCI GenBank database as indicated by accession numbers."
)
df_r04 <- df_r03 %>%
  kableExtra::kbl(caption = Tbltxt01) %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria") %>%
  kableExtra::kable_styling(latex_options = c("striped"))  %>%
  column_spec(1, italic = T) 
#and to export in a file a html file
kableExtra::save_kable(df_r04,file=pth_fl02)

# Get unique sampling times
#unique(df_r03$`Year and month collected`)



#______________________________________________________________________________
# Make table for samples collected - end
#______________________________________________________________________________

#______________________________________________________________________________
# Make a second map with pies for haplotypes -start
# Add jitter to the sample localities to see all sampling sites
#______________________________________________________________________________
df_hap_loc01 <- as.data.frame(tbl_hap_loc01)
#check the column names
#colnames(df_hap_loc01)
#df_hap_loc01$smplloca
# reshape the data frame for long to wide
df_hap_loc02 <- reshape(data=df_hap_loc01,idvar="smplloca",
                        v.names = "Freq",
                        timevar = "Var2",
                        direction="wide")
#add a column with sampling locations to be able to 
#match between data frames
df_hap_loc02$dec_loc3 <- df_clo$locality[match(df_hap_loc02$smplloca,df_clo2$locality)]
# sum for duplicated values in row
# this is to add up the multiple entries row for the same 
# localities
# https://stackoverflow.com/questions/10180132/consolidate-duplicate-rows
library(plyr)
df_hap_loc03 <- plyr::ddply(df_hap_loc02,"smplloca",numcolwise(sum))
#match between data frames
df_hap_loc03$dec_lat <- df_clo2$dec_lat[match(df_hap_loc03$smplloca,df_clo2$locality4)]
df_hap_loc03$dec_lon <- df_clo2$dec_lon[match(df_hap_loc03$smplloca,df_clo2$locality4)]
# add jitter to points
df_hap_loc03$dec_lat <- jitter(df_hap_loc03$dec_lat, factor = 20.12, amount = NULL)
df_hap_loc03$dec_lon <- jitter(df_hap_loc03$dec_lon, factor = 20.12, amount = NULL)
# limit the data frame to remove any rows that have NAs
df_hap_loc03 <-  df_hap_loc03[complete.cases(df_hap_loc03), ] 
#modify the colnames
colnames(df_hap_loc03) <- c(gsub("Freq\\.","",colnames(df_hap_loc03)))
# https://towardsdatascience.com/using-ggplot-to-plot-pie-charts-on-a-geographical-map-bb54d22d6e13
#count the number of columns, and subtract 2
enc1 <- ncol(df_hap_loc03)-2
#head(df_hap_loc03,3)
#sum up for each row
df_hap_loc03$rws <- rowSums(df_hap_loc03[,2:enc1])
rws <- as.numeric(df_hap_loc03$rws)
#count the columns and subtract 3
enc <- ncol(df_hap_loc03)-3
#install the package that allows for making pit charts in ggplot
if(!require("scatterpie")){
  install.packages("scatterpie", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("scatterpie")
}
library(scatterpie)
#make a viridis colour range
cl03 <- pals::viridis(length(unique(df_hap_loc03[,c(2:enc)])))
# https://towardsdatascience.com/using-ggplot-to-plot-pie-charts-on-a-geographical-map-bb54d22d6e13
# Using map_data()
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
if (!exists("world"))
{  
  #world <- ne_countries(scale = 10, returnclass = "sf")
  world <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
}
#
df_hap_loc03$rws2 <- df_hap_loc03$rws
#replace NAs with zeros
df_hap_loc04 <- df_hap_loc03[!is.na(df_hap_loc03$rws),]
df_hap_loc04$dec_lat <- as.numeric(df_hap_loc04$dec_lat)
df_hap_loc04$dec_lon <- as.numeric(df_hap_loc04$dec_lon)
#https://guangchuangyu.github.io/2016/12/scatterpie-for-plotting-pies-on-ggplot/
world <- ggplot2::map_data('world')
jitlvl <- 0.017
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p07 <- ggplot(data = world) +
  geom_map(map=world, aes(map_id=region), fill="grey",
           color="black") +
  #geom_sf(color = "black", fill = "azure3", lwd=0.4) +
  #geom_sf(color = "black", fill = "azure3") +
  #https://ggplot2.tidyverse.org/reference/position_jitter.html
  # use 'geom_jitter' instead of 'geom_point' 
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.10), 
                              data = df_hap_loc04, 
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  scale_color_manual(values=c(rep("black",
                                  length(unique(df_hap_loc04[,c(2:enc)]))))) +
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  geom_scatterpie_legend(df_hap_loc04$rws*0.10, x=-10, y=47) +
  # set alpha values for color intensity of fill color in point
  #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(6, 14),
                    ylim = c(52.8, 58.0),
                    expand = FALSE)
# change labels on axis
p07 <- p07 + xlab("Longitude") + ylab("Latitude")
# change label for legend
p07 <- p07 + labs(fill='haplotype')

# #https://www.statology.org/ggplot-background-color/
p07 <- p07 + theme(panel.background = element_rect(fill = 'white', color = 'black'),
panel.grid.major = element_line(color = 'grey80')) #, linetype = 'dotted'))#,
#panel.grid.minor = element_line(color = 'green', size = 2))
#p07 <- p07 + theme_bw() #white background and grey gridlines
# see the plot
p07
# Add titles
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p07t <- p07 + labs(title = "a")#,
# ------------- plot Combined figure -------------
library(patchwork)
# set a variable to TRUE to determine whether to save figures
bSaveFigures <- T
#see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
# on how to arrange plots in patchwork
p <-  p07t +
  plot_layout(nrow=1,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") +
  plot_annotation(caption=inpf01) #& theme(legend.position = "bottom")
#p
#make filename to save plot to
figname01 <- paste0("map_haplotype_pie_diagr",inpf01,"_02.png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==F){
  ggsave(p,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}

#make filename to save plot to
figname01 <- paste0("Fig04_v02_map_haplotype_pie_diagr",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}

#______________________________________________________________________________
# Make a second map with pies for haplotypes -end
# Add jitter to the sample localities to see all sampling sites
#______________________________________________________________________________
library(ggplot2)
library(sf)
library(rnaturalearth)

nSco <- length(cl03) 
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# Information on colour blind colours
#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# using only 14 colours
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                             "black","pink")
scbpl <- safe_colorblind_palette
#scales::show_col(safe_colorblind_palette)
# see how to make a number of colurs along color range
# https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
cl2 <- colorRampPalette(c(scbpl))( nSco) 
cl2 


#write.table()
#head(df_r02,4)
#head(df_clo2,4)

#__________________
#

#

#write.table()
#head(df_r02,4)
#head(df_clo2,4)

#__________________
#
mtr_dd_pip <- ape::dist.dna(pip)
tre_pip <- ape::nj(mtr_dd_pip)
#sort the branches in the tree
tre_pipr <- ape::ladderize(tre_pip, right = TRUE)

#https://joey711.github.io/phyloseq/plot_tree-examples.html
plot(tre_pipr, cex=0.4)
# Try and plot a neighbour join tree
#
#https://joey711.github.io/phyloseq/plot_tree-examples.html
#_______________________________________________________________________________
# see this website for inspiration
# https://aschuerch.github.io/posts/2017-04-24-blog-post-1
#https://bioconductor.org/packages/release/bioc/html/ggtree.html
if (!require("BiocManager", quietly = TRUE))
  if(!require(BiocManager)){
    install.packages("BiocManager")
    library(BiocManager)
  }
library(BiocManager)
if(!require(ggtree)){
  BiocManager::install("ggtree")
  library(ggtree)
}
library(ggtree)
library("ggplot2")
library("ggtree")

p <- ggtree(tre_pipr,
            # force the tree to be ladderized right
            right = TRUE,
            branch.length=0.000005,
            #yscale_mapping=0.1,
            # ignore negative branch lengths
            options(ignore.negative.edge=TRUE)) + 
  #xlim(0, 0.12) + # to allow more space for labels
  #ylim(70,0) +
  theme(aspect.ratio=6/3) +
  geom_treescale() # adds the scale

df_tiplb01 <- as.data.frame(cbind(c(tre_pipr$tip.label)))
df_tiplb01$cat <- NA
colnames(df_tiplb01) <- c("seqNm", "cat")
df_tiplb01$cat[grepl("NCBI",df_tiplb01$seqNm)] <- "NCBI"
df_tiplb01$cat[grepl("Germany",df_tiplb01$seqNm)] <- "Germany"
df_tiplb01$cat[grepl("Jylland",df_tiplb01$seqNm)] <- "Jylland"
df_tiplb01$cat[grepl("_Fyn",df_tiplb01$seqNm)] <- "Fyn"
df_tiplb01$cat[grepl("Samsoe",df_tiplb01$seqNm)] <- "Samsoe"
df_tiplb01$cat[grepl("Sjaelland",df_tiplb01$seqNm)] <- "Sjaelland"
tipcategories <- df_tiplb01

#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
colfunc <- colorRampPalette(cbbPalette2)

cl <- cbbPalette
# subset to exclude NCBI sequences
df_clo2 <- df_clo[!grepl("NCBI",df_clo$locality),]
df_clo2$locality3 <- df_clo2$locality
#Modify location names
df_clo2$locality <- gsub("FynBogense","Funen, Bogense", df_clo2$locality )
df_clo2$locality <- gsub("FynKerteminde","Funen, Kerteminde", df_clo2$locality )
df_clo2$locality <- gsub("JyllandMariagerfjord","Jutland, Mariager Fjord", df_clo2$locality )
df_clo2$locality <- gsub("NGermanyKielFjord","Germany, Kiel Fjord", df_clo2$locality )
df_clo2$locality <- gsub("NGermanyMecklenburgerBuchtWismarBucht","Germany, Wismar Bight", df_clo2$locality )
df_clo2$locality <- gsub("NJyllandLimfjord","Jutland, Limfjord", df_clo2$locality )
#df_clo2$locality <- gsub("NJyllandLimfjordLogstoer","Jutland, Limfjord-Løgstør" , df_clo2$locality )
df_clo2$locality <- gsub("NWGermanyNSeaHelgolandRds", "North Sea, Helgoland Roads", df_clo2$locality )
df_clo2$locality <- gsub("NWGermanyWaddenSeaBussumHaupstr","Germany, Büsum", df_clo2$locality )
df_clo2$locality <- gsub("SamsoeBallen","Samsøe, Ballen", df_clo2$locality )
#df_clo2$locality <- gsub("SEDenmarkMecklenburgerBucht","", df_clo2$locality )
df_clo2$locality <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", df_clo2$locality )

#identify unique localities
nloc <- length(unique(df_clo2$locality))
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#
cl <- colfunc(nloc)

dfcol01 <- as.data.frame(cbind(cl,unique(df_clo2$locality)))
colnames(dfcol01) <- c("col","loc")
dd = as.data.frame(tipcategories)
#substitute to get second part of location name
loc2 <- gsub("(.*), (.*)","\\2",dfcol01$loc)
#substitute 
loc2 <- gsub(",","",loc2)
loc2 <- gsub(" ","",loc2)
loc2 <- gsub("Roads","Rds",loc2)
loc2 <- gsub("WismarBight","WismarBucht",loc2)
loc2 <- gsub("MariagerFjord","Mariagerfjord",loc2)
loc2 <- gsub("Büsum","Bussum",loc2)
# make a new column
dd$cat2 <- NA
dd$loc2 <- NA
#count the number of rows
ncol2 <- nrow(dfcol01)
#for element in number rows, match the location, and assign the hex color
# to the new column
for (e in seq(1:ncol2)){
  dd$col2[grepl(loc2[e],dd$seqNm)] <- dfcol01$col[e]
  dd$loc2[grepl(loc2[e],dd$seqNm)] <- dfcol01$loc[e]
}
#get unique NCBI sample sets to assign them colors
unqsmplnm <- unique(gsub("(.*)_(.*)_(.*)_(.*)","\\2",dd$seqNm))
NCBIsmpl <- unqsmplnm[grepl("NCBI",unqsmplnm)]
NCBIsmpl <- c("Mediterranean", 
  "NEAtlantic",
  "CaspianSea", 
  "CentralWAtlantic", 
  "USA:WoodsHole",
  "Germany:Maasholm",
  "NWGermanyNSeaHelgolandRds", 
  "BalticSea", 
  "AtlanticOcean:NWAtlantic", 
  "Holland")
# count the number of elements
noNCBIsmpl <- length(NCBIsmpl)
# make a color range
cbbPalette3 <- c("brown4","brown3","brown2")
cbbPalette3 <- c("azure4","azure3","azure2")
cbbPalette3 <- c("cadetblue4","cadetblue3","cadetblue2")
cbbPalette3 <- cbbPalette1
colfunc2 <- colorRampPalette(cbbPalette3)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#
clNBCI <- colfunc2(noNCBIsmpl)
dfNCBIcol <- as.data.frame(cbind(NCBIsmpl,clNBCI))
#count the number of rows
ncol3 <- nrow(dfNCBIcol)
#for element in number rows, match the location, and assign the hex color
# to the new column
for (e in seq(1:ncol3)){
  dd$col2[grepl(dfNCBIcol$NCBIsmpl[e],dd$seqNm)] <- dfNCBIcol$clNBCI[e]
  dd$loc2[grepl(dfNCBIcol$NCBIsmpl[e],dd$seqNm)] <- dfNCBIcol$NCBIsmpl[e]
}
#subset dataframe to only comprise column 4 and 5
dd2 <- dd[,c(4,5)]
#find unique rows in data frame
dd3 <- dd2 %>% group_by_all %>% count
#retain only column 1 and 2
dd2 <- dd3[,c(1,2)]
#transpose the dataframe and make it a dataframe
dd2 <- as.data.frame(t(dd2))
# change the column names in the transposed data frame
colnames(dd2) <- dd2[1,]
# get the data frame excluding the first row
dd2 <- dd2[-1,]

#Begin plotting the tree
p01 <- p %<+% dd + 
  geom_tiplab(aes(fill = factor(loc2), 
                  color=factor(loc2)
  ),
  #color = "black", # color for label font
  geom = "label",  # labels not text
  #hjust=0.4,
  #align=T,
  offset=0.002,
  size=1.8,
  label.padding = unit(0.09, "lines"), # amount of padding around the labels
  label.size = 0) #+ # size of label border
  #geom_text(aes(label=cat,size=2)) +
  # theme(legend.position = c(0.5,0.6), 
  #       legend.title = element_blank(), # no title
  #       legend.key = element_blank()) # no keys
# get number of categories
nct <- length(unique(factor(dd$cat)))
#
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
colfunc <- colorRampPalette(cbbPalette2)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#
cl <- colfunc(nct)
# apply 
# Use transposed data frame with unique colors


dd4 <- dd2
dd4[1,] <- "black"
dd4[][grepl("Bogense",colnames(dd4))] <- "white"
dd4[][grepl("Kerteminde",colnames(dd4))] <- "white"
dd4[][grepl("Mariager",colnames(dd4))] <- "white"
dd4[][grepl("Kiel Fjord",colnames(dd4))] <- "white"
vl<- as.character(dd4[1,])
p01 <- p01 + scale_colour_manual(values=c(vl))

#
p01 <- p01 + scale_fill_manual(name = "col2", values = alpha(c(dd2),c(0.7) ))  

# p01 <- p01 + labs(color= NULL)
# p01 <- p01 + labs(fill='location')
p01 <- p01 + theme(legend.position = "none")
#p01 <- p01 + coord_fixed(ratio=0.06)
#
p01
#dev.off()
#
#
#make filename to save plot to
figname01 <- paste0("Fig03_v02_NJtree_",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p01,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}

#
#p01 <- p01 + scale_label_manual(name = "col2", values = dd4)  

#
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
colfunc <- colorRampPalette(cbbPalette2)
#identify unique localities
nloc <- length(unique(df_clo2$locality))
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#
cl <- colfunc(nloc)
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p04 <- 
  ggplot(data = denm_map) +
  #ggplot(st_transform(norw_map, 9122)) +
  geom_sf(color = "black", fill = "azure3") +
  theme(aspect.ratio=3/7) +
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.04),
                              data = df_hap_loc03,
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  geom_scatterpie_legend(df_hap_loc04$rws*0.04, x=14, y=57) +
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(4, 16),
                    ylim = c(53.4, 58.4),
                    expand = FALSE)
# change label for legend - Notice that you need to change for all 3 variables
# you called 'aes' in 'geom_jitter'
p04 <- p04 + labs(fill='Haplotype')
p04 <- p04 + labs(color='Haplotype')
p04 <- p04 + labs(shape='Haplotype')
p04 <- p04 + xlab("Longitude") + ylab("Latitude")
# #https://www.statology.org/ggplot-background-color/
p04 <- p04 + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                   panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,
#panel.grid.minor = element_line(color = 'green', size = 2))
# see the plot
p04
#
#
#
#
#_______________________________________________________________________________
#_______________________________________________________________________________
# start - Make haplotype circles on map second attempt 
#_______________________________________________________________________________
#make a viridis colour range
cl03 <- pals::viridis(length(unique(df_hap_loc05[,c(2:enc)])))
cl03 <- pals::inferno(length(unique(df_hap_loc05[,c(2:enc)])))
ncolHptloc03 <- length(unique(colnames(df_hap_loc05)[2:enc]))
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
colfunc <- colorRampPalette(cbbPalette2)

#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#
cl03 <- colfunc(ncolHptloc03)
# https://towardsdatascience.com/using-ggplot-to-plot-pie-charts-on-a-geographical-map-bb54d22d6e13
# Using map_data()
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
if (!exists("world3"))
{  
  #world <- ne_countries(scale = 10, returnclass = "sf")
  world3 <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
}
# grep for Mediterranean sea to check position
df_clo03[grepl("Mediterr",df_clo03$locality5),]
#replace NAs with zeros
df_hap_loc06 <- df_hap_loc05[!is.na(df_hap_loc05$rws),]
df_hap_loc06$dec_lat <- as.numeric(df_hap_loc06$dec_lat)
df_hap_loc06$dec_lon <- as.numeric(df_hap_loc06$dec_lon)
df_hap_loc06[grepl("Hl",df_hap_loc06$dec_loc2),]

#https://guangchuangyu.github.io/2016/12/scatterpie-for-plotting-pies-on-ggplot/
world3 <- ggplot2::map_data('world3')
jitlvl <- 0.017
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p08 <- 
  ggplot(data = world3) +
  #ggplot(st_transform(norw_map, 9122)) +
  geom_sf(color = "black", fill = "azure3") +
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.30), 
                              data = df_hap_loc06, 
                              cols = colnames(df_hap_loc06[,c(2:enc)])) +
  
  scale_color_manual(values=c(rep("black",
                                  length(unique(df_hap_loc06[,c(2:enc)]))))) +
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
  #scale_colour_viridis_d(breaks=11) +
  #scale_color_brewer(palette="Dark2") +
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  theme(aspect.ratio=7/8) +
  geom_scatterpie_legend(df_hap_loc06$rws*0.30, x=-10, y=47, 
                         labeller = function(ra) {ra * (1/0.10)}) +
  # set alpha values for color intensity of fill color in point
  #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(-12, 32),
                    ylim = c(34.8, 58.0),
                    expand = FALSE)
# change label for legend - Notice that you need to change for all 3 variables
# you called 'aes' in 'geom_jitter'
p08 <- p08 + labs(fill='Haplotype')
p08 <- p08 + labs(color='Haplotype')
p08 <- p08 + labs(shape='Haplotype')
p08 <- p08 + xlab("Longitude") + ylab("Latitude")
# #https://www.statology.org/ggplot-background-color/
p08 <- p08 + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                   panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,
#panel.grid.minor = element_line(color = 'green', size = 2))
# see the plot
p08


# also see : https://github.com/tidyverse/ggplot2/issues/2037
p09 <- ggplot(data = denm_map) +
  #ggplot(st_transform(norw_map, 9122)) +
  geom_sf(color = "black", fill = "azure3") +
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.05), 
                              data = df_hap_loc04, 
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  
  scale_color_manual(values=c(rep("black",
                                  length(unique(df_hap_loc04[,c(2:enc)]))))) +
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
  #scale_colour_viridis_d(breaks=11) +
  #scale_color_brewer(palette="Dark2") +
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  theme(aspect.ratio=5/8) +
  #https://stackoverflow.com/questions/71277499/adding-a-plot-legend-with-geom-scatterpie-legend-in-r
  geom_scatterpie_legend(df_hap_loc04$rws*0.05, x=13.4, y=57.0,labeller = function(ra) ra * 1/0.05) +
  #unique(df_hap_loc04$rws)[order(unique(df_hap_loc04$rws))]
  # set alpha values for color intensity of fill color in point
  #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
  #define limits of the plot
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(6, 15.4),
                    ylim = c(52.8, 58.0),
                    expand = FALSE)
#https://www.statology.org/ggplot-background-color/
#p09 <- p09 + theme_minimal() #no background annotations
# #https://www.statology.org/ggplot-background-color/
p09 <- p09 + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                   panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,

# see the plot
#p09
#dev.off()
#change labels on axis
p09 <- p09 + xlab("Longitude") + ylab("Latitude")
#change labels on legend
p09 <- p09 + labs(fill='Haplotype')
p09
# Add titles
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p08t <- p08 + labs(title = "b") +
  theme(plot.title = element_text(face="bold"))#,
# Add titles
# p09t <- p09 + labs(title = "eDNA samples attempted",
#                    subtitle = "at least approv controls and 1 or 2 pos repl")#,
p09t <- p09 + labs(title = "a") +
  theme(plot.title = element_text(face="bold"))#,


# ------------- plot Combined figure -------------
library(patchwork)
# set a variable to TRUE to determine whether to save figures
bSaveFigures <- T
#see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
# on how to arrange plots in patchwork
p <-  p09t +
  p08t +
  
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") #+
  #plot_annotation(caption=inpf01) #& theme(legend.position = "bottom")
#p

#make filename to save plot to
figname01 <- paste0("Fig04_v05_map_haplotype_pie_diagr",inpf01,".png")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}


#make filename to save plot to
figname01 <- paste0("Fig04_v06_map_haplotype_pie_diagr",inpf01,".pdf")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}
#_______________________________________________________________________________
# end - Make haplotype circles on map second attempt 
#_______________________________________________________________________________
# paste path and file together
loc04fl <- paste(wd00_wd05,"/loc04.csv",sep="")
# write csv
write_csv(df_hap_loc04,file=loc04fl)
df_hap_loc08 <- df_hap_loc04


#df_hap_loc06 <- df_hap_loc07[complete.cases(df_hap_loc07),]
# paste path and file together
loc06fl <- paste(wd00_wd05,"/loc06.csv",sep="")
# write csv
write_csv(df_hap_loc09,file=loc06fl)





#_______________________________________________________________________________
# match between abbreviated location names to get long location names
df_hap_loc06$dec_loc3 <- df_clo03$locality5[match(df_hap_loc06$dec_loc2,df_clo03$locality2 )]
# add a new colunm to add overall geographic regions to
df_hap_loc06$ov.aL <- df_hap_loc06$dec_loc3
df_hap_loc09 <- df_hap_loc06 
#df_hap_loc06 <- df_hap_loc07[complete.cases(df_hap_loc07),]
# paste path and file together
loc06fl <- paste(wd00_wd05,"/loc06.csv",sep="")
# write csv
write_csv(df_hap_loc09,file=loc06fl)

# start adding overall geographical regions based on matches
df_hap_loc06$ov.aL[grepl("Baltic",df_hap_loc06$dec_loc3)] <- "NEurope"
df_hap_loc06$ov.aL[grepl("Germany",df_hap_loc06$dec_loc3)] <- "NEurope"
df_hap_loc06$ov.aL[grepl("Jutland",df_hap_loc06$dec_loc3)] <- "NEurope"
df_hap_loc06$ov.aL[grepl("Funen",df_hap_loc06$dec_loc3)] <- "NEurope"
df_hap_loc06$ov.aL[grepl("North Sea",df_hap_loc06$dec_loc3)] <- "NEurope"
df_hap_loc06$ov.aL[grepl("NEAtlantic",df_hap_loc06$dec_loc3)] <- "NEurope"
df_hap_loc06$ov.aL[grepl("Samsoe",df_hap_loc06$dec_loc3)] <- "NEurope"
df_hap_loc06$ov.aL[grepl("Sealand",df_hap_loc06$dec_loc3)] <- "NEurope"
df_hap_loc06$ov.aL[grepl("Holland",df_hap_loc06$dec_loc3)] <- "NEurope"

df_hap_loc06$ov.aL[grepl("NWAtlantic",df_hap_loc06$dec_loc3)] <- "NWAtlantic"
df_hap_loc06$ov.aL[grepl("USA",df_hap_loc06$dec_loc3)] <- "NWAtlantic"

# sum up Hpt counts from different sampling locations within same geographic region
df_hap_loc07 <- aggregate(df_hap_loc06[,sapply(df_hap_loc06,
                                               is.numeric)],df_hap_loc06["ov.aL"],sum)
#replace the lat lon positions as they have been summed up
df_hap_loc07$dec_lat <- df_hap_loc06$dec_lat[match(df_hap_loc07$ov.aL,df_hap_loc06$ov.aL)]
df_hap_loc07$dec_lon <- df_hap_loc06$dec_lon[match(df_hap_loc07$ov.aL,df_hap_loc06$ov.aL)]

#count the columns 
enc <- ncol(df_hap_loc07)-4
#make color range for haplotypes 
cl07 <- pals::inferno(length(unique(df_hap_loc07[,c(2:enc)])))
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p10 <- 
  ggplot(data = world3) +
  #ggplot(st_transform(norw_map, 9122)) +
  geom_sf(color = "black", fill = "azure3") +
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = rws*0.20), 
                              data = df_hap_loc07, 
                              cols = colnames(df_hap_loc07[,c(2:enc)])) +
  
  scale_color_manual(values=c(rep("black",
                                  length(unique(df_hap_loc07[,c(2:enc)]))))) +
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
  #scale_colour_viridis_d(breaks=11) +
  #scale_color_brewer(palette="Dark2") +
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  #theme(aspect.ratio=7/8) +
  geom_scatterpie_legend(df_hap_loc06$rws*0.20, x=-40, y=40, 
                         labeller = function(ra) {ra * (1/0.10)}) +
  # set alpha values for color intensity of fill color in point
  #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(-80, 60),
                    ylim = c(80, 2.0),
                    expand = FALSE)
# change label for legend - Notice that you need to change for all 3 variables
# you called 'aes' in 'geom_jitter'
p10 <- p10 + labs(fill='Haplotype')
p10 <- p10 + labs(color='Haplotype')
p10 <- p10 + labs(shape='Haplotype')
p10 <- p10 + xlab("Longitude") + ylab("Latitude")
# #https://www.statology.org/ggplot-background-color/
p10 <- p10 + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                   panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,
#panel.grid.minor = element_line(color = 'green', size = 2))
# see the plot
#p10

p09 <- p09 + xlab("Longitude") + ylab("Latitude")
#change labels on legend
p09 <- p09 + labs(fill='Haplotype')
#p09
p10 <- p10 + xlab("Longitude") + ylab("Latitude")
#change labels on legend
p10 <- p10 + labs(fill='Haplotype')
#p10
# Add titles
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p10t <- p10 + labs(title = "b") +
  theme(plot.title = element_text(face="bold"))#,
# Add titles
# p09t <- p09 + labs(title = "eDNA samples attempted",
#                    subtitle = "at least approv controls and 1 or 2 pos repl")#,
p09t <- p09 + labs(title = "a") +
  theme(plot.title = element_text(face="bold"))#,


# ------------- plot Combined figure -------------
library(patchwork)
# set a variable to TRUE to determine whether to save figures
bSaveFigures <- T
#see this website: https://www.rdocumentation.org/packages/patchwork/versions/1.0.0
# on how to arrange plots in patchwork
p <-  p09t +
  p10t +
  
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") #+
#plot_annotation(caption=inpf01) #& theme(legend.position = "bottom")
#p

#make filename to save plot to
figname01 <- paste0("Fig04_v07_map_haplotype_pie_diagr",inpf01,".png")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}



#_______________________________________________________________________________
# Make principal component analysis
#_______________________________________________________________________________

# See these webpages for inspiration:
#https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html
#https://knausb.github.io/vcfR_documentation/visualization_1.html
#https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
#citation("dartR")
#Principal components analysis
#A principal components analysis (PCA) converts the observed SNP data
#into a set of values of linearly uncorrelated variables called principal 
#components that summarize the variation between samples. 
#We can perform a PCA on our genlight object by using the glPCA function.
# see barplot first

# convert genind object to genlight object with the dartR package
gl_pip <- dartR::gi2gl(gei_pip)
# Make PCA of genlight object with adegenet
h.pca <- adegenet::glPca(gl_pip, nf = 3)
# see it as a bar plot
barplot(100*h.pca$eig/sum(h.pca$eig),
        col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
# second make PCA plot
h.pca.scores <- as.data.frame(h.pca$scores)
#get pop groups from row names
popyear <- gsub("(^.*)_(.*)_(.*)_(.*)$","\\3",row.names(h.pca.scores))
poploc <- gsub("(^.*)_(.*)_(.*)_(.*)$","\\2",row.names(h.pca.scores))
#h.pca.scores$pop <- popyear
h.pca.scores$pop <- poploc
# match to get longer location names on data frame with colors
df_cll$coll_loc2 <- df_clo2$locality3[match(df_cll$coll_loc,df_clo2$locality)]
# add one long location name that is missing
df_cll$coll_loc2[df_cll$coll_loc=="Funen, Bogense"] <- "FynBogense"
# use the first 4 values in one column and add to the first 4 rows of a another
# colunm
df_cll$coll_loc2[1:4] <- df_cll$coll_loc[1:4]
# sort the unique pop locations
poplocs2 <- unique(poploc)[order(unique(poploc))]
df_cll$coll_loc3 <- df_cll$coll_loc
#poplocs2[!poplocs2 %in% df_cll$coll_loc3]
poplocs3 <- df_clo03$locality4[match(poplocs2,df_clo03$locality6)]

# now match between data frames to get colors for locations
#cl7 <- df_cll$colfcol_loc[match(poplocs2,df_cll$coll_loc2)]
cl7 <- df_cll$colfcol_loc[match(poplocs3,df_cll$coll_loc)]
col8 <- as.data.frame(cbind(poplocs3,poplocs2,cl7))

col8$cl7[grepl("Kiel",col8$poplocs2)] <- df_cll$colfcol_loc[grepl("Kiel",df_cll$coll_loc)]
col8$cl7[grepl("Bussum",col8$poplocs2)] <- df_cll$colfcol_loc[grepl("Büsum",df_cll$coll_loc)]
cl7 <- col8$cl7
# match back detailed location name
poploc3 <- df_cll$coll_loc[match(poplocs3,df_cll$coll_loc3)]

#replace location names in pca data frame
h.pca.scores$pop <- poploc3
library(ggplot2)
set.seed(9)
#make a plot with pca 
p <- ggplot(h.pca.scores, aes(x=PC1, y=PC2, colour=pop, fill=pop)) 
#p <- p + geom_point(size=2)
p <- p + geom_point(size=2,shape=21, color="black")
p <- p + stat_ellipse(level = 0.95, size = 1)
#p <- p + scale_color_manual(values = cl7) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p01 <- p
#p01
# color range from haplotype network pie diagrams for
# location sampled - used in Fig.2.
#colfh
p01 <- p01 + scale_color_manual(values = cl7) 
p01 <- p01 + scale_fill_manual(values = cl7) 
# change the heading for the legend, this must be done for all 
#settings for the points
p01 <- p01 + labs(color='Location')
p01 <- p01 + labs(fill='Location')
p01 <- p01 + labs(shape='Location')

# replace the pop name column in the pca data frame with the sampling year
h.pca.scores$pop <- popyear
#h.pca.scores$pop <- poploc

library(ggplot2)
set.seed(9)
p <- ggplot(h.pca.scores, aes(x=PC1, y=PC2, colour=pop,fill=pop)) 
#p <- p + geom_point(size=2)
p <- p + geom_point(size=2,shape=21, color="black")
p <- p + stat_ellipse(level = 0.95, size = 1)
#p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p02 <- p

# color range from haplotype network pie diagrams for
# locationyear sampled - used in Fig.2.
#clrf6
# set the color of the fill for lines and the fill for inside fill for the points
p02 <- p02 + scale_color_manual(values = clrf6) 
p02 <- p02 + scale_fill_manual(values = clrf6) 

p02 <- p02 + labs(color='Sampling year')
p02 <- p02 + labs(fill='Sampling year')
p02 <- p02 + labs(shape='Sampling year')
#p02
library(patchwork)
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p01t <- p01 + labs(title = "a", face="bold")#,
p02t <- p02 + labs(title = "b", face="bold")#,


pA <-  p01t +
  p02t +
  
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") #+
  #plot_annotation(caption=pthinf01) #& theme(legend.position = "bottom")
#p
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig05_pca_M_leyidy.png")
figname01 <- paste0("Fig05_pca_M_leyidy.jpg")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}


#_______________________________________________________________________________
#_______________________________________________________________________________
figname01 <- "Fig06_v01_network.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 5000, height = 4000, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 2), xpd=FALSE)
par(oma=c(1, 1, 0, 0))
#reset this parameter
par(mfrow = c(1, 1)) 


# Make 
# Minimum spanning networks
library(igraph)
# Add predetermined populations. 
# convert genind object to genlight object with the dartR package
gl_pip <- dartR::gi2gl(gei_pip)
gei_pip2 <- gei_pip
pop(gei_pip2) <- gsub("^(.*)_(.*)_(.*)_(.*)$","\\2",rownames(gei_pip$tab))
gl_pip <- dartR::gi2gl(gei_pip2)
poploc2 <- gsub("^(.*)_(.*)_(.*)_(.*)$","\\2",rownames(gei_pip$tab))
pop(gl_pip) <- poploc2
gl_pip$ind.names <- gsub("^(.*)_(.*)_(.*)_(.*)$","\\1",rownames(gei_pip$tab))
indnms.pip <- gl_pip$ind.names
pip.dist <- bitwise.dist(gl_pip)
pip.msn <- poppr.msn(gl_pip, pip.dist, showplot = T, 
                     include.ties = T, vertex.label = "inds")

node.size <- rep(2, times = nInd(gl_pip))
names(node.size) <- indNames(gl_pip)

# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
# make tables with PhiST values

#_______________________________________________________________________________
#
df.p.lo <- df_pw_wc_pip
df.p.ym <- df_pw_wc_pip_ym

df.p.lo <- round(df.p.lo,2)
df.p.ym <- round(df.p.ym,2)
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
# replace inm column names and row names with abbreviated names
rownames(df.p.lo) <- df_clo03$locality2[match(rownames(df.p.lo),df_clo03$locality6)]
colnames(df.p.lo) <- df_clo03$locality2[match(colnames(df.p.lo),df_clo03$locality6)]
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
  c("FBo",
"FKe",
"JMa",
"GKi",
"GMe",
"JLi",
"JLo",
"GHe",
"GWa",
"SBa",
"DMe",
"Ssk")
#combine to a data frame
df_abblocnms <- as.data.frame(cbind(llocnms,abbloccnms))
#substitue in data frame
df_abblocnms$locnm2 <- gsub(" ","",df_abblocnms$llocnms)
# match between data frames
lloc6 <- df_clo03$locality6[match(rownames(df.p.lo),df_clo03$locality2)]
ll7 <- df_abblocnms$abbloccnms[match(lloc6,df_abblocnms$locnm2)]
# combine to a data frame
df_ll08 <- as.data.frame(cbind(lloc6,ll7))
# match between the value that are NA
df_ll08$ll7[is.na(df_ll08$ll7)] <- df_clo03$locality2[match(df_ll08$lloc6[is.na(df_ll08$ll7)],df_clo03$locality4)]
df_clo03$locality7 <- df_clo03$locality2
df_clo03$locality8 <-  df_ll08$ll7[match(df_clo03$locality,df_ll08$lloc6)]
df_clo03$locality8[is.na(df_clo03$locality8)] <- df_clo03$locality7[is.na(df_clo03$locality8)]
df_clo03$locality8[df_clo03$locality8=="FB"] <- "FBo"
#replace abbreviated colunm names and row names
colnames(df.p.lo) <- df_clo03$locality8[match(colnames(df.p.lo),df_clo03$locality2)]
rownames(df.p.lo) <- df_clo03$locality8[match(rownames(df.p.lo),df_clo03$locality2)]
# paste together location and abbreviation to use for the tables
v_locabb <- paste0(df_clo03$locality4," (",df_clo03$locality8,");")
# paste all elements in vector together in one string
locat_abb <- paste(v_locabb,collapse=" ")
#_______________________________________________________________________________
#
#_______________________________________________________________________________
# Make a html table with PhiST for locations
#_______________________________________________________________________________

library(htmlTable)
d1 <- df.p.lo
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
#https://www.rdocumentation.org/packages/htmlTable/versions/2.4.0
#https://www.r-bloggers.com/2015/04/introducing-the-htmltable-package/
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm>0.3,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
#
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
capt_tbl01 <- paste0("Table 2. Comparison of population genetic fixation index (PhiST) obtained from sequences of nDNA 18S for various sampling locations. The PhiST values can range from 0 to 1, where a PhiST of 1 indicates very different genetic variations among sampled individuals. Low PhiST values close to 0 suggest that individuals sampled are genetically very similar. Abbreviations for locations are: ",
                     locat_abb)
                     
# show the table
t.HTML01 <- mdlo %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl01)
# class(t.HTML01)
#see the table
t.HTML01
# #and to export in a file
# tableHTML::write_tableHTML(t.HTML01, file = paste(wd00_wd05,"/Table02a_PhiST_locality.html",sep=""))
# #vignette("tables", package = "htmlTable")
# pth_flnm_html01 <- paste(wd00_wd05,"/Table02a_PhiST_locality.html",sep="")
#_______________________________________________________________________________
# Make a html table with PhiST for sampling year
#_______________________________________________________________________________
#_______________________________________________________________________________
library(htmlTable)
d1 <- df.p.ym
# change column and row names to ensure they are not evaluated when
# min and max values are used for coloring cells
colnames(d1) <-  paste0("y_",colnames(d1))
rownames(d1) <- paste0("y_",rownames(d1))
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm>0.3,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")
# paste together background color og font color
colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
#
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
capt_tbl02 <-        "Table 3. Comparison of PhiST for Mnemiopsis leidyi obtained from sequences of nDNA 18S for the sampling years. Sequences of nDNA 18S obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank."
# show the table
t.HTML02 <- mdlo %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02)
t.HTML02

