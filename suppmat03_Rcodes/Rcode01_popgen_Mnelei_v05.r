#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# This code is able  to run in:

#remove everything in the working environment, without a warning!!
#rm(list=ls())

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

# #read in packages and install libraries
# #install.packages("ape")
# if(!require(ape)){
#   install.packages("ape")
# }
#read in the 'ape' library
library(ape)
#read in the 'pegas' library
library(pegas)
# #install.packages("ggforce")
# if(!require(ggforce)){
#   install.packages("ggforce")
# }
library(sf)
#install package if needed
# if(!require(rnaturalearth)){
#   install.packages("rnaturalearth")
#   #install.packages("rnaturalearth", repos = "http://packages.ropensci.org", type = "source")
# }
##Get 'SNPRelate' in order to be able to install 'DartR'
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# if(!require(bioseq)){
# remotes::install_github("fkeck/bioseq")
#  }  
# library(bioseq)
#install.packages("scatterpie")
# if(!require(scatterpie)){
#   install.packages("scatterpie")
# }

library(ggforce)
library(scatterpie)

# if(!require(pegas)){
#   # make sure you have Rtools installed
#   if (!require('devtools')) install.packages('devtools')
#   # install from GitHub
#   devtools::install_github("emmanuelparadis/pegas/pegas")
# }
library(pegas)

library(devtools)
# if(!require(bioseq)){
#   # make sure you have Rtools installed
#   if (!require('devtools')) install.packages('devtools')
#   # install from GitHub
#   remotes::install_github("fkeck/bioseq")
# }
library(bioseq)
# to install 'Biostrings' the 'IRanges' package is required 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("IRanges")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
library(Biostrings)
# Get the treio package as described here : https://www.bioconductor.org/packages/release/bioc/html/treeio.html
# to be able to use the phylogenetic tree labelling
# # example with the 'BEAST' file here : https://yulab-smu.top/treedata-book/chapter5.html
# if (!require("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
#   BiocManager::install("treeio", force=T)}

library(treeio)
library(rnaturalearthhires)
#install.packages("rnaturalearthhires")
# if(!require(ggtree)){
#   BiocManager::install("ggtree")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")
# }
library(ggtree)
library(ggplot2)
# load additional packages
library(tidyr)
library(dplyr)

# load library
# if(!require(sf)){
#   install.packages("sf")
# }
library(sf)
#install package if needed
# if(!require(rnaturalearth)){
#   install.packages("rnaturalearth")
#   #install.packages("rnaturalearth", repos = "http://packages.ropensci.org", type = "source")
# }
# Get 'SNPRelate' in order to be able to install 'DartR'
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# library(rnaturalearth)
# install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
# install.packages("rnaturalearthdata", repos = "http://packages.ropensci.org", type = "source")
# # 
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
# You will need 'terra' to install 'dartR', and 'terra' requires
# 'gdal', and 'gdal' can be installed in a terminal
# $ sudo apt install libgdal-dev
#install.packages("terra")
library(ape)
library(adegenet)
library(dartR)
#install.packages('dartR')
library("StAMPP")

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
#wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
wd00 <- getwd()

wd00_wd01 <- paste0(wd00,wd01)
wd00_wd05 <- paste0(wd00,wd05)
# delete a directory -- must add recursive = TRUE
unlink(wd00_wd05, recursive = TRUE)
#create anew directory
dir.create(wd00_wd05)
#set the working dir
#setwd(wd00_wd01)
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

# paste together to get DMS coordinates, notice the separator needs to
# be a space for the  'angle2dec' - function to be able to work
df_lN02$lat01<- paste(df_lN02$ll1,
                      df_lN02$ll2,
                      0, sep=" ")
df_lN02$lon01<- paste( df_lN02$ll4,
                       df_lN02$ll5,
                       0, sep=" ")
#_______________________________________________________________________________
# function : 'angle2dec' - start
#_______________________________________________________________________________
# get function from this website
#https://stackoverflow.com/questions/30879429/how-can-i-convert-degree-minute-sec-to-decimal-in-r
angle2dec <- function(angle) {
  angle <- as.character(angle)
  x <- do.call(rbind, strsplit(angle, split=' '))
  x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    y[1] + y[2]/60 + y[3]/3600
  })
  return(x)
}
#_______________________________________________________________________________
# function : 'angle2dec' - end
#_______________________________________________________________________________
# Use the 'angle2dec' function on the DMS pasted coordinates
df_lN02$lat02 <- angle2dec(df_lN02$lat01)
df_lN02$lon02 <- angle2dec(df_lN02$lon01)
# modify to a negative lat or long if the hemisphere is not N or E
df_lN02$lat02 <- ifelse(df_lN02$ll3=="N",df_lN02$lat02, df_lN02$lat02*-1)
df_lN02$lon02 <- ifelse(df_lN02$ll6=="E",df_lN02$lon02, df_lN02$lon02*-1)
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
df_clo$DMSlat<- paste( df_clo$lat_deg,
                       df_clo$lat_min,
                       df_clo$lat_sec,
                       sep=" ")
df_clo$DMSlon<- paste( df_clo$lon_deg,
                       df_clo$lon_min,
                       df_clo$lon_sec,
                       sep=" ")
# Use the 'dms2dec' function on the DMS pasted coordinates
df_clo$dec_lat <- angle2dec(df_clo$DMSlat)
df_clo$dec_lon <- angle2dec(df_clo$DMSlon)

# modify to a negative lat or long if the hemisphere is not N or E
df_clo$dec_lat <- ifelse(df_clo$lat_sph=="N",df_clo$dec_lat, df_clo$dec_lat*-1)
df_clo$dec_lon <- ifelse(df_clo$lon_sph=="E",df_clo$dec_lon, df_clo$dec_lon*-1)

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
    "JLim",
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
            "JLim", "FBog", "GWis", 
            "FKer", "JMar", 
            "SSko", "GBus", 
            "BalS", "CasS", "CWAt", "FBog", 
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
# section 02 -  start - plot collection points
#_______________________________________________________________________________
# # Get a map, use a high number for 'scale' for a coarse resolution
# use a low number for scale for a high resolution
# if the map 'world' does not exist, then download it
if (!exists("denm_map"))
{
  denm_map <- rnaturalearth::ne_countries(country=c(
    "denmark",
    "sweden",
    "germany",
    "norway", 
    "poland"),scale = 10, returnclass = "sf")
  #wor_map <- ne_countries(country="world",scale = 10, returnclass = "sf")
}

# make a vector that holds the sampled locations
# for Denmark and Germany
DK_Germ_collc <- c(#"BalS",
  #"CasS",
  #"CWAt",
  "FBog",
  "FKer",
  "GBus",
  "GHel",
  "GKie",
  "GWis",
  "JLim",
  "JMar",
  #"Medi",
  #"NEAt",
  #"Neth",
  #"NWAt",
  "SBal",
  "SSko" )#,
#"USAG",
#"USAP",
#"USAW")

# subset to only comprise the samples collected from Denmark and Germany
df_clo04.DG <- df_clo04[(df_clo04$loc4Lett  %in% DK_Germ_collc),]
# Make some long location names , and place in a vector
# make it a matrix, and turn the matrix into a data frame
df_LfCab <- as.data.frame(matrix(c(       
  "FBog","Funen, Bogense",
  "FKer","Funen, Kerteminde",
  "GBus","Germany, Büsum",
  "GHel","North Sea, Helgoland Roads",
  "GKie","Germany, Kiel Fjord",
  "GWis","Germany, Wismar Bight",
  "JLim","Jutland, Limfjord",
  "JMar","Jutland, Mariager Fjord",
  "SBal","Samsøe, Ballen",
  "SSko","Sealand, Skovshoved") ,
  ncol=2, byrow=T))
# change the column headers
colnames(df_LfCab) <- c("loc4Lett","LngNm")
# re order the data frame to have the colors that will follow
# match the alphabetical order of the locations
df_LfCab <-  df_LfCab[order(df_LfCab$LngNm),]
# make a range of colors to use with the'colorRampPalette' 
# to get colors per sampled location
# define different combinations of colors to reflect
# the labels that are to be assigned to the propotions
# in th pie diagrams -  these colors will used in the 
# iterations over gene names here below
cbbPalette1 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
# As an alternative use this palette of grey hues
cbbPalette4 <- c("gray40","gray65",#"orange",
                 "gray80", "cyan1","cyan3","cyan4","deepskyblue4")

# As an alternative use this palette for gene variants
cbbPalette5 <- rev(c("white",
                     "yellow",#"brown3"))
                     "gold","orange",
                     "tomato","firebrick3","brown","cornflowerblue",
                     "deepskyblue3","grey34","orchid4"))

cbbPalette6 <- c("white","yellow","orange","tomato","red","brown","black")
# cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
# make a colr function  based on the color palette
colfunc <- colorRampPalette(cbbPalette2)
# get the count of locationc
nloc <- length(unique(df_LfCab$loc4Lett)) 
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# use the count of locations to make a color range to reflect the 
# sampled locations
colf.locNm <- colfunc(nloc)
# append this range of colors as a column to the data frame
df_LfCab$colf.locNm <- colf.locNm 
# match between data frames to get the long name for the collection
# locality
df_clo04.DG$LngNm <- df_LfCab$LngNm[
  match(df_clo04.DG$loc4Lett,df_LfCab$loc4Lett)]
# get the long names, and place in a vector
LngNm <- df_LfCab$LngNm
# make the color categories for the names
clfLngNm <- rbind(LngNm,colf.locNm)
colnames(clfLngNm) <- clfLngNm[1,]
clfLngNm <-clfLngNm[-1,]
mtx_clfLn<- matrix(clfLngNm)
# use ggplot2 to plot the sampled locations on a map
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p04 <- 
  ggplot(data = denm_map) +
  geom_sf(color = "black", fill = "azure3") +
  geom_point(data = df_clo04.DG, 
             aes(x = dec_lon2, 
                 y = dec_lat2, 
                 fill=LngNm),
             shape=21,
             size = 6) +
  #set the color of the points
  #use alpha to scale the intensity of the color
  scale_fill_manual(values=alpha( c(clfLngNm),   c(0.7) )) +
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(4, 16.4),
                    ylim = c(53.4, 58.4),
                    expand = FALSE)
# change label for legend - Notice that you need to change 
# for all 3 variables
p04 <- p04 + labs(fill='Location')
p04 <- p04 + labs(color='Location')
p04 <- p04 + labs(shape='Location')
p04 <- p04 + xlab("Longitude") + ylab("Latitude")
# #https://www.statology.org/ggplot-background-color/
p04 <- p04 + theme(panel.background = element_rect(fill = 'white', 
                                                   color = 'white'),
                   panel.grid.major = element_line(color = 'white')) #, linetype = 'dotted'))#,
# add border around plot
p04 <- p04 + theme(panel.border = element_rect(color = "black",
                                               fill = NA,
                                               size = 1.0))
# change background of legend
p04 <- p04 + theme(legend.key = element_rect(fill = "white"))
# see the plot
p04
#make filename to save plot to
fileNm.Fig <- paste0("Fig01_map_samples_",inpf01,".png")
# paste the path on to the file name
fileNm.Fig <- paste(wd00_wd05,"/",fileNm.Fig,sep="")
# set parameter for evaluating whether the figure should be saved
bSaveFigures <- T
# save the file if the 'bSaveFigures' is set to TRUE 
if(bSaveFigures==T){
  ggsave(p04,
         file=fileNm.Fig,
         width=210,height=297*0.5,
         units="mm",
         dpi=300)
}

# make a vector for all the other collection locations
oth_collc <- c(
  "BalS",
  "CasS",
  "CWAt",
  #"FBog",
  #"FKer",
  #"GBus",
  #"GHel",
  #"GKie",
  #"GWis",
  #"JLim",
  #"JMar",
  "Medi",
  "NEAt",
  "Neth",
  "NWAt",
  #"SBal",
  #"SSko" )#,
  "USAG",
  "USAP",
  "USAW")
# re order the vector to ensure all locations are sorted alphabetically
oth_collc <- oth_collc[order(oth_collc)]
# subset to only comprise the samples collected from Denmark and Germany
df_clo04.OC <- df_clo04[(df_clo04$loc4Lett  %in% oth_collc),]
# make a colr function  based on the color palette
colfunc <- colorRampPalette(cbbPalette4)
# get the count of locationc
nloc <- length(unique(df_clo04.OC$loc4Lett)) 
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# use the count of locations to make a color range to reflect the 
# sampled locations
colf.locNm <- colfunc(nloc)
# append this range of colors as a column to the data frame
df_clo04.OC$colf.locNm <- colf.locNm 
# make a column where the colors for locations can be stored
df_clo04$colf.locNm <- NA
df_clo04$colf.locNm[(df_clo04$loc4Lett %in% df_clo04.OC$loc4Lett)] <- df_clo04.OC$colf.locNm
df_clo04$colf.locNm[(df_clo04$loc4Lett %in% df_LfCab$loc4Lett)] <- df_LfCab$colf.locNm
colf.locNm <- df_clo04$colf.locNm
loc4Lett <- df_clo04$loc4Lett
# make the color categories for the location names
clfabNm <- rbind(loc4Lett ,colf.locNm)
colnames(clfabNm) <- clfabNm[1,]
clfabNm <-clfabNm[-1,]
mtx_clfabNm<- matrix(clfabNm)

# assuming that there will not be more than 5 gene variants
nGv <- 5
# make a color function  based on the color palette
colfunc <- colorRampPalette(cbbPalette5)
colf.Gva <- colfunc(nGv)
# make the color categories for the names
clfGvar <- rbind(paste0("var",LETTERS[seq(1,nGv,1,)]) ,colf.Gva)
clfGvar <- rbind(paste0("var",seq(1,nGv,1)) ,colf.Gva)
colnames(clfGvar) <- clfGvar[1,]
clfGvar <-clfGvar[-1,]


#_______________________________________________________________________________
# section 02 -  end - plot collection points
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 03 -  start - read in fasta file with alignment, identify coding
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
wds03Rc <- (paste0(wd00,"/suppmat03_Rcodes"))
source(paste0(wds03Rc,"/","Rfunction_DNA_sequence_subset.R"))
source(paste0(wds03Rc,"/","Rfunction_ReadFasta.R"))

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
# section 03 -  end - read in fasta file with alignment, identify coding
# regions, trim alignment, and save only coding regions alignments
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 04 -  start - make all disambiguity variants, and rename
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
    vrNms <- paste0(sqNm,"_var",sqn.disam)
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

# Get all names from both alignments, 
NmsITS2 <- names(alvrs.M.ITS2)
NmsITS1 <- names(alvrs.M.ITS1)
NmsITS.Gn <- c(NmsITS1,NmsITS2)
# and retain unique names
NmsITS.Gn <- NmsITS.Gn[!duplicated(NmsITS.Gn)]
# split by the delimter, and get each part of the sequence name
# arranged in columns
NmsITS.splt01 <- data.frame(do.call('rbind', strsplit(as.character(NmsITS.Gn),'_',fixed=TRUE)))
# assign more descriptive column names
colnames(NmsITS.splt01) <- c("seqno","smplloc","smplyer","smplmnt","GenVar")

#Get unique and order elements by using a function
get.uniq.ord <-    function(cl){c <- unique(cl)
c <- c[order(c)]
return(c)}
# use this function on the columns with sample location
# and sample year
smplloc <- get.uniq.ord(NmsITS.splt01$smplloc)
smplyer <- get.uniq.ord(NmsITS.splt01$smplyer)
GenVar <- get.uniq.ord(NmsITS.splt01$GenVar)



# make the color categories for the location names
clfabNm <- rbind(loc4Lett ,colf.locNm)
colnames(clfabNm) <- clfabNm[1,]
clfabNm <-clfabNm[-1,]
mtx_clfabNm<- matrix(clfabNm)

# Get the number of gene variants
nGv <- length(GenVar)
# make a color function  based on the color palette
colfunc <- colorRampPalette(cbbPalette6)
colf.Gva <- colfunc(nGv)
# make the color categories for the names
clfGvar <- rbind(paste0("var",LETTERS[seq(1,nGv,1,)]) ,colf.Gva)
clfGvar <- rbind(paste0("var",(seq(1,nGv,1))) ,colf.Gva)
colnames(clfGvar) <- clfGvar[1,]
clfGvar <-clfGvar[-1,]
class(clfGvar)


cbbPalette5 <- c("black","firebrick2","brown","blue","cyan","white")
# make a color function  based on the color palette
nsmplYer <- length(smplyer)
colfunc <- colorRampPalette(cbbPalette5)
# 
colfsmplYer <- colfunc(nsmplYer)
colfsmplYer <- as.data.frame(t(colfsmplYer))
colfsmplYer <- as.character(colfsmplYer)
names(colfsmplYer) <- as.character(smplyer)
class(colfsmplYer)
#


#_______________________________________________________________________________
# section 04 -  end - make all disambiguity variants, and rename
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 05 -  start - make trees from disambiguity variants
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
  #plot(tre_pip)
  #sort the branches in the tree
  tre_pipr <- ape::ladderize(tre_pip, right = TRUE)
  #https://joey711.github.io/phyloseq/plot_tree-examples.html
  #plot(tre_pipr, cex=0.4)
  
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
  # get the number of unique genetic variants in the alignment 
  nGv <- length(unique(df_tiplb01$varNm))
  # make a color function  based on the color palette
  colfunc <- colorRampPalette(cbbPalette6)
  colf.Gva <- colfunc(nGv)
  colf.Gvafg <- colf.Gva[seq(1,nGv,1)]
  # make the color categories for the names
  clfGvar <- rbind(paste0("var",(seq(1,nGv,1))) ,colf.Gvafg)
  colnames(clfGvar) <- clfGvar[1,]
  clfGvar <-clfGvar[-1,]
  class(clfGvar)
  
  #________________
  # make a vector with the column names to color lables by as categories in the
  # NJ tree
  lblvars <- c("LocNm","varNm")
  #lblvars <- c("LocNm")
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
    clfltxtlb <- c(rep("black",(n.tpcats-no.white)),rep("white",no.white))
    # make the color categories a 
    tplbcls <- rbind(tpcats,clfltxtlb)
    colnames(tplbcls) <- tplbcls[1,]
    txtlblcols <-tplbcls[-1,]
    lblNms <- colnames(tplbcls)
    # make the tipcategories a dataframe
    dd <- as.data.frame(tipcategories)
    plt_titl_LET <- LETTERS[plt.idx.n]
    # get color range dependent on variable to plot
    if (l=="LocNm")
    {
      # check out color scales here: https://sjspielman.github.io/introverse/articles/color_fill_scales.html
      # match to get colors
      cf.Nm <- clfabNm[(match(lblNms,names(clfabNm)))]
    }
    if (l=="varNm")
    {
      # check out color scales here: https://sjspielman.github.io/introverse/articles/color_fill_scales.html
      # match to get colors
      cf.Nm <- clfGvar[(match(lblNms,names(clfGvar)))]
      
    }
    
    
    plt_n <- p %<+% dd + 
      ggtree::geom_tiplab(aes(fill = cat,
                              color= cat),
                          #color = "black", # color for label font
                          geom = "label",  # labels not text
                          label.padding = unit(0.09, "lines"), # amount of padding around the labels
                          label.size = 0, # size of label border
                          size=2.4) + 
      scale_fill_manual(values = cf.Nm) +
      #scale_fill_viridis_d(option = "viridis") +
      
      #scale_color_viridis_d(option = "G", direction = -1) +
      scale_color_manual(values = txtlblcols) +
      
      labs(fill=lg.title.rpl) +
      # hide one legend
      # https://statisticsglobe.com/remove-legend-ggplot2-r
      #guides(col = FALSE) +
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
  flNm.Fig <- paste0(wd00_wd05,"/Fig02_NJ_tree_disambig_seq_",GnNm,"_v01.png")
  # make an if test to check whether the plot should be saved
  # i.e. set to FALSE if you do not want the plot to be saved
  bSaveFigures=T
  # if the 'bSaveFigures' is TRUE the plot will be saved
  if(bSaveFigures==T){
    ggsave(plt_assmbl,file=flNm.Fig,
           width=210,height=297,
           #width=210,height=(297*0.5),
           units="mm",dpi=300)
  }
  #________________
  # end iteration over alignments
}
#_______________________________________________________________________________
# section 05 -  end - make trees from disambiguity variants
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 06 -  start - make haplotype networks
#_______________________________________________________________________________
# subset the alignments, to only comprise the Denmark and Germany
# samples
alvrs.M.ITS1_DG <- alvrs.M.ITS1[grepl("Mnelei",names(alvrs.M.ITS1))]
alvrs.M.ITS2_DG <- alvrs.M.ITS2[grepl("Mnelei",names(alvrs.M.ITS2))]
#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2,
                        alvrs.M.ITS1_DG,alvrs.M.ITS2_DG)
# also make a vector with the gene names
gene.nms <- c("ITS1_all","ITS2_all",
              "ITS1_DKGer","ITS2_DKGer")

# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)

# define different combinations of colors to reflect
# the labels that are to be assigned to the propotions
# in th pie diagrams -  these colors will used in the 
# iterations over gene names here below
cbbPalette1 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
cbbPalette3 <- c("white","yellow","orange","tomato","red","brown","black")
# cbbPalette3 <- c("firebrick4","firebrick2",#"orange",
#                  "gold3")
# As an alternative use this palette of grey hues
cbbPalette4 <- c("gray40","gray65",#"orange",
                 "gray80", "cyan1","cyan3","cyan4","deepskyblue4")
#iterate over the  sequence of numbers that  reflect the genes 
for (ng in n.f.alvrs.M)
{
  print(ng)
  #ng <- 1
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm<- gene.nms[ng]
  
  #copy DNAbin object again
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  # get the names for the sequences
  Nms.seq <- names(alvrsM)
  alvrsM <- bioseq::as_DNAbin(alvrsM)
  dnb_alg <- alvrsM 
  
  # make a distance matrix 
  dst_alg <- ape::dist.dna(dnb_alg, model= "raw")
  # make it a matrix
  mtx_pip <- as.matrix(dst_alg)
  # use the sequence names
  rownames(mtx_pip) <- Nms.seq
  # make it a data frame
  df_pip5 <- as.data.frame(mtx_pip)
  # split the row name string by a character 
  lpip5 <- strsplit(as.character(row.names(mtx_pip)), "_")
  # get the 2 nd element per vector in the list - this holds the abbreviated 
  # location name and the sampling year and the gene variant 
  lbl.locNm <- sapply(lpip5, "[[", 2)
  lbl.smplY <- sapply(lpip5, "[[", 3)
  lbl.genVa <- sapply(lpip5, "[[", 5)
  #make haplotype object
  ht4 <- pegas::haplotype(dnb_alg)
  #which polyps belong to which haplotype?
  ht4_indices<-attr(ht4, "index")
  # get number of labels 
  hlab4 <- length(labels(ht4))
  #assign labels
  names(ht4_indices)<-c(1:hlab4)
  # make roman numerals arabic numerals instead
  labs.ht4 <- as.numeric(as.roman(labels(ht4)))
  # use arabian numerals instead of roman numerals for haplotypes
  ht4 <- haplotype(dnb_alg,labels=c(labs.ht4))
  # make haplonet object
  hN4 <- pegas::haploNet(ht4)
  #prepare hpt table for location sampled
  # and for the year sampled
  ind.hap.locNm<-with(
    utils::stack(setNames(attr(ht4, "index"), rownames(ht4))),
    table(hap=ind, pop=lbl.locNm[values]))
  # for years sampled
  ind.hap.smplY<-with(
    utils::stack(setNames(attr(ht4, "index"), rownames(ht4))),
    table(hap=ind, pop=lbl.smplY[values]))
  # for gene variants sampled
  ind.hap.genVa<-with(
    utils::stack(setNames(attr(ht4, "index"), rownames(ht4))),
    table(hap=ind, pop=lbl.genVa[values]))
  
  # turn the haplotype tables into dataframes
  df_ihpt.locNm <- as.data.frame(ind.hap.locNm)
  df_ihpt.smplY <- as.data.frame(ind.hap.smplY)
  df_ihpt.genVa <- as.data.frame(ind.hap.genVa)
  # get the unique labels for the genetic variants, for location 
  # and for year sampled
  rL.locNm <- unique(df_ihpt.locNm$pop)
  rL.smplY <- unique(df_ihpt.smplY$pop)
  rL.genVa <- unique(df_ihpt.genVa$pop)
  # get the number of location categories
  lrL.locNm <- length(rL.locNm)
  lrL.smplY <- length(rL.smplY)
  lrL.genVa <- length(rL.genVa)
  
  # get the number of unique genetic variants in the alignment 
  nGv <- lrL.genVa
  # make a color function  based on the color palette
  colfunc <- colorRampPalette(cbbPalette6)
  colf.Gva <- colfunc(nGv)
  colf.Gvafg <- colf.Gva[seq(1,nGv,1)]
  # make the color categories for the names
  clfGvar <- rbind(paste0("var",(seq(1,nGv,1))) ,colf.Gvafg)
  colnames(clfGvar) <- clfGvar[1,]
  clfGvar <-clfGvar[-1,]
  
  
  colf.locNm2 <- clfabNm[match(rL.locNm,names(clfabNm))]
  colf.smplY2 <- colfsmplYer[match(rL.smplY,names(colfsmplYer))]
  colf.genVa2 <- clfGvar[match(rL.genVa,names(clfGvar))]
  
  
  # make a color range with 'colorRampPalette'
  colfunc1 <- colorRampPalette(cbbPalette1)
  colfunc2 <- colorRampPalette(cbbPalette5)
  colfunc3 <- colorRampPalette(cbbPalette6)
  #https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
  #for label categories
  colf.locNm <- colfunc1(lrL.locNm)
  colf.smplY <- colfunc2(lrL.smplY)
  colf.genVa <- colfunc3(lrL.genVa)
  # instead of using the color function to create a range of colors,
  # the colors corresponding to location and to gene variants
  colf.locNm <- clfabNm[match(rL.locNm,names(clfabNm))]
  colf.genVa <- clfGvar[match(rL.genVa,names(clfGvar))]
  # bind to dataframes
  df_colf.locNm <- as.data.frame(cbind(as.character(rL.locNm),colf.locNm))
  df_colf.smplY <- as.data.frame(cbind(as.character(rL.smplY),colf.smplY))
  df_colf.genVa <- as.data.frame(cbind(as.character(rL.genVa),colf.genVa))
  
  df_colf.locNm2 <- as.data.frame(cbind(as.character(rL.locNm),colf.locNm2))
  df_colf.smplY2 <- as.data.frame(cbind(as.character(rL.smplY),colf.smplY2))
  df_colf.genVa2 <- as.data.frame(cbind(as.character(rL.genVa),colf.genVa2))
  
  # rename the column headers
  colnames(df_colf.locNm) <- c("locNm","colf.locNm")
  colnames(df_colf.smplY) <- c("smplY","colf.smplY")
  colnames(df_colf.genVa) <- c("genVa","colf.genVa")
  
  colnames(df_colf.locNm2) <- c("locNm","colf.locNm")
  colnames(df_colf.smplY2) <- c("smplY","colf.smplY")
  colnames(df_colf.genVa2) <- c("genVa","colf.genVa")
  
  # match label names to get corresponding colors
  h_colf.locNm <- df_colf.locNm$colf.locNm[match(colnames(ind.hap.locNm),df_colf.locNm$locNm)]
  h_colf.smplY <- df_colf.smplY$colf.smplY[match(colnames(ind.hap.smplY),df_colf.smplY$smplY)]
  h_colf.genVa <- df_colf.genVa$colf.genVa[match(colnames(ind.hap.genVa),df_colf.genVa$genVa)]
  
  h_colf.locNm2 <- df_colf.locNm2$colf.locNm[match(colnames(ind.hap.locNm),df_colf.locNm2$locNm)]
  h_colf.smplY2 <- df_colf.smplY2$colf.smplY[match(colnames(ind.hap.smplY),df_colf.smplY2$smplY)]
  h_colf.genVa2 <- df_colf.genVa2$colf.genVa[match(colnames(ind.hap.genVa),df_colf.genVa2$genVa)]
  
  # match label names to get corresponding colors
  h_colf.locNm <- df_colf.locNm$colf.locNm[match(colnames(ind.hap.locNm),df_colf.locNm$locNm)]
  h_colf.smplY <- df_colf.smplY$colf.smplY[match(colnames(ind.hap.smplY),df_colf.smplY$smplY)]
  h_colf.genVa <- df_colf.genVa$colf.genVa[match(colnames(ind.hap.genVa),df_colf.genVa$genVa)]
  
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  
  # make a file name to store the plot in
  flnm <- c(paste("Fig03_v",pnng,"_haplotype_network_for_",GnNm,"_01.png",  sep = ""))
  #paste output directory and filename together in a string
  outflnm <- paste(wd00_wd05,"/",flnm,sep="")
  #  make the output file a 'png' file
  png(outflnm,
      width=(210),height=(297),res=150, units = "mm")
  #define plot arrangement
  tbt.par <- par(mfrow=c(3, 1),
                 oma=c(0,0,0,0), #define outer margins
                 mai=c(0,0,0,0), #define inner margins
                 mar=c(0,0,3,0))
  # to prevent the subfigure letter being overset by the next subfigure plot
  # then adjust the 'mar=c(0,0,2,0))' part to 'mar=c(0,0,3,0))'
  # or to 'mar=c(0,0,4,0))' - I worked this out by iterating over a sequence
  # of numbers and then storing the different plots and then inspecting
  # the resulting plots
  #______________ make plot a __________________________________________________
  #make a haplotype network
  #make the plot
  plot(hN4, 
       size = sqrt(attr(hN4,"freq")/pi), 
       #size = log10(attr(hN4,"freq")), 
       scale.ratio = 0.6, 
       cex = 1.1, # set size of roman numerals on circles for haplotype ID
       bg= h_colf.locNm2, 
       pie = ind.hap.locNm, 
       show.mutation = 2, threshold = 0, labels(T))
  
  #add a legend to the plot
  legend("topright",colnames(ind.hap.locNm), 
         pt.bg=h_colf.locNm2,box.col=NA,
         #col=rainbow(ncol(new.hap.smplloc)), 
         pt.lwd=0.4,
         pch=21, ncol=2, cex=1.4)
  title(main = "a",
        cex.main = 4.0,   font.main= 2, col.main= "black",
        adj = 0.01, line = 0.1)
  #______________ make plot b __________________________________________________
  #make the plot
  plot(hN4, 
       size = sqrt(attr(hN4,"freq")/pi), 
       #size = log10(attr(hN4,"freq")), 
       scale.ratio = 0.6, 
       cex = 1.1, # set size of roman numerals on circles for haplotype ID
       bg= h_colf.smplY2, 
       pie = ind.hap.smplY, 
       show.mutation = 2, threshold = 0, labels(T))
  
  #add a legend to the plot
  legend("topright",colnames(ind.hap.smplY), 
         #col=rainbow(ncol(new.hap.smplye)), 
         pt.bg=h_colf.smplY2,box.col=NA,pt.lwd=0.4,
         pch=21, ncol=1, cex=1.4)
  title(main = "b",
        cex.main = 4.0,   font.main= 2, col.main= "black",
        adj = 0.01, line = 0.1)
  #______________ make plot c __________________________________________________
  #make the plot
  plot(hN4, 
       size = sqrt(attr(hN4,"freq")/pi), 
       #size = log10(attr(hN4,"freq")), 
       scale.ratio = 0.6, 
       cex = 1.1, # set size of roman numerals on circles for haplotype ID
       bg= h_colf.genVa2, 
       pie = ind.hap.genVa, 
       show.mutation = 2, threshold = 0, labels(T))
  
  #add a legend to the plot
  legend("topright",colnames(ind.hap.genVa), 
         #col=rainbow(ncol(new.hap.smplye)), 
         pt.bg=h_colf.genVa2,box.col=NA,pt.lwd=0.4,
         pch=21, ncol=1, cex=1.4)
  title(main = "c",
        cex.main = 4.0,   font.main= 2, col.main= "black",
        adj = 0.01, line = 0.1)
  #______________ stop making plot a , b and c _________________________________
  par(tbt.par)
  # end svg file to save as
  dev.off()  
  #reset this parameter, to get the default traditional plot 
  # parameters back again
  par(mfrow = c(1, 1)) 
  # end the iteration over genes
}
#

#_______________________________________________________________________________
# section 06 -  end - make haplotype networks
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 07 -  start - make ML trees
#_______________________________________________________________________________


#_______________________________________________________________________________
# start  -  make ML tree with alignment included for all the samples, 
# including the NCBI Genbank obtained samples
#_______________________________________________________________________________

# I had some issue with 'phangorn::pml_bb' when one alignment 
# was to prepared as a ML tree, perhaps setting the seed 
# might help avoid this problem reoccuring
set.seed(124)

# https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html
library(phangorn)

# make a vector that categorizes the overall sampling locations for :
# North East Atlantic, Caspian Sea, Central Western Atlantic and Mediterranean
ov.cll <- c(
  "BalS","NEA",
  "CasS", "CS",
  "CWAt","CWA",
  "FBog","NEA",
  "FKer","NEA",
  "GBus","NEA",
  "GHel","NEA",
  "GKie","NEA",
  "GWis","NEA",
  "JLim","NEA",
  "JMar","NEA",
  "Medi","M",
  "Neth","NEA",
  "SBal","NEA",
  "SSko","NEA",
  "USAP","NWA",
  "USAW","NWA"
)
# arrange it as a data frame
df_ovl <- as.data.frame(matrix(ov.cll ,ncol=2, byrow=T))
# change the column headers
colnames(df_ovl) <- c("uloc", "ovloc")
# subset to only comprise the NEA overall location names
df_ovlNEA <- subset(df_ovl, ovloc=="NEA")
# get unique NEA locations and all full sample names
NEAloc <- df_ovlNEA$uloc
# this data frame prepared above will be used for identifying 
# the North East Atlantic samples (the NEA samples) in the 
#make a list to collect plots in
lst.MLplts  <- list()
# iteration below
# iterate over the numbers for the alignments in the list of alignments 
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
  dnb_alsM <- bioseq::as_DNAbin(alvrsM)
  names(dnb_alsM) <- LalvrsM
  # make the dnabin object a phydata object
  phyd_alsM <- phangorn::as.phyDat(dnb_alsM)
  # strsplit names and turn into characters
  # and rowbind the nested lists to a dataframe
  df_ppNm <- data.frame(do.call
                        ('rbind', 
                          strsplit(as.character(names(phyd_alsM)),
                                   "_")))
  # modify columns names
  colnames(df_ppNm) <- c("smplNm","loc","yer","mtn","GenVar")
  # replace the Mnelei abbr numbers with NCBI accession numbers obtained
  # after depositing sequences on NCBI
  # then get the names
  lbN <- df_ppNm$smplNm
  # get unique location names
  uloc <- unique(df_ppNm$loc)
  # order the unique location names
  uloc <- uloc[order(uloc)]
  
  # Get the names from the phydat object
  Nmsloc <- names(phyd_alsM)
  # use grepl instead of dplyr, to get all NEA full sample names
  # but first collapse and paste with '|' sign to be able to grep for
  # multiple chracteristics
  
  # https://stackoverflow.com/questions/45357806/dplyr-select-and-starts-with-on-multiple-values-in-a-variable-list?noredirect=1&lq=1
  NmslocNEA <- Nmsloc[grepl(paste(NEAloc, collapse="|"),Nmsloc)]
  # now subset the phydat object by the full names that match bein NEA samples
  phyd_pip <- subset(phyd_alsM, NmslocNEA)
  # make the alignment a data matrix
  df_phdp <- as.data.frame(as.matrix(phyd_alsM))
  # count rows and columns
  nrd4 <- nrow(df_phdp)
  ncd4 <- ncol(df_phdp)
  
  library(phangorn)
  # calculate distances
  dm4 <- phangorn::dist.p(phyd_alsM)
  #dm5 <- phangorn::dist.p(phyd_pip)
  
  dm4 <- phangorn::dist.ml(phyd_alsM, "F81")
  #dm5 <- phangorn::dist.ml(phyd_pip, "F81")
  tree_NJ4 <- NJ(dm4)
  # try the model test function in phangorn 
  #Modl.f.phyd <- phangorn::modelTest(phyd_alsM)
  # extract best model
  #(best_model <- as.pml(Modl.f.phyd))
  
  #(fit4 <- pml_bb(phyd_alsM, model="GTR+G"))
  # (fit4 <- phangorn::pml_bb(phyd_alsM, model="JC",
  #                 method	="unrooted" ) )#, "ultrametric" or "tipdated"))
  # Make a model fitted - this creates a 'pml'-object 
  fit_alsM <- phangorn::pml(tree_NJ4, phyd_alsM)
  #fit4 <- fit_alsM
  #prepare a "multiPhylo"-object
  #bs_alsM <- phangorn::bootstrap.pml(fit4, bs=100, optNni=F, multicore=F)
  bs_alsM <- phangorn::bootstrap.pml(fit_alsM, bs=100, optNni=F, multicore=F)
  # make the 'pml'-object and the "multiPhylo"-object a 'phylo object
  phytree_ml_alsM <- phangorn::plotBS(fit_alsM$tree, bs_alsM)
  # https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
  # also see : https://stackoverflow.com/questions/52956181/how-do-i-annotate-a-ggtree-object-with-bootstrap-values
  gt04 <- ggtree::ggtree(phytree_ml_alsM, right =T, #layout="slanted",
                         options(ignore.negative.edge=TRUE)) +
    ggtree::geom_tiplab(align=F, size=1.82) +
    ggtree::geom_nodelab(aes(label = label),color="red",hjust=-.3, size=4.2) +
    ggtree::geom_treescale() # adds the scale
  gt04
  # try the groupClade function, this function needs a number that
  # defines a group but it can also just be feed a number that is higher 
  # than the number of tips in the tree
  nfGcl <- length(phytree_ml_alsM$tip.label)+1  
  # specify the clade number for the long clade
  phytree_ml_alsM.wG <- tidytree::groupClade(phytree_ml_alsM, c(nfGcl))
  # Split the sequence names to get the variables to assign colors to
  # assign tip points by. This requires a match back to the 'clfabNm'
  # data frame with the colors for the abbreviated names and for the
  # sample locations
  # arrange the sequence names as a data frame
  df_Nms <- data.frame(do.call('rbind', strsplit(as.character(phytree_ml_alsM$tip.label),'_',fixed=TRUE)))
  colnames(df_Nms) <- c("smplNm","LocNm","smplYer","smplMn","varNm")
  
  df_Nms$lngseqNm <- phytree_ml_alsM$tip.label
  # substitute to get location name
  Nmsloc <-  df_Nms[,2]
  # match with data frame that has colours for locations 
  cr <- clfabNm[(match(Nmsloc,names(clfabNm)))]
  cr1 <- cr
  names(cr) <- phytree_ml_alsM$tip.label
  
  # From the 'phylo' object a 'tibble tree object' must be prepared
  # to be able to plot the tree with labels
  # see here: https://yulab-smu.top/treedata-book/chapter2.html
  # Under ' 2.1.1 The phylo object '
  # where the phylo tree obejct can be used for getting a 'tibble tree object'
  # with the information aobut the nodes
  tbt_ml_alsM.wG <- as_tibble(phytree_ml_alsM.wG)
  class(tbt_ml_alsM.wG) # it is a 'tibble tree object'
  # the ggtree plotting will not allow you to plot the 'tibble tree object'
  # use the 'tidytree::as.treedata' to convert it to a "treedata tidytree" object 
  tdtt_ml_alsM.wG <- tidytree::as.treedata(tbt_ml_alsM.wG)
  class(tdtt_ml_alsM.wG)
  
  
  # https://www.rdocumentation.org/packages/geiger/versions/2.0.10/topics/rescale
  # plot a tree with colors on the long clade
  p4.01 <- ggtree(phytree_ml_alsM.wG, right =T, aes(color=group)) + 
    ggtree::geom_tiplab(align=F, size=1.82) +
    geom_tippoint(size = 2, color = cr) +
    
    scale_color_manual(values=c("black", "firebrick", "steelblue"))
  p4.01
  
  
  # plot a tree without colors on the long clade
  p4.02 <- ggtree(tdtt_ml_alsM.wG, #it has to be the "treedata tidytree" object 
                  right =T, ladderize = T) + 
    geom_tippoint(aes(color=label),
                  size=2.2 ) +
    #?geom_tiplab(size=1.4, align=T) +
    scale_color_manual(values = cr) +
    guides(col = "none") +
    
    #ggtree::geom_text2(aes(label=node,x=branch), size=2.2) +
    ggtree::geom_text2(aes(label=round(as.numeric(label)),
                           subset=as.numeric(label)< 100, 
                           x=branch), vjust=0, color="black",
                       hjust=0.8, size=3.2) #+
  #scale_color_manual(values=c("black", "firebrick", "steelblue"))
  p4.02
  
  # replace all NAs with '-' since the 'ggtree::msaplot' function
  # cannot handle NAs in  the sequence
  alvrsM.1 <- alvrsM %>% replace(is.na(.), "-")
  #pad with zeros to three characters for  to get a running number 
  # to assign files
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  # make the 'bioseq alignment' a data frame instead
  df_alsM <- as.data.frame(as_tibble(alvrsM.1))
  # ensure there are no NAs in the sequences
  df_alsM[,1] <- df_alsM[,1] %>% replace(is.na(.), "-")
  # the data frame cannot be made a bioseq tibble directly
  # use only the columns with the sequences to make a bioseq
  # object with disambiguity sequences 
  bsq_disam <- bioseq::as_dna(df_alsM[,1])
  # then make this a DNAbin object - however this DNAbin object need names
  dnb_disam <- bioseq::as_DNAbin(bsq_disam)
  # get the names from column 1 in the dataframe
  names(dnb_disam) <- names(alvrsM.1)
  # define file name to write to using the zero-padded number
  outF.algn <- paste0(wd00_wd05,"/algn_Mnelei_pip_v",pnng,".fasta")
  #write the DNAbin object as a fasta file - this fasta file is required 
  # for the 'ggtree::msaplot' function
  ape::write.FASTA(dnb_disam,outF.algn)
  #https://stackoverflow.com/questions/68925906/set-x-axis-on-ggtree-heatmap-in-r
  # plot alignment together with tree, the offset defines the space between aligment and tree
  p4.03 <- ggtree::msaplot(p4.02, #offset=0.05, 
                           fasta=outF.algn)
  p4.03 <- p4.03 + scale_fill_manual(values=c(#"gray70",
    "violetred3",
    "blue",
    "orange",
    "green1"))
  p4.03 <- p4.03 + ggplot2::labs(fill='nt')
  
  lst.MLplts[[ng]] <- p4.03
  bSaveFigures=T
  # make an if test to check whether the plot should be saved
  # i.e. set to FALSE if you do not want the plot to be saved
  if(bSaveFigures==T){
    ggsave(p4.03,file=paste0(wd00_wd05,"/Fig04_v",pnng,"_MLtree_w_align_for_",GnNm,".png"),
           #width=210,height=297,
           width=210*0.5,height=(297),
           units="mm",dpi=300)
  }
  # end iteration over the numbers for the alignments in the list of alignments 
}



library(patchwork)
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p01t <- lst.MLplts[[1]] #+ labs(title = "a", face="bold")#,
p02t <- lst.MLplts[[2]] #+ labs(title = "b", face="bold")#,
#p03t <- p03 +   ggtitle('Plot 4', size =16) # (title = "c", face="bold")#,
#p03t <- p03 + theme(plot.title = element_text(size = 12, face = "bold"))
#p03t <- p03 + theme(plot.title = element_text(face="bold", size=18))

pA <-   p01t + 
        p02t +
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  # see : https://stackoverflow.com/questions/63434957/how-to-adjust-the-font-style-of-tags-with-plot-annotation-in-figures-assembled-w
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold', size =22)) +
  plot_layout(guides = "collect") #+
#plot_annotation(caption=pthinf01) #& theme(legend.position = "bottom")
#p
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig04_v05_MLtree_all_ITS1_and_ITS2.png")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,
         width=210,
         height=297,
         units="mm",dpi=300)
}

#_______________________________________________________________________________
# end  -  make ML tree with alignment included for all the samples, 
# including the NCBI Genbank obtained samples
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 07 -  end - make ML trees
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 08 -  start - plot genotypes on map  for Denmark, Germany and all NCBI
# samples
#_______________________________________________________________________________

library(mapplots)
library(maps)
library(mapdata)

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
  #}
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
  #_______________________________________________________________________________
  # start - plot the DNAbin object with all  haplotypes 
  # for both Denmark- Germany samples collected in 2017-2018 and 
  # NCBI GenBank samples
  #_______________________________________________________________________________
  # use the alignment with modified variants as DNAbin object
  dnb_pip4 <- alvrsM
  # make the dnabin object a phydata object
  phyd_pip4 <- phangorn::as.phyDat(dnb_pip4)
  
  # make the phy-dat object a DNAbin object instead
  dnb_pip7 <- as.DNAbin(phyd_pip4)
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
  
  #df_colf.locNm
  
  # get latitude and longtiude
  df_ihpt07$dec_lat <-  df_clo04$dec_lat2[match(df_ihpt07$pop.loc,df_clo04$loc4Lett)]
  df_ihpt07$dec_lon <-  df_clo04$dec_lon2[match(df_ihpt07$pop.loc,df_clo04$loc4Lett)]
  
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
  cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
  # make a colour ramp
  colfunc <- colorRampPalette(cbbPalette2)
  #https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
  # use the colour function
  cl07 <- colfunc(nHpt7)
  # Replace the colour range with rainbow colours
  colra7 <- cl07
  
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  
  # define name for output plot file
  figname01 <- paste0("Fig05_v",pnng,"_map_genotype_pie_",GnNm,".png")
  pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
  # set to save plot as pdf file with dimensions 8.26 to 2.9
  # 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
  # and 210 mm and 74.25 mm matches 1/4 of a A4 page
  png(pthfignm01,
      width = 297, 
      height = 210, 
      units="mm",
      res =300)#,width=(1.6*2.9),height=(0.8*8.26))
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
  lrgst_pie <- 1.1*sqrt(max(rowSums(xyz.Hpt7$z))/pi)
  #add pies to map
  mapplots::draw.pie(xyz.Hpt7$x, xyz.Hpt7$y, sqrt(xyz.Hpt7$z/pi), 
                     radius = lrgst_pie, col=scales::alpha(colra7,0.7))
  #_______add histogram bars to positions on the map
  # add a legend
  legend.z <- round(max(rowSums(xyz.Hpt7$z,na.rm=TRUE)),0)
  #write.csv2(as.data.frame(as.matrix(xyz.Hpt7$z)),paste0(wd00_wd05,"/xyz.Hpt7.csv"))
  # also adjust the circle size for the legend , so that it is based on the
  # square root of the max size divided by pi
  legend.bubble(-40,30,z=legend.z,round=0,
                maxradius=lrgst_pie,bty="n",txt.cex=1.1)
  text(-40,34.2,"samples",cex=1.1) 
  #add legend to plot
  legend("bottomright", inset=c(0.03, 0), legend=c(sHpt7), pch=c(rep(22,nHpt7)),
         bg="white",title=paste0("genotype ",GnNm), pt.bg=colra7, cex=0.96, pt.cex = 1.2, ncol=6,
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
  # modify the total of haplotype numbers
  tot.hp[1] <- "Total"
  # bind the row with the total
  df_i07 <- rbind(df_i07,tot.hp)
  # replcae 0 with nothing
  df_i07[(df_i07==0)] <- ""
  # rename the column header
  colnames(df_i07)[grepl("hap",colnames(df_i07))] <- "genotype"
  # write out a csv file
  write.csv(df_i07,
            file = paste0(wd00_wd05,"/Table02a_v",pnng,"_genotype_freq_global_",GnNm,".csv"))
  
  # end iteration over ITS gene variant alignments
}
#_______________________________________________________________________________
# end - plot the DNAbin object with all  haplotypes 
# for both Denmark- Germany samples collected in 2017-2018 and 
# NCBI GenBank samples
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 08 -  end - plot genotypes on map  for Denmark, Germany and all NCBI
# samples
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 09 -  start - plot genotypes on map  for Denmark, Germany  samples
#_______________________________________________________________________________

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
  #}
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm<- gene.nms[ng]
  # Now the alignment has been picked from the list,
  # Then prepare the NJ trees, and the label categories
  # and a color range to assign the labels in the trees
  # also get the labels on the sequences
  LalvrsM <- names(alvrsM)
  
  # strsplit names and turn into characters
  # and rowbind the nested lists to a dataframe
  df_ppNm <- data.frame(do.call
                        ('rbind', 
                          strsplit(as.character(names(alvrsM)),
                                   "_")))
  # modify columns names
  colnames(df_ppNm) <- c("smplNm","loc","yer","mtn","GenVar")
  # replace the Mnelei abbr numbers with NCBI accession numbers obtained
  # after depositing sequences on NCBI
  # then get the names
  lbN <- df_ppNm$smplNm
  # get unique location names
  uloc <- unique(df_ppNm$loc)
  # order the unique location names
  uloc <- uloc[order(uloc)]
  # Get the names from the phydat object
  Nmsloc <- names(alvrsM)
  # https://stackoverflow.com/questions/45357806/dplyr-select-and-starts-with-on-multiple-values-in-a-variable-list?noredirect=1&lq=1
  NmslocNEA <- Nmsloc[grepl(paste(NEAloc, collapse="|"),Nmsloc)]
  
  # now subset the bioseq object by the full names that match the NEA samples
  NmslocNEA.tgrp <- paste(NmslocNEA,collapse = "|")
  bsq_pip5 <- alvrsM[grepl(NmslocNEA.tgrp,names(alvrsM))]
  # limit to only include samples collected in 2017 and 2018
  Nmswy17_18 <- NmslocNEA[grepl(paste(c("2017","2018"), collapse="|"),NmslocNEA)]
  # collapse a string of sequences names to grep for
  Nmswy17_18.tgrp <- paste(Nmswy17_18,collapse = "|")
  bsq_pip5 <- bsq_pip5[grepl(Nmswy17_18.tgrp,names(bsq_pip5))]
  # get the names for the sequences
  Lbsq5 <- names(bsq_pip5)
  #make the bioseq object a DNAbin object
  dnb_alsM <-  bioseq::as_DNAbin(bsq_pip5)
  #also assign the names to the sequences
  names(dnb_alsM) <- Lbsq5
  dnb_pip5 <- dnb_alsM
  
  
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
  # use arabian numerals instead of roman numerals for genotypes
  ht5 <- pegas::haplotype(dnb_pip5,labels=c(labs.ht5))
  # make haplonet object
  hN5 <- pegas::haploNet(ht5)
  
  #names(dnb_pip5)
  #prepare hpt table
  ind.hap5 <- with(
    utils::stack(setNames(attr(ht5, "index"), rownames(ht5))),
    table(hap=ind, pop=names(dnb_pip5)[values]))
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
  df_ihpt05$dec_lat <-  df_clo04$dec_lat2[match(df_ihpt05$pop.loc,df_clo04$loc4Lett)]
  df_ihpt05$dec_lon <-  df_clo04$dec_lon2[match(df_ihpt05$pop.loc,df_clo04$loc4Lett)]
  
  #______________________________________________________________
  #___________start plot map genotype _________________
  #______________________________________________________________
  # make xyz list 
  xyz.Hpt5 <- mapplots::make.xyz(df_ihpt05$dec_lon,df_ihpt05$dec_lat,df_ihpt05$Freq,df_ihpt05$hap.ab)
  length(xyz.Hpt5$x)==length(xyz.Hpt5$y) & length(xyz.Hpt5$y)==nrow(xyz.Hpt5$z)
  #sort the genotypes increasingly by number
  sHpt5 <- unique(df_ihpt05$hap.ab)[order(unique(df_ihpt05$hap.ab))]
  # count number of genotypes
  nHpt5 <- length(unique(df_ihpt05$hap.ab))
  # make a colour range to reflect the number of genotypes
  colra5 <- rainbow(nHpt5)
  # try making a different colour range
  cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                   "yellow","white")
  cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
  # make a colour ramp
  colfunc <- colorRampPalette(cbbPalette2)
  #https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
  # use the colour function
  cl05 <- colfunc(nHpt5)
  # Replace the colour range with rainbow colours
  colra5 <- cl05
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  
  # define name for output plot file
  figname01 <- paste0("Fig06_v",pnng,"_map_genotype_pie_",GnNm,".png")
  pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
  # set to save plot as pdf file with dimensions 8.26 to 2.9
  # 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
  # and 210 mm and 74.25 mm matches 1/4 of a A4 page
  png(pthfignm01,
      width = 297, 
      height = 210,
      units= "mm",
      res =300)#,width=(1.6*2.9),height=(0.8*8.26))
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
  lrgst_pie <- 0.2*sqrt(max(rowSums(xyz.Hpt5$z))/pi)
  #add pies to map
  draw.pie(xyz.Hpt5$x, xyz.Hpt5$y, sqrt(xyz.Hpt5$z/pi), 
           radius = lrgst_pie, col=scales::alpha(colra5,0.7))
  #_______add histogram bars to positions on the map
  # add a legend
  legend.z <- round(max(rowSums(xyz.Hpt5$z,na.rm=TRUE)),0)
  legend.bubble(5,54.4,legend.z,round=0,
                maxradius=lrgst_pie,bty="n",txt.cex=1.1)
  text(5,55.4,"samples",cex=1.1) 
  #add legend to plot
  legend("topright", inset=c(0.03, 0), legend=c(sHpt5), pch=c(rep(22,nHpt5)), 
         bg="white",title=paste0("genotype \n",GnNm), pt.bg=colra5, cex=0.96, pt.cex = 1.2, ncol=3,
         x.intersp = 0.6,y.intersp = 0.9)
  
  # end plot
  dev.off()
  #______________________________________________________________
  #___________end plot map genotype __________________
  #______________________________________________________________
  #define columns to keep
  cke <- c("hap","Freq","pop.loc")
  df_i05 <- df_ihpt05[cke]
  # summarize genotypes
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
  # rename the column header
  colnames(df_i05)[grepl("hap",colnames(df_i05))] <- "genotype"
  # write out a csv file
  write.csv(df_i05,
            file = paste0(wd00_wd05,"/Table02b_v",pnng,"_genotypes_freq_",GnNm,"_DK_Germany.csv"))
  
  # end iteration over alignments
}

#_______________________________________________________________________________
# section 09 -  end - plot genotypes on map  for Denmark, Germany  samples
#_______________________________________________________________________________



#_______________________________________________________________________________
# section 10 -  start - make Fst tables
#_______________________________________________________________________________

library(ape)
library(adegenet)
library(dartR)
library("StAMPP")

nms.algnm <- names(alvrs.M.ITS2)
alM.ITS2.DG <- alvrs.M.ITS2[grepl("Mnelei",nms.algnm)]

nms.algnm <- names(alvrs.M.ITS1)
alM.ITS1.DG <- alvrs.M.ITS1[grepl("Mnelei",nms.algnm)]

#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
#lst_alvrs.M.ITS <- list(alM.ITS1.DG,alM.ITS2.DG)
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
  
  #}
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm<- gene.nms[ng]
  # Now the alignment has been picked from the list,
  # Then prepare the NJ trees, and the label categories
  # and a color range to assign the labels in the trees
  # also get the labels on the sequences
  Lbsq5 <- names(alvrsM)
  #make the bioseq object a DNAbin object
  dnb_alsM <-  bioseq::as_DNAbin(alvrsM)
  #also assign the names to the sequences
  names(dnb_alsM) <- Lbsq5
  dnb_pip5 <- dnb_alsM
  dnb_02 <- dnb_pip5
  #copy the DNAbin object again
  dnb_pip02 <- dnb_02
  # check the label names
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
  # use 'StAMPP' package to calculate Fst and p-values on populations
  Mnl.fst <- StAMPP::stamppFst(sta_mnelei, 100, 95, 1)
  
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
  capt_tbl02 <- paste0("Table 3. Comparison of FST for Mnemiopsis leidyi obtained from sequences of nDNA ",GnNm," for the sampling locations. The lower triangle shows the Fst values, the corresponding probability values are in the upper triangle. Sequences of nDNA ",GnNm," obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A low FST index indicates no variation among the individuals sampled, a high FST index indicates there is high variation among sequences compared. The color gradient goes from yellow for low FST values, and dark blue for high FST values.")
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
  #
  #t.HTML03
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  
  # ensure the package 'htmltools' is loaded to be able to 
  # save the html table
  library(htmltools)
  # Make a filename to store the html table with Fst vaules
  filNm.for_html <- paste0(wd00_wd05,"/Table03_v",pnng,"_html_table_Fst_val_for_",GnNm,".html")
  
  htmltools::save_html(t.HTML03,file=filNm.for_html)
  # end iteration over alignments in list
}

#_______________________________________________________________________________
# section 10 -  end - make Fst tables
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 11 -  start - plot genotypes on map  for Denmark, Germany  samples
#_______________________________________________________________________________

#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)


# #pad with zeros to three characters for own Mnelei samples
# #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
# pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)

# define name for output plot file
figname01 <- paste0("Fig07_v01_map_genotype_pie_DK_Germ.png")
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
png(pthfignm01,
    width = 210, 
    height = 297*0.36,
    units= "mm",
    res =300)#,width=(1.6*2.9),height=(0.8*8.26))
#add extra space to the right of the plot
# the arrangement is c('bottom,'left','top','right')
par(mar=c(4.2, 4.2, 1.8, 0.4), xpd=FALSE)
par(oma=c(0, 0, 0.2, 0))
#reset this parameter
par(mfrow = c(1, 2)) 

#iterate over the  sequence of numbers that  reflect the genes 
for (ng in n.f.alvrs.M)
{
  print(ng)
  #}
  subfglts <- letters[n.f.alvrs.M]
  sbfigl <- subfglts[ng]
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm<- gene.nms[ng]
  # Now the alignment has been picked from the list,
  # Then prepare the NJ trees, and the label categories
  # and a color range to assign the labels in the trees
  # also get the labels on the sequences
  LalvrsM <- names(alvrsM)
  
  # strsplit names and turn into characters
  # and rowbind the nested lists to a dataframe
  df_ppNm <- data.frame(do.call
                        ('rbind', 
                          strsplit(as.character(names(alvrsM)),
                                   "_")))
  # modify columns names
  colnames(df_ppNm) <- c("smplNm","loc","yer","mtn","GenVar")
  # replace the Mnelei abbr numbers with NCBI accession numbers obtained
  # after depositing sequences on NCBI
  # then get the names
  lbN <- df_ppNm$smplNm
  # get unique location names
  uloc <- unique(df_ppNm$loc)
  # order the unique location names
  uloc <- uloc[order(uloc)]
  # Get the names from the phydat object
  Nmsloc <- names(alvrsM)
  # https://stackoverflow.com/questions/45357806/dplyr-select-and-starts-with-on-multiple-values-in-a-variable-list?noredirect=1&lq=1
  NmslocNEA <- Nmsloc[grepl(paste(NEAloc, collapse="|"),Nmsloc)]
  
  # now subset the bioseq object by the full names that match the NEA samples
  NmslocNEA.tgrp <- paste(NmslocNEA,collapse = "|")
  bsq_pip5 <- alvrsM[grepl(NmslocNEA.tgrp,names(alvrsM))]
  # limit to only include samples collected in 2017 and 2018
  Nmswy17_18 <- NmslocNEA[grepl(paste(c("2017","2018"), collapse="|"),NmslocNEA)]
  # collapse a string of sequences names to grep for
  Nmswy17_18.tgrp <- paste(Nmswy17_18,collapse = "|")
  bsq_pip5 <- bsq_pip5[grepl(Nmswy17_18.tgrp,names(bsq_pip5))]
  # get the names for the sequences
  Lbsq5 <- names(bsq_pip5)
  #make the bioseq object a DNAbin object
  dnb_alsM <-  bioseq::as_DNAbin(bsq_pip5)
  #also assign the names to the sequences
  names(dnb_alsM) <- Lbsq5
  dnb_pip5 <- dnb_alsM
  
  
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
  # use arabian numerals instead of roman numerals for genotypes
  ht5 <- pegas::haplotype(dnb_pip5,labels=c(labs.ht5))
  # make haplonet object
  hN5 <- pegas::haploNet(ht5)
  
  #names(dnb_pip5)
  #prepare hpt table
  ind.hap5 <- with(
    utils::stack(setNames(attr(ht5, "index"), rownames(ht5))),
    table(hap=ind, pop=names(dnb_pip5)[values]))
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
  df_ihpt05$dec_lat <-  df_clo04$dec_lat2[match(df_ihpt05$pop.loc,df_clo04$loc4Lett)]
  df_ihpt05$dec_lon <-  df_clo04$dec_lon2[match(df_ihpt05$pop.loc,df_clo04$loc4Lett)]
  
  #______________________________________________________________
  #___________start plot map genotype _________________
  #______________________________________________________________
  # make xyz list 
  xyz.Hpt5 <- mapplots::make.xyz(df_ihpt05$dec_lon,df_ihpt05$dec_lat,df_ihpt05$Freq,df_ihpt05$hap.ab)
  length(xyz.Hpt5$x)==length(xyz.Hpt5$y) & length(xyz.Hpt5$y)==nrow(xyz.Hpt5$z)
  #sort the genotypes increasingly by number
  sHpt5 <- unique(df_ihpt05$hap.ab)[order(unique(df_ihpt05$hap.ab))]
  # count number of genotypes
  nHpt5 <- length(unique(df_ihpt05$hap.ab))
  # make a colour range to reflect the number of genotypes
  colra5 <- rainbow(nHpt5)
  # try making a different colour range
  cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                   "yellow","white")
  cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
  # make a colour ramp
  colfunc <- colorRampPalette(cbbPalette2)
  #https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
  # use the colour function
  cl05 <- colfunc(nHpt5)
  # Replace the colour range with rainbow colours
  colra5 <- cl05
  
  # PLOT A
  # begin plot, with defined borders
  plot(NA,NA, xlim = c(4, 16), ylim = c(53, 59),
       las=1,
       #      # use 'asp' to change the aspect ratio: https://statisticsglobe.com/asp-r-plot
       asp=1.6,
       #surpress tickmarks
       xaxt="n", yaxt="n",
       #xlab="Longitude", ylab="Latitude",
       xlab="", ylab="",
       cex.lab = 0.6)
  title(xlab = "Longitude", ylab = "Latitude", line = 3.2,cex=2)
  
  title(main = sbfigl,
        cex.main = 2.0,   font.main= 2, col.main= "black",
        adj = 0.01, line = 0.4)
  # add internal map inside coastline
  maps::map('worldHires', add=TRUE, fill=TRUE, 
            xlim = c(4, 16), ylim = c(53, 59),
            xlab = "Longitude", ylab = "Latitude",
            col="azure3", bg=transp_col)
  #https://stackoverflow.com/questions/51799118/writing-the-symbol-degrees-celsius-in-axis-titles-with-r-plotly
  longwE <- paste(5*(0:4),"\u00B0E",sep="")
  axis(1, at=c(5*(0:4)), labels=longwE,line=0,cex.axis = 1.0)
  lattwN <- paste((53:59),"\u00B0N",sep="")
  axis(2, at=c(53:59), labels=lattwN, las=1,line=0,cex.axis = 1.0)
  #add pies to map
  # the 'draw.pie' function requires the largest pie diameter as input value for
  # the radius of the largest pie
  lrgst_pie <- 0.2*sqrt(max(rowSums(xyz.Hpt5$z))/pi)
  #add pies to map
  draw.pie(xyz.Hpt5$x, xyz.Hpt5$y, sqrt(xyz.Hpt5$z/pi), 
           radius = lrgst_pie, col=scales::alpha(colra5,0.7))
  #_______add histogram bars to positions on the map
  # add a legend
  legend.z <- round(max(rowSums(xyz.Hpt5$z,na.rm=TRUE)),0)
  legend.bubble(5,54.4,legend.z,round=0,
                maxradius=lrgst_pie,bty="n",txt.cex=0.82)
  text(5,55.5,"samples",cex=1) 
  #add legend to plot
  legend("topright", inset=c(0.03, 0), legend=c(sHpt5), pch=c(rep(22,nHpt5)), 
         bg="white",title=paste0("genotype \n",GnNm), 
         pt.bg=colra5, cex=0.8, pt.cex = 1.2, ncol=1,
         x.intersp = 0.6,y.intersp = 0.7)
  
  #______________________________________________________________
  #___________end plot map genotype __________________
  #______________________________________________________________
  # end iteration over alignments
}

# end plot
dev.off()

#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
# section 11 -  end - plot genotypes on map  for Denmark, Germany  samples
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 12 -  start - plot genotypes on map  for all  samples
#_______________________________________________________________________________

library(mapplots)
library(maps)
library(mapdata)

#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)



# define name for output plot file
figname01 <- paste0("Fig08_v01_map_genotype_pie_all.png")
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
png(pthfignm01,
    width = 210, 
    height = 297*0.33, 
    units="mm",
    res =300)#,width=(1.6*2.9),height=(0.8*8.26))
#add extra space to the right of the plot
# the arrangement is c('bottom,'left','top','right')
par(mar=c(6.2, 6.2, 1.8, 0.4), xpd=FALSE)
par(oma=c(0, 0, 0.2, 0))
#reset this parameter
par(mfrow = c(1, 2)) 

#iterate over the  sequence of numbers that  reflect the genes 
for (ng in n.f.alvrs.M)
{
  print(ng)
  subfglts <- letters[n.f.alvrs.M]
  sbfigl <- subfglts[ng]
  #}
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
  #_______________________________________________________________________________
  # start - plot the DNAbin object with all  haplotypes 
  # for both Denmark- Germany samples collected in 2017-2018 and 
  # NCBI GenBank samples
  #_______________________________________________________________________________
  # use the alignment with modified variants as DNAbin object
  dnb_pip4 <- alvrsM
  # make the dnabin object a phydata object
  phyd_pip4 <- phangorn::as.phyDat(dnb_pip4)
  
  # make the phy-dat object a DNAbin object instead
  dnb_pip7 <- as.DNAbin(phyd_pip4)
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
  
  #df_colf.locNm
  
  # get latitude and longtiude
  df_ihpt07$dec_lat <-  df_clo04$dec_lat2[match(df_ihpt07$pop.loc,df_clo04$loc4Lett)]
  df_ihpt07$dec_lon <-  df_clo04$dec_lon2[match(df_ihpt07$pop.loc,df_clo04$loc4Lett)]
  
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
  cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
  # make a colour ramp
  colfunc <- colorRampPalette(cbbPalette2)
  #https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
  # use the colour function
  cl07 <- colfunc(nHpt7)
  # Replace the colour range with rainbow colours
  colra7 <- cl07
  
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  
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
  title(main = sbfigl,
        cex.main = 1.6,   font.main= 2, col.main= "black",
        adj = 0.01, line = 0.4)
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
  
  
  axis(1, at=c(seq(-80,60, by=20)), labels=longwW0E,cex.axis = 0.8)
  lattwN <- paste(round(seq(10, 60,10),0),"\u00B0N",sep="")
  axis(2, at=c(round(seq(6, 60,10),0)), labels=lattwN, las=1,cex.axis = 0.8)
  # the 'draw.pie' function requires the largest pie diameter as input value for
  # the radius of the largest pie
  lrgst_pie <- 2*sqrt(max(rowSums(xyz.Hpt7$z))/pi)
  #add pies to map
  mapplots::draw.pie(xyz.Hpt7$x, xyz.Hpt7$y, sqrt(xyz.Hpt7$z/pi), 
                     radius = lrgst_pie, col=scales::alpha(colra7,0.7))
  #_______add histogram bars to positions on the map
  # add a legend
  legend.z <- round(max(rowSums(xyz.Hpt7$z,na.rm=TRUE)),0)
  #write.csv2(as.data.frame(as.matrix(xyz.Hpt7$z)),paste0(wd00_wd05,"/xyz.Hpt7.csv"))
  # also adjust the circle size for the legend , so that it is based on the
  # square root of the max size divided by pi
  legend.bubble(-40,30,z=legend.z,round=0,
                maxradius=lrgst_pie,bty="n",txt.cex=0.8)
  text(-40,38.2,"samples",cex=0.8) 
  #add legend to plot
  legend("bottomright", inset=c(0.03, 0), legend=c(sHpt7), pch=c(rep(22,nHpt7)),
         bg="white",title=paste0("genotype ",GnNm), pt.bg=colra7, cex=0.7, pt.cex = 1.2, ncol=6,
         x.intersp = 0.6,y.intersp = 0.7)
  
  
  # end iteration over ITS gene variant alignments
}


# end plot
dev.off()



#reset this parameter
par(mfrow = c(1, 1)) 


#_______________________________________________________________________________
# section 12 -  end - plot genotypes on map  for all  samples
#_______________________________________________________________________________
#_______________________________________________________________________________
# section 13 -  start - ANOSIM plot  for all  samples
#_______________________________________________________________________________

library(vegan)
#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)
# turn off previous plots
dev.off()
# plot the test
# define name for outpu plot file
figname08 <- "Fig09_v01_ANOSIMplot_location_all.jpg"
pthfignm08 <- paste(wd00_wd05,"/",figname08,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
png(pthfignm08,
    width = 210, 
    height = 297*0.4,
    units = "mm",
    res =300)#,width=(1.6*2.9),height=(0.8*8.26))
# the arrangement is c('bottom,'left','top','right')
par(mar=c(4, 5, 2, 1), xpd=FALSE,
    oma=c(0,0,1,0),
    mfrow = c(1, 2))
# iterate over numbers representing sequence alignments stored in a list
for (ng in n.f.alvrs.M)
{
  print(ng)
  #}
  subfglts <- letters[n.f.alvrs.M]
  sbfigl <- subfglts[ng]
  #}
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
  # use the alignment with modified variants as DNAbin object
  dnb_pip4 <- alvrsM
  
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
  # 'grp.loc5' needs to hold location names
  grp.loc5 <- sapply(lpip5, "[[", 2)
  # order locations alphabetically
  grp.l5 <- unique(grp.loc5)[order(unique(grp.loc5))]
  # make a sequence of numbers to use as index numbers for the locations
  nof.l5 <- seq(1,length(grp.l5))
  # bind the columns in a data frame that has index numbers per location
  df_nfl5 <- as.data.frame(cbind(nof.l5,grp.l5))
  # attach the group location back to the data frame
  df_pip5$grp.loc <- grp.loc5
  # match the location name to get an index number, and ensure this number is numeric
  # each location needs to have number assigned, in order to make the ANOSIM test work
  df_pip5$grp.loc <- as.numeric(df_nfl5$nof.l5[match(df_pip5$grp.loc,df_nfl5$grp.l5)])
  #make community matrix
  m_comp5 <- sapply(df_pip5[,1:ncol(df_pip5)],as.numeric)
  # ensure the row names in the matrix are as in the data frame
  row.names(m_comp5) <- row.names(df_pip5)
  # group by site
  df_comp5 <- as.data.frame(m_comp5)
  # remove all non duplicated rows - the ANOSIM test requires that there are
  # multiple representations per group
  df_comp5 <- df_comp5[duplicated(df_comp5[c('grp.loc')]), ]
  #row.names(df_comp5)
  slpip5 <- strsplit(as.character(row.names(df_comp5)), "_")
  # get second elements
  grp5 <- as.factor(sapply(slpip5, "[[", 2))
  # make it a matrix
  m_comp5 <- as.matrix(df_comp5)
  # make the ANOSIM test
  pip5.ano <- vegan::anosim(m_comp5,grp5)
  # plot the figure
  plot(pip5.ano,
       xlab="Location",
       ylab="Dissimilarity ranks between \n and within classes",
       las=2) # 'las=2' rotates the axis tick labels
  # add the corresponding subfigure letter
  title(main = sbfigl,
        cex.main = 1.6,   font.main= 2, col.main= "black",
        adj = 0.01, line = 0.4)
}
# end plot
dev.off()

#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
# section 13 -  end - ANOSIM plot  for all  samples
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 14 -  start - make PCA plots - all samples
#_______________________________________________________________________________


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

#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)
# turn off previous plots

# iterate over numbers representing sequence alignments stored in a list
for (ng in n.f.alvrs.M)
{
  print(ng)
  # }
  # ng <- 1
  # 
  subfglts <- letters[n.f.alvrs.M]
  sbfigl <- subfglts[ng]
  #}
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm<- gene.nms[ng]
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  
  # Now the alignment has been picked from the list,
  # Then prepare the NJ trees, and the label categories
  # and a color range to assign the labels in the trees
  # also get the labels on the sequences
  LalvrsM <- names(alvrsM)
  alvrsM <- bioseq::as_DNAbin(alvrsM)
  names(alvrsM) <- LalvrsM
  # use the alignment with modified variants as DNAbin object
  dnb_pip4 <- alvrsM
  
  #Make the dnabin a genind object
  gei_pip3 <- adegenet::DNAbin2genind(dnb_pip4)
  # convert genind object to genlight object with the dartR package
  gl_pip <- dartR::gi2gl(gei_pip3)
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
  # split the row name string by a character 
  h.pca.scrsplt <- strsplit(as.character(row.names(h.pca.scores)), "_")
  # get the 2 nd element per vector in the list - this holds the abbreviated 
  # location name, the 3rd element is the sampling year, the 5th is the 
  # gene variant
  poploc <- sapply(h.pca.scrsplt, "[[", 2)
  popyear <- sapply(h.pca.scrsplt, "[[", 3)
  popGvar <- sapply(h.pca.scrsplt, "[[", 5)
  
  # add a column with location names
  h.pca.scores$pop <- poploc
  # order the data frame by the column with location names
  h.pca.scores <- h.pca.scores[order(h.pca.scores$pop),]
  # get unique location names
  ppll3 <- unique(h.pca.scores$pop)
  
  library(ggplot2)
  set.seed(9)
  #make a plot with pca 
  p <- ggplot(h.pca.scores, aes(x=PC1, y=PC2, colour=pop, fill=pop)) 
  p <- p + geom_point(size=2,shape=21, color="black")
  p <- p + stat_ellipse(level = 0.95, linewidth = 1)
  p <- p + geom_hline(yintercept = 0) 
  p <- p + geom_vline(xintercept = 0) 
  p <- p + theme_bw()
  p <- p + guides(fill=guide_legend(ncol=2))
  p01 <- p
  #p01
  # color range from haplotype network pie diagrams for
  # location sampled - used in Fig.2.
  p01 <- p01 + scale_color_manual(values = clfabNm) 
  p01 <- p01 + scale_fill_manual(values = clfabNm) 
  # change the heading for the legend, this must be done for all 
  #settings for the points
  p01 <- p01 + labs(color='Location')
  p01 <- p01 + labs(fill='Location')
  p01 <- p01 + labs(shape='Location')
  p01
  
  # replace the pop name column in the pca data frame with the sampling year
  h.pca.scores$pop <- as.character(popyear)
  #h.pca.scores$pop <- poploc
  clfyr01 <- colfsmplYer[match(get.uniq.ord(h.pca.scores$pop),
                               names(colfsmplYer))]
  clfyr01 <- as.character(unlist(clfyr01))
  
  library(ggplot2)
  set.seed(9)
  p <- ggplot(h.pca.scores, aes(x=PC1, y=PC2, 
                                colour=pop,
                                fill=pop)) 
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
  #
  # set the color of the fill for lines and the fill for inside fill for the points
  p02 <- p02 + scale_color_manual(values = clfyr01) 
  p02 <- p02 + scale_fill_manual(values = clfyr01) 
  
  p02 <- p02 + labs(color='Sampling year')
  p02 <- p02 + labs(fill='Sampling year')
  p02 <- p02 + labs(shape='Sampling year')
  p02
  
  #_____
  # replace the pop name column in the pca data frame with the sampling year
  h.pca.scores$pop <- as.character(popGvar)
  
  #h.pca.scores$pop <- poploc
  
  library(ggplot2)
  set.seed(9)
  p <- ggplot(h.pca.scores, aes(x=PC1, y=PC2, 
                                colour=pop,
                                fill=pop)) 
  #p <- p + geom_point(size=2)
  p <- p + geom_point(size=2,shape=21, color="black")
  p <- p + stat_ellipse(level = 0.95, size = 1)
  #p <- p + scale_color_manual(values = cols) 
  p <- p + geom_hline(yintercept = 0) 
  p <- p + geom_vline(xintercept = 0) 
  p <- p + theme_bw()
  p03 <- p
  # color range from haplotype network pie diagrams for
  # locationyear sampled - used in Fig.2.
  #
  # get the number of unique genetic variants in the alignment 
  nGv <- length(unique(h.pca.scores$pop))
  # make a color function  based on the color palette
  colfunc <- colorRampPalette(cbbPalette6)
  colf.Gva <- colfunc(nGv)
  colf.Gvafg <- colf.Gva[seq(1,nGv,1)]
  # make the color categories for the names
  clfGvar <- rbind(paste0("var",(seq(1,nGv,1))) ,colf.Gvafg)
  colnames(clfGvar) <- clfGvar[1,]
  clfGvar <-clfGvar[-1,]
  
  # set the color of the fill for lines and the fill for inside fill for the points
  p03 <- p03 + scale_color_manual(values = clfGvar) 
  p03 <- p03 + scale_fill_manual(values = clfGvar) 
  
  p03 <- p03 + labs(color='Gene variant')
  p03 <- p03 + labs(fill='Gene variant')
  p03 <- p03 + labs(shape='Gene variant')
  p03
  
  
  library(patchwork)
  # see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
  #caption = "Data source: ToothGrowth")
  p01t <- p01 + labs(title = "a", face="bold", linewidth=14)#,
  p02t <- p02 + labs(title = "b", face="bold")#,
  #p03t <- p03 +   ggtitle('Plot 4', size =16) # (title = "c", face="bold")#,
  #p03t <- p03 + theme(plot.title = element_text(size = 12, face = "bold"))
  p03t <- p03 + theme(plot.title = element_text(face="bold", size=18))
  
  p01t <- p01
  p02t <- p02
  p03t <- p03
  
  pA <-  p01t +
    p02t +
    p03t +
    plot_layout(nrow=3,byrow=T) + #xlab(xlabel) +
    # see : https://stackoverflow.com/questions/63434957/how-to-adjust-the-font-style-of-tags-with-plot-annotation-in-figures-assembled-w
    plot_annotation(tag_levels = 'a') &
    theme(plot.tag = element_text(face = 'bold', size =22)) +
    
    
    plot_layout(guides = "collect") #+
  #plot_annotation(caption=pthinf01) #& theme(legend.position = "bottom")
  #p
  bSaveFigures=T
  #make filename to save plot to
  figname01 <- paste0("Fig10_v",pnng,"_pca_M_leyidy_",GnNm,".png")
  figname02 <- paste(wd00_wd05,"/",figname01,sep="")
  if(bSaveFigures==T){
    ggsave(pA,file=figname02,
           width=210*0.9,
           height=297,
           units="mm",dpi=300)
  }
  # end iteration over alignments
}
#

#_______________________________________________________________________________
# section 14 -  end - make PCA plots - all samples
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 15 -  start - Make haplotype circles on map second attempt 
# # on all samples
#_______________________________________________________________________________
# #install.packages("scatterpie")
# if(!require(scatterpie)){
#   install.packages("scatterpie")
# }

library(scatterpie)

#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)

lst_plts.G <-list()
#iterate over the  sequence of numbers that  reflect the genes 
for (ng in n.f.alvrs.M)
{
  print(ng)
  #}
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
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  # make a name for the legend in the ggplot
  lgn.Nm <- paste0("genotype ",GnNm)
  # use the alignment with modified variants as DNAbin object
  dnb_pip4 <- alvrsM
  pip <- dnb_pip4
  #create haplotypes from dna.bin
  pipHaps <- pegas::haplotype(pip)
  #prepare hpt table
  ih<-with(
    utils::stack(setNames(attr(pipHaps, "index"), rownames(pipHaps))),
    table(hap=ind, pop=names(pip)[values]))
  # make it a data frame
  dih01 <- as.data.frame(ih)
  #limit to include only 'Freq' that equals 1
  dih02 <- dih01[dih01$Freq == 1,]	
  #get pop groups from row names
  # split the row name string by a character 
  dihpg <- strsplit(as.character(dih02$pop), "_")
  # get the 2 nd element per vector in the list - this holds the abbreviated 
  # location name, the 3rd element is the sampling year, the 5th is the 
  # gene variant
  poploc <- sapply(dihpg, "[[", 2)
  popyear <- sapply(dihpg, "[[", 3)
  popGvar <- sapply(dihpg, "[[", 5)
  #make it a table
  hsll3 <- table(dih02$hap, poploc)
  # make the roman numerals arabic numerals
  rwNm.hsll3 <- as.numeric(as.roman(rownames(hsll3)))
  rwNm.hsll3 <- paste0("G",rwNm.hsll3)
  # get the names
  rownames(hsll3) <- rwNm.hsll3
  # transpose
  thl1.3 <- t(hsll3)
  # make it a data frame
  dhl01 <- as.data.frame(thl1.3)
  # reshape the data frame for long to wide
  dhl02 <- reshape(data=dhl01,idvar="poploc",
                   v.names = "Freq",
                   timevar = "Var2",
                   direction="wide")
  enc <- ncol(dhl02)
  # substitute in the column names
  clNdh02 <- colnames(dhl02)
  clNdh02 <- gsub("Freq\\.","",clNdh02)
  colnames(dhl02) <- clNdh02
  #make a viridis colour range
  cl03 <- pals::viridis(length(unique(dhl02[,c(2:enc)])))
  cl03 <- pals::inferno(length(unique(dhl02[,c(2:enc)])))
  ncolHptloc03 <- length(unique(colnames(dhl02)[2:enc]))
  cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                   "yellow","white")
  cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
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
  # sum up rows , avoid confusion with function 'rowsum'
  dhl02$rws <- rowSums( dhl02[,c(2:enc)])
  # Get latitude and longitude positions
  dhl02$dec_lat <- df_clo04$dec_lat2[match(dhl02$poploc, df_clo04$loc4Lett)]
  dhl02$dec_lon <- df_clo04$dec_lon2[match(dhl02$poploc, df_clo04$loc4Lett)]
  dhl02$dec_lat <- as.numeric(dhl02$dec_lat)
  dhl02$dec_lon <- as.numeric(dhl02$dec_lon)
  # get limits for the latitudes and longitudes to use in the ggplot
  mnlat <- floor(min(dhl02$dec_lat)-6)
  mxlat <- ceiling(max(dhl02$dec_lat)+6)
  mnlon <- floor(min(dhl02$dec_lon)-6)
  mxlon <- ceiling(max(dhl02$dec_lon)+6)
  # make a factor for enlarging the pies
  ftenlpies <- 4.8
  # set a raius for the pies, where the area reflects the number 
  dhl02$rafc <- (sqrt(dhl02$rws)/pi)*ftenlpies
  #
  ((dhl02$rafc/ftenlpies)*pi)^2
  #
  dhl02$rafc/max(dhl02$rafc)*ftenlpies
  #https://guangchuangyu.github.io/2016/12/scatterpie-for-plotting-pies-on-ggplot/
  #world3 <- ggplot2::map_data('world3')
  jitlvl <- 0.017
  
  # also see : https://github.com/tidyverse/ggplot2/issues/2037
  g_plt <- 
    ggplot(data = world3) +
    geom_sf(color = "black", fill = "azure3") +
    scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                    r = rafc),
                                data = dhl02, 
                                cols = colnames(dhl02[,c(2:enc)]),
                                linewidth=0.2) +
    
    scale_color_manual(values=c(rep("black",
                                    length(unique(dhl02[,c(2:enc)]))))) +
    scale_fill_manual(values=alpha(
      c(cl03),
      c(0.7)
    )) +
    # set a blank theme
    theme_bw() +
    # see : https://www.r-bloggers.com/2020/12/visualizing-geospatial-data-in-r-part-2-making-maps-with-ggplot2/
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        color = "grey23", linetype = "dashed", linewidth = 0.1
      ),
      panel.background = ggplot2::element_rect(fill = "grey99")
    ) +
    
    theme(aspect.ratio=1/2.6) +
    geom_scatterpie_legend(max(dhl02$rafc),
                           labeller = function(ra) {c(max(dhl02$rws))},
                           x=-40, y=20 , n=2) +
    
    geom_scatterpie_legend(max(dhl02$rafc)/2,
                           labeller = function(ra) {c(max(dhl02$rws)/2)},
                           x=-40, y=20 , n=2) +
    # set alpha values for color intensity of fill color in point
    #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
    #define limits of the plot
    ggplot2::coord_sf(xlim = c(mnlon, mxlon),
                      ylim = c(mnlat, mxlat),
                      expand = FALSE)
  # change label for legend - Notice that you need to change for all 3 variables
  # you called 'aes' in 'geom_jitter'
  g_plt <- g_plt + labs(fill=lgn.Nm)
  g_plt <- g_plt + labs(color=lgn.Nm)
  g_plt <- g_plt + labs(shape=lgn.Nm)
  g_plt <- g_plt + xlab("Longitude") + ylab("Latitude")
  # #https://www.statology.org/ggplot-background-color/
  g_plt <- g_plt + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                         panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,
  # see the plot
  #g_plt
  lst_plts.G[[ng]] <- g_plt
  # end iteration over alignments
}
# collect plots in to one plot
pA <-   lst_plts.G[[1]] +
  lst_plts.G[[2]] +
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  # see : https://stackoverflow.com/questions/63434957/how-to-adjust-the-font-style-of-tags-with-plot-annotation-in-figures-assembled-w
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold', size =18)) +
  plot_layout(guides = "collect") #+
# evaluate whether to store the plot
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig11_v01_map_genotypes_all.png")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,
         width=210,
         height=297*0.6,
         units="mm",dpi=300)
}
#


#_______________________________________________________________________________
# section 15 -  end - Make haplotype circles on map second attempt
# on all samples
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 16 - start - Make haplotype circles on map second attempt 
# on DK- Germany samples only
#_______________________________________________________________________________
# #install.packages("scatterpie")
# if(!require(scatterpie)){
#   install.packages("scatterpie")
# }

library(scatterpie)

#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)

lst_plts.G <-list()
#iterate over the  sequence of numbers that  reflect the genes 
for (ng in n.f.alvrs.M)
{
  print(ng)
  #}
  subfglts <- letters[n.f.alvrs.M]
  sbfigl <- subfglts[ng]
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm<- gene.nms[ng]
  # Now the alignment has been picked from the list,
  # Then prepare the NJ trees, and the label categories
  # and a color range to assign the labels in the trees
  # also get the labels on the sequences
  LalvrsM <- names(alvrsM)
  
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  # make a name for the legend in the ggplot
  lgn.Nm <- paste0("genotype ",GnNm)
  
  
  #::::
  
  # strsplit names and turn into characters
  # and rowbind the nested lists to a dataframe
  df_ppNm <- data.frame(do.call
                        ('rbind', 
                          strsplit(as.character(names(alvrsM)),
                                   "_")))
  # modify columns names
  colnames(df_ppNm) <- c("smplNm","loc","yer","mtn","GenVar")
  # replace the Mnelei abbr numbers with NCBI accession numbers obtained
  # after depositing sequences on NCBI
  # then get the names
  lbN <- df_ppNm$smplNm
  # get unique location names
  uloc <- unique(df_ppNm$loc)
  # order the unique location names
  uloc <- uloc[order(uloc)]
  # Get the names from the phydat object
  Nmsloc <- names(alvrsM)
  # https://stackoverflow.com/questions/45357806/dplyr-select-and-starts-with-on-multiple-values-in-a-variable-list?noredirect=1&lq=1
  NmslocNEA <- Nmsloc[grepl(paste(NEAloc, collapse="|"),Nmsloc)]
  
  # now subset the bioseq object by the full names that match the NEA samples
  NmslocNEA.tgrp <- paste(NmslocNEA,collapse = "|")
  bsq_pip5 <- alvrsM[grepl(NmslocNEA.tgrp,names(alvrsM))]
  # limit to only include samples collected in 2017 and 2018
  Nmswy17_18 <- NmslocNEA[grepl(paste(c("2017","2018"), collapse="|"),NmslocNEA)]
  # collapse a string of sequences names to grep for
  Nmswy17_18.tgrp <- paste(Nmswy17_18,collapse = "|")
  bsq_pip5 <- bsq_pip5[grepl(Nmswy17_18.tgrp,names(bsq_pip5))]
  # get the names for the sequences
  Lbsq5 <- names(bsq_pip5)
  #make the bioseq object a DNAbin object
  dnb_alsM <-  bioseq::as_DNAbin(bsq_pip5)
  #also assign the names to the sequences
  names(dnb_alsM) <- Lbsq5
  pip <- dnb_alsM
  #::::
  #create haplotypes from dna.bin
  pipHaps <- pegas::haplotype(pip)
  #prepare hpt table
  ih<-with(
    utils::stack(setNames(attr(pipHaps, "index"), rownames(pipHaps))),
    table(hap=ind, pop=names(pip)[values]))
  # make it a data frame
  dih01 <- as.data.frame(ih)
  #limit to include only 'Freq' that equals 1
  dih02 <- dih01[dih01$Freq == 1,]	
  #get pop groups from row names
  # split the row name string by a character 
  dihpg <- strsplit(as.character(dih02$pop), "_")
  # get the 2 nd element per vector in the list - this holds the abbreviated 
  # location name, the 3rd element is the sampling year, the 5th is the 
  # gene variant
  poploc <- sapply(dihpg, "[[", 2)
  popyear <- sapply(dihpg, "[[", 3)
  popGvar <- sapply(dihpg, "[[", 5)
  #make it a table
  hsll3 <- table(dih02$hap, poploc)
  # make the roman numerals arabic numerals
  rwNm.hsll3 <- as.numeric(as.roman(rownames(hsll3)))
  rwNm.hsll3 <- paste0("G",rwNm.hsll3)
  # get the names
  rownames(hsll3) <- rwNm.hsll3
  # transpose
  thl1.3 <- t(hsll3)
  # make it a data frame
  dhl01 <- as.data.frame(thl1.3)
  # reshape the data frame for long to wide
  dhl02 <- reshape(data=dhl01,idvar="poploc",
                   v.names = "Freq",
                   timevar = "Var2",
                   direction="wide")
  enc <- ncol(dhl02)
  # substitute in the column names
  clNdh02 <- colnames(dhl02)
  clNdh02 <- gsub("Freq\\.","",clNdh02)
  colnames(dhl02) <- clNdh02
  #make a viridis colour range
  cl03 <- pals::viridis(length(unique(dhl02[,c(2:enc)])))
  cl03 <- pals::inferno(length(unique(dhl02[,c(2:enc)])))
  ncolHptloc03 <- length(unique(colnames(dhl02)[2:enc]))
  cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                   "yellow","white")
  cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
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
  
  # scandiWorld <- rnaturalearth::ne_countries(country=c(
  #   'Denmark',
  #   'Sweden',
  #   'Norway',
  #   'Germany'
  #   ),scale = 10, returnclass = "sf")
  # sum up rows , avoid confusion with function 'rowsum'
  dhl02$rws <- rowSums( dhl02[,c(2:enc)])
  # Get latitude and longitude positions
  dhl02$dec_lat <- df_clo04$dec_lat2[match(dhl02$poploc, df_clo04$loc4Lett)]
  dhl02$dec_lon <- df_clo04$dec_lon2[match(dhl02$poploc, df_clo04$loc4Lett)]
  dhl02$dec_lat <- as.numeric(dhl02$dec_lat)
  dhl02$dec_lon <- as.numeric(dhl02$dec_lon)
  # get limits for the latitudes and longitudes to use in the ggplot
  mnlat <- floor(min(dhl02$dec_lat)-0)
  mxlat <- ceiling(max(dhl02$dec_lat)+1)
  mnlon <- floor(min(dhl02$dec_lon)-1)
  mxlon <- ceiling(max(dhl02$dec_lon)+2)
  # make a factor for enlarging the pies
  ftenlpies <- 0.6
  # set a raius for the pies, where the area reflects the number 
  dhl02$rafc <- (sqrt(dhl02$rws)/pi)*ftenlpies
  #
  ((dhl02$rafc/ftenlpies)*pi)^2
  #
  dhl02$rafc/max(dhl02$rafc)*ftenlpies
  #https://guangchuangyu.github.io/2016/12/scatterpie-for-plotting-pies-on-ggplot/
  #world3 <- ggplot2::map_data('world3')
  jitlvl <- 0.017
  
  # also see : https://github.com/tidyverse/ggplot2/issues/2037
  g_plt <- 
    ggplot(data = world3) +
    #ggplot(data = scandiWorld) +
    
    geom_sf(color = "black", fill = "azure3") +
    scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                    r = rafc),
                                data = dhl02, 
                                cols = colnames(dhl02[,c(2:enc)]),
                                linewidth=0.2) +
    
    scale_color_manual(values=c(rep("black",
                                    length(unique(dhl02[,c(2:enc)]))))) +
    scale_fill_manual(values=alpha(
      c(cl03),
      c(0.7)
    )) +
    # set a blank theme
    theme_bw() +
    # see : https://www.r-bloggers.com/2020/12/visualizing-geospatial-data-in-r-part-2-making-maps-with-ggplot2/
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        color = "grey23", linetype = "dashed", linewidth = 0.1
      ),
      panel.background = ggplot2::element_rect(fill = "grey99")
    ) +
    
    theme(aspect.ratio=1/1.8) +
    #coord_equal() +
    geom_scatterpie_legend(max(dhl02$rafc),
                           labeller = function(ra) {c(max(dhl02$rws))},
                           x=7, y=57 , n=2) +
    
    geom_scatterpie_legend(max(dhl02$rafc)/2,
                           labeller = function(ra) {c(max(dhl02$rws)/2)},
                           x=7, y=57 , n=2) +
    # set alpha values for color intensity of fill color in point
    #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
    #define limits of the plot
    ggplot2::coord_sf(xlim = c(mnlon, mxlon),
                      ylim = c(mnlat, mxlat),
                      expand = FALSE)
  # change label for legend - Notice that you need to change for all 3 variables
  # you called 'aes' in 'geom_jitter'
  g_plt <- g_plt + labs(fill=lgn.Nm)
  g_plt <- g_plt + labs(color=lgn.Nm)
  g_plt <- g_plt + labs(shape=lgn.Nm)
  g_plt <- g_plt + xlab("Longitude") + ylab("Latitude")
  # #https://www.statology.org/ggplot-background-color/
  g_plt <- g_plt + theme(panel.background = element_rect(fill = 'white', color = 'black'),
                         panel.grid.major = element_line(color = 'azure3')) #, linetype = 'dotted'))#,
  # see the plot
  #g_plt
  lst_plts.G[[ng]] <- g_plt
  # end iteration over alignments
}
# collect plots in to one plot
pA <-   lst_plts.G[[1]] +
  lst_plts.G[[2]] +
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  # see : https://stackoverflow.com/questions/63434957/how-to-adjust-the-font-style-of-tags-with-plot-annotation-in-figures-assembled-w
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold', size =18)) +
  plot_layout(guides = "collect") #+
# evaluate whether to store the plot
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig11_v02_map_genotypes_DK_Germ.png")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,
         width=210*0.9,
         height=297*0.6,
         units="mm",dpi=300)
}
#


#_______________________________________________________________________________
# section 16 -  end - Make haplotype circles on map second attempt 
# on DK- Germany samples only
#_______________________________________________________________________________


#_______________________________________________________________________________
# section 17 -  start - Make NMDS plots - all samples
#_______________________________________________________________________________


#_______________________________________________________________________________
# Make MDS plot
#_______________________________________________________________________________
#https://www.statmethods.net/advstats/mds.html

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)

lst.plt12b <- list()
lst.plt12c <- list()
lst_plts.G <- list()
#iterate over the  sequence of numbers that  reflect the genes 
for (ng in n.f.alvrs.M)
{
  print(ng)
  #}
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
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  # make a name for the legend in the ggplot
  lgn.Nm <- paste0("genotype ",GnNm)
  # strsplit names and turn into characters
  # and rowbind the nested lists to a dataframe
  df_ppNm <- data.frame(do.call
                        ('rbind', 
                          strsplit(as.character(names(alvrsM)),
                                   "_")))
  # modify columns names
  colnames(df_ppNm) <- c("smplNm","loc","yer","mtn","GenVar")
  # replace the Mnelei abbr numbers with NCBI accession numbers obtained
  # after depositing sequences on NCBI
  # then get the names
  lbN <- df_ppNm$smplNm
  # get unique location names
  uloc <- unique(df_ppNm$loc)
  # order the unique location names
  uloc <- uloc[order(uloc)]
  pip <- alvrsM
  #create haplotypes from dna.bin
  pipHaps <- pegas::haplotype(pip)
  #prepare hpt table
  ih<-with(
    utils::stack(setNames(attr(pipHaps, "index"), rownames(pipHaps))),
    table(hap=ind, pop=names(pip)[values]))
  # make it a data frame
  dih01 <- as.data.frame(ih)
  #limit to include only 'Freq' that equals 1
  dih02 <- dih01[dih01$Freq == 1,]	
  #get pop groups from row names
  # split the row name string by a character 
  dihpg <- strsplit(as.character(dih02$pop), "_")
  # get the 2 nd element per vector in the list - this holds the abbreviated 
  # location name, the 3rd element is the sampling year, the 5th is the 
  # gene variant
  poploc <- sapply(dihpg, "[[", 2)
  popyear <- sapply(dihpg, "[[", 3)
  popGvar <- sapply(dihpg, "[[", 5)
  #make it a table
  hsll3 <- table(dih02$hap, poploc)
  # make the roman numerals arabic numerals
  rwNm.hsll3 <- as.numeric(as.roman(rownames(hsll3)))
  rwNm.hsll3 <- paste0("G",rwNm.hsll3)
  # get the names
  rownames(hsll3) <- rwNm.hsll3
  # transpose
  thl1.3 <- t(hsll3)
  # make it a data frame
  dhl01 <- as.data.frame(thl1.3)
  # reshape the data frame for long to wide
  dhl02 <- reshape(data=dhl01,idvar="poploc",
                   v.names = "Freq",
                   timevar = "Var2",
                   direction="wide")
  enc <- ncol(dhl02)
  # substitute in the column names
  clNdh02 <- colnames(dhl02)
  clNdh02 <- gsub("Freq\\.","",clNdh02)
  colnames(dhl02) <- clNdh02
  #count number of columns in data frame
  nc02 <-ncol(dhl02)
  #get columns from the second the 7th last column 
  dhl03 <- dhl02[,2:(nc02)]
  # assign row names to the data frame
  rownames(dhl03) <- dhl02[,1]
  d <- dist(dhl03) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  # plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
       main="Metric MDS", type="n")
  text(x, y, labels = row.names(dhl03), cex=.7) 
  # Nonmetric MDS
  # N rows (objects) x p columns (variables)
  # each row identified by a unique row name
  library(MASS)
  library(vegan)
  d <- dist(dhl03) # euclidean distances between the rows
  #fit <- MASS::isoMDS(d, k=2) # k is the number of dim
  mmdsfit  <- vegan::monoMDS(d, k=2)
  #fit # view results
  mmdsfit$points[,1]
  # end plot
  dev.off()
  # define output file name for plot
  figname01 <- paste0("Fig12a_v",pnng,"_",GnNm,"_NMDS_plot_all_smpls.png")
  pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
  # set to save plot as png 
  png(pthfignm01,
      width = 210, 
      height = 297,
      units="mm",
      res =300)
  # try this example code
  # From this website: https://stackoverflow.com/questions/12302366/positioning-axes-labels/12302557#12302557
  plot(1:100, cumsum(rnorm(100)), type="l", mgp=c(2.4,0.2,.5), las=1)
  #add extra space to the right of the plot
  par(mar=c(4, 4, 2, 2), xpd=FALSE)
  par(oma=c(0, 0, 0, 0))
  #reset this parameter
  par(mfrow = c(1, 1)) 
  # begin plot, with defined borders
  x <- mmdsfit$points[,1]
  y <- mmdsfit$points[,2]
  plot(x, y, xlab="NMDS1", ylab="NMDS2",
       type="n",
       mgp=c(2.4,1.1,0.00001), las=1
  )
  text(x, y, labels = row.names(dhl03), cex=1.8) 
  # end plot
  dev.off()
  #reset this parameter
  par(mfrow = c(1, 1)) 
  # get the x and y values
  NMDS1 <- x
  NMDS2 <- y
  dhl03.lbls <- row.names(dhl03)
  df_NMDS_01 <- as.data.frame(cbind(NMDS1,NMDS2,dhl03.lbls))
  ggplot(df_NMDS_01) +
    geom_point(aes(x=NMDS1,
                   y=NMDS2))
  #https://www.r-statistics.com/2016/01/multidimensional-scaling-with-r-from-mastering-data-analysis-with-r/
  #make genind and hierfstat objects by region
  gi_pip_local <- adegenet::DNAbin2genind(pip,pop=poploc)
  #make hierfstat objects - Note that it is required to set 'pop' to the 
  # vector that matches the population to which the indivudal belongs
  hfst_pip_lo <- hierfstat::genind2hierfstat(gi_pip_local,pop=poploc)
  hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
  #replace all NAs with 0
  hfst_pip_lo <- hfst_pip_lo %>% replace(is.na(.), 0)
  # Remove rows where the pop name is unique - 
  # as the 'hierfstat::pairwise.WCfst'
  # function cannot handle the unique pop rows only  represented by only a single individual
  hfst_pip_lo <- hfst_pip_lo[hfst_pip_lo$pop %in% hfst_pip_lo$pop[duplicated(hfst_pip_lo$pop)],]
  hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
  #pw_Nei_pip_un <- hierfstat::pairwise.neifst(hfst_pip_un, diploid = FALSE)
  pw_Nei_pip_lo <- hierfstat::pairwise.neifst(hfst_pip_lo, diploid = FALSE)
  # make the data frame as a heat map
  df_pw_pip_lo <- as.data.frame(pw_Nei_pip_lo)
  #replace NAs with zeroes
  df_pw_pip_lo[is.na(df_pw_pip_lo)] = 0
  # assign it a different name 
  dlo <- df_pw_pip_lo
  #copy the dna.bin object
  pip3 <- pip
  #make the DNAbin object a genind object
  geni_pip3 <- adegenet::DNAbin2genind(pip3)
  #make the genind object a dataframe
  df_pip3 <- adegenet::genind2df(geni_pip3)
  #get the row names
  orig_rwnm <- row.names(df_pip3)
  # replace all NAs
  df_pip3 <- df_pip3 %>% replace(is.na(.), "-")
  #make the date frame a matrix and a DNAbin object again
  dnb_pip3 <- as.DNAbin(as.matrix(df_pip3))
  # get abbreviated location name
  ablo <- poploc
  library(MASS)
  d <- dist(dlo) # euclidean distances between the rows
  fit <- MASS::isoMDS(d, k=2) # k is the number of dim
  # make the fitted NMDS a data frame
  fmds <- as.data.frame(fit)
  # find the stress level and the moin and the max
  minstr <- min(fmds$stress)
  maxstr <- max(fmds$stress)
  strslvl <- ifelse(minstr==maxstr,maxstr, paste0(minstr,"-",maxstr) )
  strslvl <- round(strslvl,2)
  # make data frame with overall location names for the 4 letter location codes
  lscd <- c("BalS", "CasS", "CWAt", "FBog", "FKer", "GBus", "GHel", "GKie", 
            "GWis", "JLim", "JMar", "Medi", "NEAt", "Neth", "NWAt", "SBal", 
            "SSko", "USAG", "USAP", "USAW")
  lncd <- c("NE Atlantic", "Caspian Sea", "Central W Atlantic",
            "NE Atlantic", "NE Atlantic", "NE Atlantic", "NE Atlantic", "NE Atlantic", 
            "NE Atlantic", "NE Atlantic", "NE Atlantic", "Mediterranean", 
            "NE Atlantic", "NE Atlantic",
            "NW Atlantic", "NE Atlantic", 
            "NE Atlantic", "NW Atlantic", "NW Atlantic", "NW Atlantic")
  # and assemble in a data frame
  df_lnNm <- as.data.frame(cbind(lscd,lncd))
  # match to get overall location
  fmds$ov.loc <- df_lnNm$lncd[match(rownames(fmds),df_lnNm$lscd)]
  # get unique over locations
  locs3 <- unique(fmds$ov.loc)
  # count unique overall locations
  nlocs3 <- length(locs3)
  # make series of numbers to use for pch shapes
  pchnmbs <- rep(c(21:24),2)
  pchnmbs <- pchnmbs[1:nlocs3]
  # make series of colors to use for coloring fill
  colnmbs <- rep(c("white","black","black"),4)
  colnmbs <- colnmbs[1:nlocs3]
  #bind together in a data frame
  df_pcl04 <- as.data.frame(cbind(locs3,pchnmbs,colnmbs))
  fmds$pchnmb <- df_pcl04$pchnmbs[match(fmds$ov.loc, df_pcl04$locs3)]
  fmds$flcnmb <- df_pcl04$colnmbs[match(fmds$ov.loc, df_pcl04$locs3)]
  # load the library
  library(ggplot2)
  library(ggrepel)
  nrwsfmds <- length(rownames(fmds))
  # plot it with ggplot
  p14 <- ggplot(fmds, aes(points.1, -points.2, label = rownames(fmds))) +
    geom_point(aes(shape= ov.loc, 
                   fill=ov.loc),size=3.0) +  
    geom_text(check_overlap = TRUE,
              vjust = -0.75) #+ theme_minimal() + xlab('') + ylab('') +
  # add text to the plot to indicate stress level
  p14 <- p14 +                               # Add text element to plot
    annotate("text", x = 0, y = -5.5, label = paste0("stress level: ",strslvl))
  p14 <- p14 +  scale_shape_manual(values = c(as.integer(pchnmbs))) +
    scale_fill_manual(values = c(colnmbs))
  p14 <- p14 + theme(panel.background = element_rect(fill = 'white', 
                                                     color = 'white'))#,
  # change label for legend - Notice that you need to change for all 3 variables
  # you called 'aes' in 'geom_jitter'
  p14 <- p14 + labs(fill='Location')
  p14 <- p14 + labs(color='Location')
  p14 <- p14 + labs(shape='Location')
  p14 <- p14 + xlab("NMDS1") + ylab("NMDS2")
  # add border around plot
  p14 <- p14 + theme(panel.border = element_rect(
    color = "black",
    fill = NA,
    size = 1.0))
  # change background of legend
  p14 <- p14 + theme(legend.key = element_rect(fill = "white"))
  # https://stackoverflow.com/questions/66102618/ggplot-increasing-the-distance-between-axis-labels-and-axis-ticks
  p14 <- p14 + theme(axis.text.x =
                       element_text(color="black", 
                                    size=12.4,
                                    margin = margin(t = 2, 
                                                    unit = "mm")),
                     axis.text.y = 
                       element_text(color="black", 
                                    size=12.4,
                                    margin = margin(l = 0, r=2,
                                                    unit = "mm")))
  # collect plot in a list
  lst.plt12b[[ng]] <- p14
  #make filename to save plot to
  figname14 <- paste0("Fig12b_v",pnng,"_smpl_location_NMDS_",GnNm,".png")
  figname02 <- paste(wd00_wd05,"/",figname14,sep="")
  if(bSaveFigures==T){
    ggsave(p14,file=figname02,
           #width=210,height=297,
           width=210,height=(297*0.5),
           units="mm",dpi=300)
  }
  #_______________________________________________________________________________
  #make genind and hierfstat objects by year
  gi_pip_local <- adegenet::DNAbin2genind(pip,pop=popyear)
  #make hierfstat objects - Note that it is required to set 'pop' to the 
  # vector that matches the population to which the indivudal belongs
  hfst_pip_lo <- hierfstat::genind2hierfstat(gi_pip_local,pop=popyear)
  hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
  #replace all NAs with 0
  hfst_pip_lo <- hfst_pip_lo %>% replace(is.na(.), 0)
  # Remove rows where the pop name is unique - 
  # as the 'hierfstat::pairwise.WCfst'
  # function cannot handle the unique pop rows only  represented by only a single individual
  hfst_pip_lo <- hfst_pip_lo[hfst_pip_lo$pop %in% hfst_pip_lo$pop[duplicated(hfst_pip_lo$pop)],]
  hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
  #pw_Nei_pip_un <- hierfstat::pairwise.neifst(hfst_pip_un, diploid = FALSE)
  pw_Nei_pip_lo <- hierfstat::pairwise.neifst(hfst_pip_lo, diploid = FALSE)
  # make the data frame as a heat map
  df_pw_pip_lo <- as.data.frame(pw_Nei_pip_lo)
  #replace NAs with zeroes
  df_pw_pip_lo[is.na(df_pw_pip_lo)] = 0
  # assign it a different name 
  dlo <- df_pw_pip_lo
  # copy the object
  dly <- dlo
  # exclude the years that are too variable
  yrs_to_excl <- as.character(paste(c("unknown",2007),collapse = "|"))
  indx.r.tk <- which(!grepl(yrs_to_excl,row.names(dly)))
  indx.c.tk <- which(!grepl(yrs_to_excl,colnames(dly)))
  dly <- dly[indx.r.tk,indx.c.tk]
  library(MASS)
  d <- dist(dly) # euclidean distances between the rows
  # fit a NMDS
  fit <- MASS::isoMDS(d, k=2, p=80, maxit=100, tol=1e-3) # k is the number of dim
  # make the fitted NMDS a data frame
  fmds <- as.data.frame(fit)
  # get the stress level
  minstr <- min(fmds$stress)
  maxstr <- max(fmds$stress)
  strslvl <- ifelse(minstr==maxstr,maxstr, paste0(minstr,"-",maxstr) )
  strslvl <- round(strslvl,2)
  # add row names as a column
  fmds$smplyear <- row.names(fmds)
  # get the number of unique sampling years
  nyers3 <- length(unique(fmds$smplyear))
  smpl_yrs <- unique(fmds$smplyear)
  # make series of numbers to use for pch shapes
  pchnmbs <- rep(c(21:25),4)
  pchnmbs <- pchnmbs[1:nyers3]
  # make series of colors to use for coloring fill
  colnmbs <- rep(c("white","white","white","white","white"
                   ,"black","black","black","black","black"
                   ,"grey","grey","grey"),6)
  colnmbs <- colnmbs[1:nyers3]
  #bind together in a data frame
  df_pcl04 <- as.data.frame(cbind(smpl_yrs,pchnmbs,colnmbs))
  #check if there are an equal number of colored symbols
  nyers3==length(unique(paste0(df_pcl04$pchnmbs,df_pcl04$colnmbs)))
  # match to fitted MDS data frame
  fmds$pchnmb <- df_pcl04$pchnmbs[match(fmds$smplyear, df_pcl04$smpl_yrs)]
  fmds$flcnmb <- df_pcl04$colnmbs[match(fmds$smplyear, df_pcl04$smpl_yrs)]
  # get the number of rows
  nrwsfmds <- length(rownames(fmds))
  # plot it with ggplot
  p15 <- ggplot(fmds, aes(points.1, -points.2, label = rownames(fmds))) +
    geom_point(aes(shape= smplyear, 
                   fill=smplyear),size=3.0) +
    ggrepel::geom_text_repel()
  # add text to the plot to indicate stress level
  p15 <- p15 +                               # Add text element to plot
    annotate("text", x = 0, y = -5.5, label = paste0("stress level: ",strslvl))
  p15 <- p15 +  scale_shape_manual(values = c(as.integer(pchnmbs))) +
    scale_fill_manual(values = c(colnmbs))
  p15 <- p15 + theme(panel.background = element_rect(fill = 'white', color = 'white'))#,
  p15 <- p15 + labs(fill='Year')
  p15 <- p15 + labs(color='Year')
  p15 <- p15 + labs(shape='Year')
  p15 <- p15 + xlab("NMDS1") + ylab("NMDS2")
  # add border around plot
  p15 <- p15 + theme(panel.border = element_rect(color = "black",
                                                 fill = NA,
                                                 size = 1.0))
  # change background of legend
  p15 <- p15 + theme(legend.key = element_rect(fill = "white"))
  # https://stackoverflow.com/questions/66102618/ggplot-increasing-the-distance-between-axis-labels-and-axis-ticks
  p15 <- p15 + theme(axis.text.x = 
                       element_text(color="black", 
                                    size=12.4,
                                    margin = margin(t = 2, 
                                                    unit = "mm")),
                     
                     axis.text.y = 
                       element_text(color="black", 
                                    size=12.4,
                                    margin = margin(l = 0, r=2,
                                                    unit = "mm")))
  # collect plot in a list
  lst.plt12c[[ng]] <- p15
  #make filename to save plot to
  figname15 <- paste0("Fig12c_v",pnng,"_smpl_year_NMDS_",GnNm,".png")
  # define pathe and file together in a string
  figname02 <- paste(wd00_wd05,"/",figname15,sep="")
  # make an if test to check whether the plot should be saved
  # i.e. set to FALSE if you do not want the plot to be saved
  if(bSaveFigures==T){
    ggsave(p15,file=figname02,
           #width=210,height=297,
           width=210,
           height=(297*0.5),
           units="mm",dpi=300)
  }
  #_______________________________________________________________________________
  pip3 <- pip
  dnb_pip3 <- pip3
  #make the date frame a matrix and a DNAbin object again
  dnb_pip3 <- as.DNAbin(as.matrix(df_pip3))
  #
  library(MASS)
  # set pair wise deletion to TRUE to remove identical sequences
  d <- ape::dist.dna(dnb_pip3, model = "raw", pairwise.deletion = TRUE)
  #replace NAs with zero
  d[is.na(d)] <- 0
  #fit <- isoMDS(d, k=2) # k is the number of dim
  fmMDS<- metaMDS(d, k=2)
  # get the label categories for the fitted MDS
  lcNm0 <- labels(dnb_pip3)
  
  lcNm0 <- strsplit(lcNm0,"_")
  lcNm1 <- sapply(lcNm0, "[[", 1)
  lcNm2 <- sapply(lcNm0, "[[", 2)
  lcNm3 <- sapply(lcNm0, "[[", 3)
  # get the short location name
  shNm3 <-  lcNm2
  #fmMDS$species <- lcNm4
  fmMDS$species <- shNm3
  #row.names(fmMDS$points) <- lcNm4
  row.names(fmMDS$points) <- shNm3
  strslvl2 <- fmMDS$stress
  strslvl2 <- round(strslvl2,3)
  # define output file name for plot
  figname01 <- paste0("Fig12d_v",pnng,"_smpl_location_NMDS_plot",GnNm,".png")
  pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
  # set to save plot as png file 
  png(pthfignm01,
      width = 210, 
      height = 297,
      units="mm",
      res =300)
  #add extra space to the right of the plot
  par(mar=c(5, 5, 2, 2), xpd=FALSE)
  par(oma=c(0, 0, 0, 0))
  #reset this parameter
  par(mfrow = c(1, 1)) 
  #dev.off()
  ordiplot(fmMDS,type="n",cex.axis=1.8,cex.lab = 1.8)
  text(x = 0.2, y = -0.4, cex=1.8,                # Add text element
       paste0("stress level: ",strslvl2))
  orditorp(fmMDS,display="sites",
           air=0.2,cex=1.25)
  # end plot
  dev.off()
  #reset this parameter
  par(mfrow = c(1, 1)) 
  #_______________________________________________________________________________
  library(MASS)
  # set pair wise deletion to TRUE to remove identical sequences
  d <- ape::dist.dna(dnb_pip3, model = "raw", pairwise.deletion = TRUE)
  #replace NAs with zero
  d[is.na(d)] <- 0
  #fit <- isoMDS(d, k=2) # k is the number of dim
  fmMDS<- metaMDS(d, k=2)
  # make the fitted NMDS a data frame
  lcNm0 <- labels(dnb_pip3)
  lcNm0 <- strsplit(lcNm0,"_")
  lcNm1 <- sapply(lcNm0, "[[", 1)
  lcNm2 <- sapply(lcNm0, "[[", 2)
  lcNm3 <- sapply(lcNm0, "[[", 3)
  #combine to a data frame
  df_smp4 <- as.data.frame(cbind(lcNm1,lcNm3))
  #match to get sampling year
  df_smp4$smply <- lcNm3
  #
  m.pnt.NMDS1 <- mean(fmMDS$points[,1])
  m.pnt.NMDS2 <- min(fmMDS$points[,2])+((max(fmMDS$points[,2]) - min(fmMDS$points[,2]))*(1/8))
  # # get sampling year and unqiue sample number from sequence name
  fmMDS$species <- df_smp4$smply
  lcNm5 <- df_smp4$smply
  row.names(fmMDS$points) <- lcNm5
  strslvl2 <- fmMDS$stress
  strslvl2 <- round(strslvl2,3)
  # define output file name for plot
  figname01 <- paste0("Fig12e_v",pnng,"_",GnNm,"_smpl_year_NMDS_plot.png")
  pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
  # set to save plot as png 
  png(pthfignm01,
      width = 210, 
      height = 297,
      units="mm",
      res =300)
  #add extra space to the right of the plot
  par(mar=c(5, 5, 2, 2), xpd=FALSE)
  par(oma=c(0, 0, 0, 0))
  #reset this parameter
  par(mfrow = c(1, 1)) 
  ordiplot(fmMDS,type="n",cex.axis=1.8,cex.lab = 1.8)
  text(x = m.pnt.NMDS1, y = m.pnt.NMDS2, cex=0.8,                # Add text element
       paste0("stress level: ",strslvl2))
  orditorp(fmMDS,display="sites",
           air=0.2,cex=1.25)
  # end plot
  dev.off()
  #reset this parameter
  par(mfrow = c(1, 1)) 
  #_______________________________________________________________________________
  # end iteration over alignments
}


library(patchwork)
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p01t <- lst.plt12b[[1]] #+ labs(title = "a", face="bold")#,
p02t <- lst.plt12b[[2]] #+ labs(title = "b", face="bold")#,
#p03t <- p03 +   ggtitle('Plot 4', size =16) # (title = "c", face="bold")#,
#p03t <- p03 + theme(plot.title = element_text(size = 12, face = "bold"))
#p03t <- p03 + theme(plot.title = element_text(face="bold", size=18))

pA <-   p01t + 
  p02t +
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  # see : https://stackoverflow.com/questions/63434957/how-to-adjust-the-font-style-of-tags-with-plot-annotation-in-figures-assembled-w
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold', size =22)) +
  plot_layout(guides = "collect") #+
#plot_annotation(caption=pthinf01) #& theme(legend.position = "bottom")
#p
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig12b_v05_all_ITS1_and_ITS2_NMDS_all_smpl_regions.png")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,
         width=210,
         height=297,
         units="mm",dpi=300)
}


library(patchwork)
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p01t <- lst.plt12c[[1]] #+ labs(title = "a", face="bold")#,
p02t <- lst.plt12c[[2]] #+ labs(title = "b", face="bold")#,
#p03t <- p03 +   ggtitle('Plot 4', size =16) # (title = "c", face="bold")#,
#p03t <- p03 + theme(plot.title = element_text(size = 12, face = "bold"))
#p03t <- p03 + theme(plot.title = element_text(face="bold", size=18))

pA <-   p01t + 
  p02t +
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  # see : https://stackoverflow.com/questions/63434957/how-to-adjust-the-font-style-of-tags-with-plot-annotation-in-figures-assembled-w
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold', size =22)) +
  plot_layout(guides = "collect") #+
#plot_annotation(caption=pthinf01) #& theme(legend.position = "bottom")
#p
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig12c_v05_all_ITS1_and_ITS2_NMDS_all_smpl_regions.png")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,
         width=210,
         height=297,
         units="mm",dpi=300)
}
#_______________________
#_______________________________________________________________________________
#_______________________________________________________________________________
# section 17 -  end - Make NMDS plots - all samples
#_______________________________________________________________________________

#_______________________________________________________________________________
# section 18 -  start - Make NMDS plots - only DK and Germany samples
#_______________________________________________________________________________


#combine the alignments in a list
lst_alvrs.M.ITS <- list(alvrs.M.ITS1,alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1","ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1,n.alvrs.M,1)

lst_plts.G <-list()
lst.plt13b <- list()
#iterate over the  sequence of numbers that  reflect the genes 
for (ng in n.f.alvrs.M)
{
  print(ng)
  #}
  subfglts <- letters[n.f.alvrs.M]
  sbfigl <- subfglts[ng]
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm<- gene.nms[ng]
  # Now the alignment has been picked from the list,
  # Then prepare the NJ trees, and the label categories
  # and a color range to assign the labels in the trees
  # also get the labels on the sequences
  LalvrsM <- names(alvrsM)
  
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng)<2,stringr::str_pad(ng, 2, pad = "0"),ng)
  # make a name for the legend in the ggplot
  lgn.Nm <- paste0("genotype ",GnNm)
  
  
  #::::
  
  # strsplit names and turn into characters
  # and rowbind the nested lists to a dataframe
  df_ppNm <- data.frame(do.call
                        ('rbind', 
                          strsplit(as.character(names(alvrsM)),
                                   "_")))
  # modify columns names
  colnames(df_ppNm) <- c("smplNm","loc","yer","mtn","GenVar")
  # replace the Mnelei abbr numbers with NCBI accession numbers obtained
  # after depositing sequences on NCBI
  # then get the names
  lbN <- df_ppNm$smplNm
  # get unique location names
  uloc <- unique(df_ppNm$loc)
  # order the unique location names
  uloc <- uloc[order(uloc)]
  # Get the names from the phydat object
  Nmsloc <- names(alvrsM)
  # https://stackoverflow.com/questions/45357806/dplyr-select-and-starts-with-on-multiple-values-in-a-variable-list?noredirect=1&lq=1
  NmslocNEA <- Nmsloc[grepl(paste(NEAloc, collapse="|"),Nmsloc)]
  
  # now subset the bioseq object by the full names that match the NEA samples
  NmslocNEA.tgrp <- paste(NmslocNEA,collapse = "|")
  bsq_pip5 <- alvrsM[grepl(NmslocNEA.tgrp,names(alvrsM))]
  # limit to only include samples collected in 2017 and 2018
  Nmswy17_18 <- NmslocNEA[grepl(paste(c("2017","2018"), collapse="|"),NmslocNEA)]
  # collapse a string of sequences names to grep for
  Nmswy17_18.tgrp <- paste(Nmswy17_18,collapse = "|")
  bsq_pip5 <- bsq_pip5[grepl(Nmswy17_18.tgrp,names(bsq_pip5))]
  # get the names for the sequences
  Lbsq5 <- names(bsq_pip5)
  #make the bioseq object a DNAbin object
  dnb_alsM <-  bioseq::as_DNAbin(bsq_pip5)
  #also assign the names to the sequences
  names(dnb_alsM) <- Lbsq5
  pip <- dnb_alsM
  
  #create haplotypes from dna.bin
  pipHaps <- pegas::haplotype(pip)
  #prepare hpt table
  ih<-with(
    utils::stack(setNames(attr(pipHaps, "index"), rownames(pipHaps))),
    table(hap=ind, pop=names(pip)[values]))
  # make it a data frame
  dih01 <- as.data.frame(ih)
  #limit to include only 'Freq' that equals 1
  dih02 <- dih01[dih01$Freq == 1,]	
  #get pop groups from row names
  # split the row name string by a character 
  dihpg <- strsplit(as.character(dih02$pop), "_")
  # get the 2 nd element per vector in the list - this holds the abbreviated 
  # location name, the 3rd element is the sampling year, the 5th is the 
  # gene variant
  poploc <- sapply(dihpg, "[[", 2)
  popyear <- sapply(dihpg, "[[", 3)
  popGvar <- sapply(dihpg, "[[", 5)
  #make it a table
  hsll3 <- table(dih02$hap, poploc)
  # make the roman numerals arabic numerals
  rwNm.hsll3 <- as.numeric(as.roman(rownames(hsll3)))
  rwNm.hsll3 <- paste0("G",rwNm.hsll3)
  # get the names
  rownames(hsll3) <- rwNm.hsll3
  # transpose
  thl1.3 <- t(hsll3)
  # make it a data frame
  dhl01 <- as.data.frame(thl1.3)
  # reshape the data frame for long to wide
  dhl02 <- reshape(data=dhl01,idvar="poploc",
                   v.names = "Freq",
                   timevar = "Var2",
                   direction="wide")
  enc <- ncol(dhl02)
  # substitute in the column names
  clNdh02 <- colnames(dhl02)
  clNdh02 <- gsub("Freq\\.","",clNdh02)
  colnames(dhl02) <- clNdh02
  #count number of columns in data frame
  nc02 <-ncol(dhl02)
  #get columns from the second the 7th last column 
  dhl03 <- dhl02[,2:(nc02)]
  # assign row names to the data frame
  rownames(dhl03) <- dhl02[,1]
  d <- dist(dhl03) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  # plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
       main="Metric MDS", type="n")
  text(x, y, labels = row.names(dhl03), cex=.7) 
  # Nonmetric MDS
  # N rows (objects) x p columns (variables)
  # each row identified by a unique row name
  library(MASS)
  library(vegan)
  d <- dist(dhl03) # euclidean distances between the rows
  #fit <- MASS::isoMDS(d, k=2) # k is the number of dim
  mmdsfit  <- vegan::monoMDS(d, k=2)
  #fit # view results
  mmdsfit$points[,1]
  # end plot
  dev.off()
  # define output file name for plot
  figname01 <- paste0("Fig13a_v",pnng,"_",GnNm,"_NMDS_plot_DK_Germ_smpls_DK_Germ_smpls.png")
  pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
  # set to save plot as png 
  png(pthfignm01,
      width = 210, 
      height = 297,
      units="mm",
      res =300)
  # try this example code
  # From this website: https://stackoverflow.com/questions/12302366/positioning-axes-labels/12302557#12302557
  plot(1:100, cumsum(rnorm(100)), type="l", mgp=c(2.4,0.2,.5), las=1)
  #add extra space to the right of the plot
  par(mar=c(4, 4, 2, 2), xpd=FALSE)
  par(oma=c(0, 0, 0, 0))
  #reset this parameter
  par(mfrow = c(1, 1)) 
  # begin plot, with defined borders
  x <- mmdsfit$points[,1]
  y <- mmdsfit$points[,2]
  plot(x, y, xlab="NMDS1", ylab="NMDS2",
       type="n",
       mgp=c(2.4,1.1,0.00001), las=1
  )
  text(x, y, labels = row.names(dhl03), cex=1.8) 
  # end plot
  dev.off()
  #reset this parameter
  par(mfrow = c(1, 1)) 
  # get the x and y values
  NMDS1 <- x
  NMDS2 <- y
  dhl03.lbls <- row.names(dhl03)
  df_NMDS_01 <- as.data.frame(cbind(NMDS1,NMDS2,dhl03.lbls))
  ggplot(df_NMDS_01) +
    geom_point(aes(x=NMDS1,
                   y=NMDS2))
  #https://www.r-statistics.com/2016/01/multidimensional-scaling-with-r-from-mastering-data-analysis-with-r/
  #make genind and hierfstat objects by region
  gi_pip_local <- adegenet::DNAbin2genind(pip,pop=poploc)
  #make hierfstat objects - Note that it is required to set 'pop' to the 
  # vector that matches the population to which the indivudal belongs
  hfst_pip_lo <- hierfstat::genind2hierfstat(gi_pip_local,pop=poploc)
  hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
  #replace all NAs with 0
  hfst_pip_lo <- hfst_pip_lo %>% replace(is.na(.), 0)
  # Remove rows where the pop name is unique - 
  # as the 'hierfstat::pairwise.WCfst'
  # function cannot handle the unique pop rows only  represented by only a single individual
  hfst_pip_lo <- hfst_pip_lo[hfst_pip_lo$pop %in% hfst_pip_lo$pop[duplicated(hfst_pip_lo$pop)],]
  hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
  #pw_Nei_pip_un <- hierfstat::pairwise.neifst(hfst_pip_un, diploid = FALSE)
  pw_Nei_pip_lo <- hierfstat::pairwise.neifst(hfst_pip_lo, diploid = FALSE)
  # make the data frame as a heat map
  df_pw_pip_lo <- as.data.frame(pw_Nei_pip_lo)
  #replace NAs with zeroes
  df_pw_pip_lo[is.na(df_pw_pip_lo)] = 0
  # assign it a different name 
  dlo <- df_pw_pip_lo
  #copy the dna.bin object
  pip3 <- pip
  #make the DNAbin object a genind object
  geni_pip3 <- adegenet::DNAbin2genind(pip3)
  #make the genind object a dataframe
  df_pip3 <- adegenet::genind2df(geni_pip3)
  #get the row names
  orig_rwnm <- row.names(df_pip3)
  # replace all NAs
  df_pip3 <- df_pip3 %>% replace(is.na(.), "-")
  #make the date frame a matrix and a DNAbin object again
  dnb_pip3 <- as.DNAbin(as.matrix(df_pip3))
  # get abbreviated location name
  ablo <- poploc
  library(MASS)
  d <- dist(dlo) # euclidean distances between the rows
  fit <- MASS::isoMDS(d, k=2) # k is the number of dim
  # make the fitted NMDS a data frame
  fmds <- as.data.frame(fit)
  # find the stress level and the moin and the max
  minstr <- min(fmds$stress)
  maxstr <- max(fmds$stress)
  strslvl <- ifelse(minstr==maxstr,maxstr, paste0(minstr,"-",maxstr) )
  strslvl <- round(strslvl,2)
  # make data frame with overall location names for the 4 letter location codes
  lscd <- c("BalS", "CasS", "CWAt", "FBog", "FKer", "GBus", "GHel", "GKie", 
            "GWis", "JLim", "JMar", "Medi", "NEAt", "Neth", "NWAt", "SBal", 
            "SSko", "USAG", "USAP", "USAW")
  lncd <- c("NE Atlantic", "Caspian Sea", "Central W Atlantic",
            "NE Atlantic", "NE Atlantic", "NE Atlantic", "NE Atlantic", "NE Atlantic", 
            "NE Atlantic", "NE Atlantic", "NE Atlantic", "Mediterranean", 
            "NE Atlantic", "NE Atlantic",
            "NW Atlantic", "NE Atlantic", 
            "NE Atlantic", "NW Atlantic", "NW Atlantic", "NW Atlantic")
  # and assemble in a data frame
  df_lnNm <- as.data.frame(cbind(lscd,lncd))
  # match to get overall location
  # fmds$ov.loc <- df_lnNm$lncd[match(rownames(fmds),df_lnNm$lscd)]
  fmds$ov.loc <- rownames(fmds)
  # 
  # get unique over locations
  locs3 <- unique(fmds$ov.loc)
  # count unique overall locations
  nlocs3 <- length(locs3)
  # make series of numbers to use for pch shapes
  pchnmbs <- rep(c(21:24),3)
  pchnmbs <- pchnmbs[1:nlocs3]
  # make series of colors to use for coloring fill
  colnmbs <- rep(c("white","black","black"),4)
  colnmbs <- colnmbs[1:nlocs3]
  
  colnmbs <- df_colf.locNm2$colf.locNm[match(locs3,df_colf.locNm2$locNm)]
  #bind together in a data frame
  df_pcl04 <- as.data.frame(cbind(locs3,pchnmbs,colnmbs))
  fmds$pchnmb <- df_pcl04$pchnmbs[match(fmds$ov.loc, df_pcl04$locs3)]
  fmds$flcnmb <- df_pcl04$colnmbs[match(fmds$ov.loc, df_pcl04$locs3)]
  # load the library
  library(ggplot2)
  library(ggrepel)
  nrwsfmds <- length(rownames(fmds))
  # plot it with ggplot
  p14 <- ggplot(fmds, aes(points.1, -points.2, label = rownames(fmds))) +
    geom_point(aes(shape= ov.loc, 
                   fill=ov.loc),
               color="black",
               size=3.0) +  
    geom_text(check_overlap = TRUE,
              vjust = -0.75) #+ theme_minimal() + xlab('') + ylab('') +
  # add text to the plot to indicate stress level
  p14 <- p14 +                               # Add text element to plot
    annotate("text", x = 0, y = -5.5, label = paste0("stress level: ",strslvl))
  p14 <- p14 +  scale_shape_manual(values = c(as.integer(pchnmbs))) +
    scale_fill_manual(values = c(colnmbs))
  p14 <- p14 + theme(panel.background = element_rect(fill = 'white', 
                                                     color = 'white'))#,
  # change label for legend - Notice that you need to change for all 3 variables
  # you called 'aes' in 'geom_jitter'
  p14 <- p14 + labs(fill='Location')
  p14 <- p14 + labs(color='Location')
  p14 <- p14 + labs(shape='Location')
  p14 <- p14 + xlab("NMDS1") + ylab("NMDS2")
  # add border around plot
  p14 <- p14 + theme(panel.border = element_rect(
    color = "black",
    fill = NA,
    size = 1.0))
  # change background of legend
  p14 <- p14 + theme(legend.key = element_rect(fill = "white"))
  # https://stackoverflow.com/questions/66102618/ggplot-increasing-the-distance-between-axis-labels-and-axis-ticks
  p14 <- p14 + theme(axis.text.x = 
                       element_text(color="black", 
                                    size=12.4,
                                    margin = margin(t = 2, 
                                                    unit = "mm")),
                     axis.text.y = 
                       element_text(color="black", 
                                    size=12.4,
                                    margin = margin(l = 0, r=2,
                                                    unit = "mm")))
  # collect the plots in a list
  lst.plt13b[[ng]] <- p14 
  
  #make filename to save plot to
  figname14 <- paste0("Fig13b_v",pnng,"_smpl_location_NMDS_",GnNm,"_DK_Germ_smpls.png")
  figname02 <- paste(wd00_wd05,"/",figname14,sep="")
  if(bSaveFigures==T){
    ggsave(p14,file=figname02,
           #width=210,height=297,
           width=210,height=(297*0.5),
           units="mm",dpi=300)
  }
  
  
  #_______________________________________________________________________________
  
  #_______________________________________________________________________________
  pip3 <- pip
  dnb_pip3 <- pip3
  #make the date frame a matrix and a DNAbin object again
  dnb_pip3 <- as.DNAbin(as.matrix(df_pip3))
  #
  library(MASS)
  # set pair wise deletion to TRUE to remove identical sequences
  d <- ape::dist.dna(dnb_pip3, model = "raw", pairwise.deletion = TRUE)
  #replace NAs with zero
  d[is.na(d)] <- 0
  #fit <- isoMDS(d, k=2) # k is the number of dim
  fmMDS<- metaMDS(d, k=2)
  # get the label categories for the fitted MDS
  lcNm0 <- labels(dnb_pip3)
  lcNm0 <- strsplit(lcNm0,"_")
  lcNm1 <- sapply(lcNm0, "[[", 1)
  lcNm2 <- sapply(lcNm0, "[[", 2)
  lcNm3 <- sapply(lcNm0, "[[", 3)
  # get the short location name
  shNm3 <-  lcNm2
  #fmMDS$species <- lcNm4
  fmMDS$species <- shNm3
  #row.names(fmMDS$points) <- lcNm4
  m.pnt.NMDS1 <- mean(fmMDS$points[,1])
  m.pnt.NMDS2 <- min(fmMDS$points[,2])+((max(fmMDS$points[,2]) - min(fmMDS$points[,2]))*(1/8))
  row.names(fmMDS$points) <- shNm3
  strslvl2 <- fmMDS$stress
  strslvl2 <- round(strslvl2,3)
  # define output file name for plot
  figname01 <- paste0("Fig13d_v",pnng,"_smpl_location_NMDS_plot",GnNm,"_DK_Germ_smpls.png")
  pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
  # set to save plot as png file 
  png(pthfignm01,
      width = 210, 
      height = 297,
      units="mm",
      res =300)
  #add extra space to the right of the plot
  par(mar=c(5, 5, 2, 2), xpd=FALSE)
  par(oma=c(0, 0, 0, 0))
  #reset this parameter
  par(mfrow = c(1, 1)) 
  #dev.off()
  ordiplot(fmMDS,type="n",cex.axis=1.8,cex.lab = 1.8)
  text(x = m.pnt.NMDS1, y = m.pnt.NMDS2, cex=0.8,                # Add text element
       paste0("stress level: ",strslvl2))
  orditorp(fmMDS,display="sites",
           air=0.2,cex=1.25)
  # end plot
  dev.off()
  #reset this parameter
  par(mfrow = c(1, 1)) 
  #_______________________________________________________________________________
  library(MASS)
  # set pair wise deletion to TRUE to remove identical sequences
  d <- ape::dist.dna(dnb_pip3, model = "raw", pairwise.deletion = TRUE)
  #replace NAs with zero
  d[is.na(d)] <- 0
  #fit <- isoMDS(d, k=2) # k is the number of dim
  fmMDS<- metaMDS(d, k=2)
  # make the fitted NMDS a data frame
  lcNm0 <- labels(dnb_pip3)
  lcNm0 <- strsplit(lcNm0,"_")
  lcNm1 <- sapply(lcNm0, "[[", 1)
  lcNm2 <- sapply(lcNm0, "[[", 2)
  lcNm3 <- sapply(lcNm0, "[[", 3)
  #combine to a data frame
  df_smp4 <- as.data.frame(cbind(lcNm1,lcNm3))
  #match to get sampling year
  df_smp4$smply <- lcNm3
  # # get sampling year and unqiue sample number from sequence name
  fmMDS$species <- df_smp4$smply
  
  m.pnt.NMDS1 <- mean(fmMDS$points[,1])
  m.pnt.NMDS2 <- min(fmMDS$points[,2])+((max(fmMDS$points[,2]) - min(fmMDS$points[,2]))*(1/8))
  lcNm5 <- df_smp4$smply
  row.names(fmMDS$points) <- lcNm5
  strslvl2 <- fmMDS$stress
  strslvl2 <- round(strslvl2,3)
  # define output file name for plot
  figname01 <- paste0("Fig13e_v",pnng,"_",GnNm,"_smpl_year_NMDS_plot_DK_Germ_smpls.png")
  pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
  # set to save plot as png 
  png(pthfignm01,
      width = 210, 
      height = 297,
      units="mm",
      res =300)
  #add extra space to the right of the plot
  par(mar=c(5, 5, 2, 2), xpd=FALSE)
  par(oma=c(0, 0, 0, 0))
  #reset this parameter
  par(mfrow = c(1, 1)) 
  ordiplot(fmMDS,type="n",cex.axis=1.8,cex.lab = 1.8)
  text(x = m.pnt.NMDS1, y = m.pnt.NMDS2, cex=0.8,                # Add text element
       paste0("stress level: ",strslvl2))
  orditorp(fmMDS,display="sites",
           air=0.2,cex=1.25)
  # end plot
  dev.off()
  #reset this parameter
  par(mfrow = c(1, 1)) 
  #_______________________________________________________________________________
  # end iteration over alignments
}

library(patchwork)
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p01t <- lst.plt13b[[1]] #+ labs(title = "a", face="bold")#,
p02t <- lst.plt13b[[2]] #+ labs(title = "b", face="bold")#,
#p03t <- p03 +   ggtitle('Plot 4', size =16) # (title = "c", face="bold")#,
#p03t <- p03 + theme(plot.title = element_text(size = 12, face = "bold"))
#p03t <- p03 + theme(plot.title = element_text(face="bold", size=18))

pA <-   p01t + 
  p02t +
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  # see : https://stackoverflow.com/questions/63434957/how-to-adjust-the-font-style-of-tags-with-plot-annotation-in-figures-assembled-w
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold', size =22)) +
  plot_layout(guides = "collect") #+
#plot_annotation(caption=pthinf01) #& theme(legend.position = "bottom")
#p
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig13b_v05_all_ITS1_and_ITS2_NMDS.png")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,
         width=210,
         height=297,
         units="mm",dpi=300)
}
#_______________________________________________________________________________
#_______________________________________________________________________________
# section 18 -  start - Make NMDS plots - only DK and Germany samples
#_______________________________________________________________________________