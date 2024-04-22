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
#read in packages and install libraries
#read in the 'ape' library
library(ape)
#read in the 'pegas' library
library(pegas)

#install packages if required
# gaston had issues installing. Returning the error: “failed to create lock directory”
# I looked here: https://stackoverflow.com/questions/14382209/r-install-packages-returns-failed-to-create-lock-directory
if(!require("gaston")){
  install.packages("Rcpp", dependencies = TRUE, INSTALL_opts = '--no-lock')
  install.packages("RcppParallel", dependencies = TRUE, INSTALL_opts = '--no-lock')
  install.packages("gaston", dependencies = TRUE, INSTALL_opts = '--no-lock')

}
library("gaston")
# The 'hierfstat' package has many functions for various pop gen 
# calculations
if(!require("hierfstat")){
  install.packages("hierfstat", dependencies = TRUE, INSTALL_opts = '--no-lock')
}
library("hierfstat")
library("ggrepel")

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

#https://www.rdocumentation.org/packages/strataG/versions/2.4.905
# if(!require(strataG)){
#   # make sure you have Rtools installed
#   if (!require('devtools')) install.packages('devtools')
#   # install from GitHub
#   devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
# }
# get ips packge to write dnabin to fasta files
if(!require("ips")){
  install.packages("ips", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("ips")
}
library("ips")

# to be able to install 'dartR' 'rjags' must be installed
# first. 'rjags' requires that 'jags' is installed.
# 'jags' can be installed in the terminal
# see this problem
# https://github.com/r-hub/rhub/issues/296
# in a terminal type:
# $ sudo apt install jags
if(!require("rjags")){
  install.packages("rjags", dependencies = TRUE, INSTALL_opts = '--no-lock')
}
library("rjags")
# to install 'dartR' it is rquired that 'SNPRelate' is
# installed.
# this website explains how "SNPRelate" can be installed
# https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html
if(!require("SNPRelate")){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", force =T)
  BiocManager::install("SNPRelate", force =T)
}
library("SNPRelate")
#install the 'dartR' package
if(!require("dartR")){
  install.packages("dartR", dependencies = TRUE, INSTALL_opts = '--no-lock')
}
library("dartR")


if(!require("ecodist")){
  install.packages("ecodist", dependencies = TRUE, 
                   force=T,INSTALL_opts = '--no-lock')
}
library("ecodist")

if(!require("hierfstat")){
  install.packages("hierfstat", dependencies = TRUE, INSTALL_opts = '--no-lock')
}
library("hierfstat")
# get textclean packge to replace non ASCII characters
if(!require("textclean")){
  install.packages("textclean", dependencies = TRUE, INSTALL_opts = '--no-lock')
  library("textclean")
}
library("textclean")
#The 'tidyverse' package is required for rearranging the data frames before
# making the colored tables
if(!require("tidyverse")){
  install.packages("tidyverse", dependencies = TRUE, INSTALL_opts = '--no-lock')
}
library("tidyverse")
# the 'pals' enables you to make color gradients - here used for colored tables
if(!require("pals")){
  install.packages("pals", dependencies = TRUE, INSTALL_opts = '--no-lock')

}
library("pals")
# strataG package can calculate Tajimas D,

#But the 'strata' package is not available for R v3.6.2
#The strataG package works for R v4.0.2
if(!require("strataG")){
  #install.packages("strataG", dependencies = TRUE, INSTALL_opts = '--no-lock')
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  library("devtools")
  # 'apex' package is required for 'strataG'
  install_version("apex", "1.0.4")
  library(apex)
  # install from GitHub
  devtools::install_github('ericarcher/strataG', build_vignettes = F)
}
library("strataG")

#install mmod package
if(!require("mmod")){
  install.packages("mmod", dependencies = TRUE)
}
library("mmod")

if(!require(sf)){
  install.packages("sf")
}
library(sf)

# you need the package 'scico' for the package 'paletteer
if(!require(scico)){
  install.packages("scico")
  }
library(scico)
# get the 'paletteer' package to make color palettes
if(!require(paletteer)){
  install.packages("paletteer")
  }
library(paletteer)
#install package if needed
if(!require(rnaturalearth)){
  install.packages("rnaturalearth")
  #install.packages("rnaturalearth", repos = "http://packages.ropensci.org", type = "source")
  
}
library(rnaturalearth)

#read in libraries
library("apex")
library("adegenet")
library("pegas")
library("mmod")
library("hierfstat")
library("tidyverse")
library("pals")
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
# # get the biogeo package 
if(!require(biogeo)){
  # check if the package is available
  ap <- available.packages()
  "biogeo" %in% rownames(ap)
  # https://stackoverflow.com/questions/25721884/how-should-i-deal-with-package-xxx-is-not-available-for-r-version-x-y-z-wa
  # it is out of date. So install it from the archive
  # install.packages("biogeo")
  library(remotes)
  install_version("biogeo", "1.0")
}
library(biogeo)
# get the BiocManager package
if (!require("BiocManager", quietly = TRUE))
  if(!require(BiocManager)){
    install.packages("BiocManager")
  }
library(BiocManager)
# get the htmlTable package
if (!require("htmlTable", quietly = TRUE))
  if(!require(htmlTable)){
    install.packages("htmlTable", force=T)
  }
library(htmlTable)

# get the scatterpie package
if (!require("scatterpie", quietly = TRUE))
  if(!require(scatterpie)){
    install.packages("scatterpie", force=T)
  }
library(scatterpie)

# get the tableHTML package
if (!require("tableHTML", quietly = TRUE))
  if(!require(tableHTML)){
    install.packages("tableHTML", force=T)
  }
library(tableHTML)


# get the ggtree package
if(!require(ggtree)){
  BiocManager::install("ggtree")
}
library(ggtree)
# get the textclean package
if(!require(textclean)){
  install.packages("textclean")
}
library(textclean)
library(ggtree)
library("ggplot2")
library("ggtree")
# get package
if(!require("kableExtra")){
  # see https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html#Table_Styles
  #For dev version
  install.packages("devtools")
  devtools::install_github("haozhu233/kableExtra")
  
}
library("kableExtra")


 # I installed the Bioconductor package  DECIPHER
# some of the dependencies that came with the DECIPHER
# package appear to block out some of the base functions
# even though I tried to call the namespace for each for the 
# functions I could not call the base functions, and I 
# could not uninstall single Bioconductor packages.
# Instead I found this webiste
#https://support.bioconductor.org/p/7071/
# That suggested that I can remove everything
# so in my R library I decided to run 'rm -rf' in a terminal
# to get rid of everything
# Unfortunately, I need to reinstall everything. 
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
inpf01 <- "algn_Mnelei_v13.fa"
inpf01 <- "algn_Mnelei_18s_16.fas.aligned.fasta"
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
# define a filename to write the data frame to as a csv file
wd00_wd05_flnm3 <- paste(wd00_wd05,"/df_lN02.csv",sep="")
# write the data frame as a csv file
write.csv(df_lN02,file=wd00_wd05_flnm3)

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
  c("FBo",
    "FKe",
    "JMa",
    "GKi",
    "GWi",
    "JLi",
    "JLo",
    "GHe",
    "GBu",
    "SBa",
    "GWi",
    "Ssk")
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


df_lN02$accession_nmb

# write the df_lN02 data frame to a csv file
write.csv(df_lN02,paste0(wd00_wd05,"/df_lN02.csv"))
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
nrow(pip)
ncol(pip)

#make the DNAbin object a genind object
geni_pip <- adegenet::DNAbin2genind(pip)
#make the genind object a dataframe
df_pip <- adegenet::genind2df(geni_pip)
# replace with gsub to get unique rows only
rwNm <- gsub("_$","",rownames(df_pip))
unqrwNm <- unique(rwNm)
df_pip <- df_pip[!is.na(match(rownames(df_pip),unqrwNm)),]

#get the row names
orig_rwnm <- row.names(df_pip)
#substitute in the row names
mdf_rnm01 <- gsub("_consensus_sequence","",orig_rwnm)
mdf_rnm02 <- gsub("Mnelei","Mne_lei", mdf_rnm01)
mdf_rnm03 <- gsub("Mnemiopsis_leidyi","Mne_lei", mdf_rnm02)
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
df_pip$rwnm04 <- gsub("_","",df_pip$rwnm04)
df_pip$rwnm05 <- gsub("_","",df_pip$rwnm05)
df_pip$rwnm06 <- gsub("_","",df_pip$rwnm06)
df_pip$rwnm07 <- gsub("_","",df_pip$rwnm07)
df_pip$rwnm09 <- df_pip$rwnm05
df_pip$rwnm10 <- df_pip$rwnm06
#match accesion numbers with sample location and sample year
df_pip$rwnm09[grepl("[0-9]{5}",df_pip$rwnm09)] <- df_lN02$location[match(df_pip$rwnm09[grepl("[0-9]{5}",df_pip$rwnm09)],df_lN02$accession_nmb)]
df_pip$rwnm10[grepl("[0-9]{5}",df_pip$rwnm10)] <- df_lN02$smplyear2[match(df_pip$rwnm10[grepl("[0-9]{5}",df_pip$rwnm10)],df_lN02$accession_nmb)]
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
df_pip$rwnm08 <- paste(df_pip$rwnm07,df_pip$rwnm05,df_pip$rwnm06,df_pip$rwnm04,sep="_")
df_pip$rwnm08 <- gsub("NJyllandLimfjordLogstoer","NJyllandLimfjord",df_pip$rwnm08)
#replace the row names in the data frame
row.names(df_pip) <- mdf_rnm03
row.names(df_pip) <- df_pip$rwnm08
df_pip$rwnm02 <- row.names(df_pip)
# subset to only include the rows that match
df_pip02 <- subset(df_pip, grepl("Mnelei", rwnm02))
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

nrow(df_pip02)
ncol(df_pip02)
#make the date frame a matrix and a DNAbin object again
dnb_pip02 <- as.DNAbin(as.matrix(df_pip02))
# define a filename to write the data frame to as a csv file
wd00_wd05_flnm2 <- paste(wd00_wd05,"/dnb_pip02_df.csv",sep="")
# write the data frame as a csv file
write.csv(df_pip02,file=wd00_wd05_flnm2)
# copy the data frame to a new object
pip <- dnb_pip02
#________________________________________________________________
#start - Plot haplotype network:
#________________________________________________________________

#create haplotypes from dna.bin
pipHaps <- pegas::haplotype(pip)
#view haplotype 
#as.matrix(pipHaps)
#prepare hpt table
ind.hap<-with(
	utils::stack(setNames(attr(pipHaps, "index"), rownames(pipHaps))),
	table(hap=ind, pop=rownames(pip)[values]))
#make it a dataframe
df_ihpt01 <- as.data.frame(ind.hap)
#limit to include only 'Freq' that equals 1
df_ihpt02 <- df_ihpt01[df_ihpt01$Freq == 1,]	
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
#unique(smplloca)
accnms <- gsub("Mnelei","",sapply(locations, "[[", 1))
#df_lN02$smplyear[match(accnms,df_lN02$accession_nmb)]

# only year
smplye <- sapply(locations, "[[", 3)
# year and month
smplym <- sapply(locations, "[[", 4)
#class(smplym)
# see the first header of this list of characters
#head(smplym)
#replace "Germany:Maasholm" with,"Germany, Kiel Fjord"
smplloca <- gsub("Germany:Maasholm","Germany, Kiel Fjord" ,smplloca)
smplloca <- gsub("AtlanticOcean:NWAtlantic","NW Atlantic" ,smplloca)
smplloca <- gsub("BalticSea","Baltic Sea" ,smplloca)
smplloca <- gsub("NEAtlantic" ,"NE Atlantic",smplloca)
smplloca <- gsub("CentralWAtlantic" ,"CW Atlantic",smplloca)
smplloca <- gsub("CaspianSea" ,"Caspian Sea",smplloca)
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
# add together the 2 sampling locations that are close to each other
#df_nhpt02.loc$`Germany, Kiel Fjord` <- df_nhpt02.loc$`Germany, Kiel Fjord`+df_nhpt02.loc$`Germany:Maasholm`
# and then remove the redundant column
#df_nhpt02.loc$`Germany:Maasholm` <- NULL

#dev.off()
#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
# turn of previous plot
#dev.off()
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 1.2, cex = 0.1, pie = new.hap.smplloc, 
     show.mutation = 0, threshold = 0, labels(FALSE))
#add a legend to the plot
legend("topright",colnames(new.hap.smplloc), 
       col=rainbow(ncol(new.hap.smplloc)), 
       pch=19, ncol=1)


plot(pipNet, size = (sqrt((attr(pipNet,"freq"))/pi)), 
     scale.ratio = 1.2, cex = 1.1, pie = new.hap.smplloc, 
     show.mutation = 0, threshold = 0, labels(T))

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
     scale.ratio = 1.6, cex = 1.1, pie = new.hap.smplloc, 
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
     scale.ratio = 1.6, cex = 1.1, pie = new.hap.smplye, 
     show.mutation = 0, threshold = 0, labels(FALSE))
#add a legend to the plot
legend("topleft",colnames(new.hap.smplye), 
       col=rainbow(ncol(new.hap.smplye)), 
       pch=19, ncol=1, cex=0.6)

#end plot area
par(op)
# end svg file to save as
dev.off()  

plot(pipNet, size = (sqrt((attr(pipNet,"freq"))/pi)), 
     scale.ratio = 2.6, cex = 1.1, pie = new.hap.smplye, 
     show.mutation = 2, threshold = 0, labels(TRUE))
#________________________________________________________________
# end - Plot haplotype network:
#________________________________________________________________
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
plot(pipNet, size = sqrt(attr(pipNet,"freq")/pi), 
     scale.ratio = 1.6, cex = 1.1, 
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
plot(pipNet, size = sqrt(attr(pipNet,"freq")/pi), 
     scale.ratio = 1.6, cex = 1.1, pie = new.hap.smplye, 
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
mtr_dd_pip <- ape::dist.dna(pip)
dmtr_dd_pip <- as.data.frame(as.matrix(mtr_dd_pip))
# replace NAn with zero
mtr_dd_pip[is.na(mtr_dd_pip)] <- 0
tre_pip <- ape::njs(mtr_dd_pip)
#tre_pip <- ape::njs(mtr_dd_pip)

#plot(tre_pip)
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

p <- ggtree(tre_pipr, right = TRUE,
            options(ignore.negative.edge=TRUE)) + 
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
flnm <- c(paste("Fig02_NJ_tree_",inp.f.fnm,".pdf",  sep = ""))
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
pth_inpf05 <- paste(wd00_wd05,"/",inpf04,sep="")
#write the output file
ips::write.fas(dnb_pip02,pth_inpf05)

#________________________________________________________________


#________________________________________________________________
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#define input file as variable
#inpf01 <- "algn_Mnelei_18s_09.fas"
pth_inpf01 <- paste(wd00_wd01,"/",inpf01,sep="")
#Test the 'apex' package :
#read FASTA as multi.dna
mltdn_pip <- apex::read.multiFASTA(pth_inpf01)
#check what kind of object it is
#class(mltdn_pip)
#view as plotted alignment
plot(mltdn_pip, cex = 0.2)

#dev.off()
#string split the character labels by underscore
# in the multidna object
# rowbind the nested lists to a dataframe
df_pip02 <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(mltdn_pip@labels),
                                  "_")))
#check object
#class(df_pip02)
if(length(colnames(df_pip02))==7){
#rename columns names
colnames(df_pip02) <- c("gnnm1","spnm1","spnm2",
                        "gnnm2","spnm3","spnm3",
                        "gnnm3")}
# see the header rows of df
#subset multidna object - to exclude species that are not Mnemiopsis leidyii
mltdn_pip02 <- mltdn_pip[grepl("Mnelei",mltdn_pip@labels)]
#view as plotted alignment
plot(mltdn_pip02, cex = 0.2)
#dev.off()
#_______________________________________________________________________________
#_______________________________________________________________________________

# ------------- Options -------------
#set fraction OK for all other files 
FracOK <- 0.90
# if the input file is a sequence file with ets2 gene fragment, then set
# the fraction higher
if (grepl("ets2",pth_inpf01)==T)
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
#_______________________________________________________________________________
#define input file as variable
#Test the 'apex' package :
#read FASTA as multi.dna
mltdn_pip03 <- apex::read.multiFASTA(pth_inpf01)
# see the trimmed alignment
plot(mltdn_pip03, cex = 0.2)
# see the names of the sequences
#mltdn_pip03@labels
# rowbind the nested lists to a dataframe
df_pip03 <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(mltdn_pip03@labels),
                                  "_")))
#check object
#class(df_pip03)
df_pip03$X1[grepl("[0-9]{5}",df_pip03$X3)] <- df_pip03$X3[grepl("[0-9]{5}",df_pip03$X3)]
loc_for_accnno <- df_lN02$location[match(df_pip03$X3[grepl("[0-9]{5}",df_pip03$X3)],df_lN02$accession_nmb)]
year_for_accnno <- df_lN02$smplyear2[match(df_pip03$X3[grepl("[0-9]{5}",df_pip03$X3)],df_lN02$accession_nmb)]
df_pip03$X4[grepl("[0-9]{5}",df_pip03$X1)] <- year_for_accnno
df_pip03$X2[grepl("[0-9]{5}",df_pip03$X1)] <- loc_for_accnno
df_pip03$X3[grepl("[0-9]{5}",df_pip03$X1)] <- loc_for_accnno
#head(df_pip03,9)
#rename columns names
colnames(df_pip03) <- c("unqnm","locality2","locality","year","month")
# replace the 'Germany:Maasholm' location name with "NGermanyKielFjord"
# since these two locations are very close to each other
df_pip03$locality[df_pip03$locality=="Germany:Maasholm"] <- "NGermanyKielFjord"
#df_pip03$locality3 <- paste0(df_pip03$locality,"_",df_pip03$locality2)
#make the multidna object a genind object
gi_pip03 <- apex::multidna2genind(mltdn_pip03)
#Make the genind object a data frame
df_geipip03 <- genind2df(gi_pip03)
#replace all NAs with "-"
df_geipip04 <- df_geipip03 %>% replace(is.na(.), "-")
df_geipip03 <- df_geipip04
#make the df object a matrix and make this a dnabin object
dnb_pip03 <- as.DNAbin(as.matrix(df_geipip04))
dp3lo <- df_pip03$locality
dp3lo <- gsub("MecklenburgerBuchtWismarBucht","MecklenburgerBucht",dp3lo)
dp3lo <- gsub("USA:GalvestonBay","USA:WoodsHole",dp3lo)
dp3lo <- gsub("USA:Panacea","USA:WoodsHole",dp3lo)
dp3lo <- gsub("NGermanyKielFjord","KielFjord",dp3lo)
dp3lo <- gsub("Germany:KielFjord","KielFjord",dp3lo)
dp3lo <- gsub("Germany:Maasholm","KielFjord",dp3lo)
dp3lo <- gsub("USA:WoodsHole","AtlanticOcean:NWAtlantic",dp3lo)
dp3lo <- gsub("Germany:Helgoland","NSeaHelgolandRds",dp3lo)
dp3lo <- gsub("Netherlands","NEAtlantic",dp3lo)
df_pip03$locality3 <- dp3lo
# copy column in data frame
df_clo03$locality15 <- df_clo03$locality4
# make two vectors. One with old names and one with new names
olNm <- c("Samsøe, Ballen","Sealand, Skovshoved","Jutland, Mariager Fjord","Funen, Kerteminde","Funen, Bogense","Jutland, Limfjord" ,"North Sea, Helgoland Roads","Germany, Büsum","Germany, Wismar Bight","Germany, Kiel Fjord","Netherlands")
#olNm <- c("Samsøe, Ballen","Sealand, Skovshoved","Jutland, Mariager Fjord","Funen, Kerteminde","Funen, Bogense","Jutland, Limfjord" ,"North Sea, Helgoland Roads","Germany, Büsum","Germany, Wismar Bight","Germany, Kiel Fjord","Netherlands")
neNm <- c("Ballen","Skovshoved" ,"Mariagerfjord", "Kerteminde" ,"Bogense","Loegstoer","NSeaHelgolandRds","WaddenSeaBussumHaupstr","MecklenburgerBucht","KielFjord" ,"NEAtlantic")
# bind the vectors together as columns in a dataframe
df_olNm <- as.data.frame(cbind(olNm,neNm))
#iterate over rows in this data frame, and call elements from rows in the data frame
# and use this to substitute in the column
for (i in seq(1:nrow(df_olNm)))
{
  olr <- df_olNm$olNm[i]
  ner <- df_olNm$neNm[i]
  df_clo03$locality15 <- gsub(olr,ner,df_clo03$locality15)
}
# match to get abbreviated location name
df_pip03$locality4 <- df_clo03$locality8[match(df_pip03$locality3,df_clo03$locality15)]
#paste together to get a new string to use as sequence names
df_pip03$nwseqNm <- paste0(df_pip03$unqnm,"_",df_pip03$locality4,"_",df_pip03$year)
# original sequence names
df_pip03$oriseqNm <- as.character(mltdn_pip03@labels)

library(tidyverse)
library(ape)
#install.packages("ape")
if(!require(ape)){
  install.packages("ape")
  library(ape)
}
# get the working dir
prwd<- getwd()
# set working directory to be able to read in functions
setwd(paste0(wd00,"/suppmat03_Rcodes"))
source("Rfunction_DNA_sequence_subset.R")
source("Rfunction_ReadFasta.R")
setwd(prwd)
# ------------- Load data using ReadFasta function -------------
df <- ReadFasta(pth_inpf01)
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


filename<-paste0(substr(inpf01,1,nchar(inpf01)-4),"_cropped.csv")
folder_out <- wd00_wd05
write.table(df_out,file=paste0(folder_out,filename),row.names=F,col.names=T,sep=",",quote=F)
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
fnm_toread <- pth_inpf01
#fnm_toread <- inpf01
# read in the alignment
al_dt01 <- ape::read.dna(fnm_toread,format="fasta", as.matrix=TRUE)
# get the upper and lower limits to cut the alignment by
uppbp <- as.numeric(gsub(".*-","",sListOK))
lowbp <- as.numeric(gsub("-.*","",sListOK))
if (uppbp==Inf){
  uppbp <- ncol(al_dt01)
lowbp <- 1}
#class(al_dt02)
#trim the alignment
al_dt02 <- al_dt01[,lowbp:uppbp]

# Notice that this line 'labels(al_dt02) <- df_pip03$nwseqNm[match(labels(al_dt02),df_pip03$oriseqNm)]'
# will not work
row.names(al_dt02) <- df_pip03$nwseqNm[match(row.names(al_dt02),df_pip03$oriseqNm)]
# define output directory
outdir <- wd00_wd05
# define a file to write to
fnm01_1 <- gsub("_fastafiles\\.","_",inpf01)
fnm01_2 <- gsub("fas\\.","trimmed\\.",fnm01_1)
fnm01_2 <- gsub(".fasta.fas","\\.trimmed\\.fas",fnm01_1)
fnm01_2 <- gsub("18s_16","18s_17",fnm01_2)
#fnm01_3 <- gsub("triplefin_out04_","triplefin_out06_",fnm01_2)
fnm02 <- gsub("\\.txt","",fnm01_2)
#check if row names are as they are supposed to be
#row.names(al_dt02)
fnm_towrite2 <- paste0(outdir,"/",fnm02)
fnm_towrite <- paste0(outdir,"/",fnm02)
#write the trimmed alignment as a fasta file
ape::write.dna(al_dt02, file=fnm_towrite, format="fasta", nbcol=-1, colsep="")

#get the local sampled label names
ownsmplN2<- row.names(al_dt01)[grepl("Mnelei",row.names(al_dt01))]
lngNm <- ownsmplN2
# use substitute to get hte different elements
MnsmNo <- gsub("^(Mnelei[0-9]{+})_.*","\\1",ownsmplN2)
ownsmplN3 <- gsub("^(Mnelei[0-9]{+})_.*_([0-9]{4}.*)$","\\2",ownsmplN2)
smplyear4 <- gsub("^(.*)_(.*)$","\\1",ownsmplN3)
smplmnt4 <- gsub("^(.*)_(.*)$","\\2",ownsmplN3)

# make data frame with dates and months and years for own samples
df_ownsmpldt <- as.data.frame(cbind(MnsmNo,smplyear4,smplmnt4,lngNm))

# also make a data frame that holds the new names and the old names
# the abbreviations
smplnwNm <- row.names(al_dt02)
sAbbrv <- gsub("^(.*)_(.*)_(.*)$","\\2",row.names(al_dt02))
sNoAbbrv <- gsub("^(.*)_(.*)_(.*)$","\\1",row.names(al_dt02))
smplolNm <- row.names(al_dt01)
df_olnwNm <- as.data.frame(cbind(smplolNm,smplnwNm,sAbbrv,sNoAbbrv))

# also write a faste file with only the sequences that are to be submitted to
# NCBI GenBank
fnm_towrite <- gsub("18s_17","18s_18",fnm_towrite)
#fnm01_3 <- gsub("triplefin_out04_","triplefin_out06_",fnm01_2)
fnm02 <- gsub("\\.txt","",fnm01_2)
#check if row names are as they are supposed to be
#row.names(al_dt02)

# make it a bioseq tibble
sbt_al_dt02 <- bioseq::as_tibble.DNAbin(al_dt02)
# get the index number for the rows with Mneli in the label name
rw.indx_sbs_al_dt02 <- which(grepl("Mnelei",sbt_al_dt02$label) )
# use this index numbering to subset the tibble data frame
sbt_al_dt03 <- sbt_al_dt02[rw.indx_sbs_al_dt02,]

# get the number of sequences to be able to iterate over them
nosq <- nrow(sbt_al_dt03)

wd_NCBI <- "NCBI_seq_submission_for_Mnemiopsis"
outdir_NCBI <- paste0(wd00,"/",wd_NCBI)
outfile_NCBI <- paste0(outdir_NCBI,"/",fnm02)
# remove any previous versions of the file you are about to save
unlink(fnm_towrite)
unlink(outfile_NCBI)



infl01 <- paste0(outdir_NCBI,"/Mnemiopsis_18s_2024apr15_v01.fas")
bsq_Mnelei_18s <- bioseq::read_fasta(infl01)
# make it a tible bioseq data frame
tib.bsq_Mnelei_18s <- bioseq::as_tibble.DNAbin(bsq_Mnelei_18s)


# strsplit names and turn into characters
# and rowbind the nested lists to a dataframe
df_lbl05 <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(tib.bsq_Mnelei_18s$label),
                                  "_")))

# modify columns names
colnames(df_lbl05) <- c("smplNm",
                        "loc01",
                        "loc02",
                        "yer",
                        "mnth")
# rpelace if there is no sample month
df_lbl05$mnth[grepl("Mnel",df_lbl05$mnth)] <- "Sep"

# the sequence name identifier needs to be shorter in order for NCBI GenBank 
# BankIt to accept the fasta file
tib.bsq_Mnelei_18s$label <- df_lbl05$smplNm
NCBIBankit.fnm_towrite <- paste0(wd00 ,"/",
                          wd_NCBI,"/",
                          "Bankit_NCBI_seqs_submiss_file_Mnelei_18s.fasta")
# remove any previous versions of the file you are about to save to
unlink(NCBIBankit.fnm_towrite)

NCBIBankit.fnm_feature.f <- paste0(wd00 ,"/",
                                 wd_NCBI,"/",
                                 "Bankit_NCBI_feature_file_Mnelei_18s.txt")

# remove any previous versions of the file you are about to save to
unlink(NCBIBankit.fnm_feature.f)
#fileConn <- file(NCBIBankit.fnm_feature.f)
# write("NA_line", 
#       file=fileConn,
#       append=TRUE)
# close(fileConn)

fileConn <- file(NCBIBankit.fnm_feature.f,open="a")
# iterate over every single line to write out a fasta file
for (e in seq(1,nrow(tib.bsq_Mnelei_18s),1))
{print(e)
  
  lblN <- tib.bsq_Mnelei_18s$label[e]
  seqN <- tib.bsq_Mnelei_18s$sequence[e]
  # use sub to replace the first occurence of multiple Ns
  seqN.noN<- 
    sub("^[N]+{1}", "",seqN)
  # then use sub to replace the  occurence of multiple Ns from the other end
  seqN.noN <- sub("[N]+{1}$", "",seqN.noN)
  # use this funciotn to reverse the seq toi check how it ended up looking
  # https://stackoverflow.com/questions/13612967/how-to-reverse-a-string-in-r
  strReverse <- function(x)
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  strReverse(c(seqN.noN))
  # replace the previous version of the sequence
  seqN <- seqN.noN
  seqinr::write.fasta( seqN.noN,
                        lblN,
                       NCBIBankit.fnm_towrite,
                        open="a")
  # get end and start nucleotide position of the annotations
  idntf_18s_end <- (stringr::str_locate(seqN, "GATCATTA"))[2]
  idntfITS1_str <- (stringr::str_locate(seqN, "ACGAATCCAA"))[1]
  idntfITS1_end <- (stringr::str_locate(seqN, "CTAAAAGCGAA"))[2]
  idntf5.8s_str <- (stringr::str_locate(seqN, "CAACTTTAAACGG"))[1]
  idntf5.8s_end <- (stringr::str_locate(seqN, "GAGCGTCGTTT"))[2]
  idntfITS2_str <- (stringr::str_locate(seqN, "CTCACATCCCAT"))[1]
  
  if (T==is.na(idntf_18s_end)){idntf_18s_end <- idntfITS1_str-1}
  
  if (T==is.na(idntfITS1_str)){idntfITS1_str <- idntf_18s_end+1}
  if (T==is.na(idntfITS1_end)){idntfITS1_end <- idntf5.8s_str-1}
  
  if (T==is.na(idntf5.8s_str)){idntf5.8s_str <- idntfITS1_end+1}
  if (T==is.na(idntf5.8s_end)){idntf5.8s_end <- idntfITS2_str-1}
  
  if (T==is.na(idntfITS2_str)){idntfITS2_str <- idntf5.8s_end+1}
  
  # get the length of the seq
  seqL <- nchar(seqN)
  # make a feature text
Ft_tx <- paste0(">Feature ",lblN,"
  
<1    ",idntf_18s_end,"     gene
                            gene          18Sr
<1    ",idntf_18s_end,"     mRNA
                            product       18Sr ribosomal RNA

",idntfITS1_str,"    ",idntfITS1_end,"    gene
                                          gene          its1
",idntfITS1_str,"    ",idntfITS1_end,"    CDS
                                          product       internal transcribed spacer 1
",idntfITS1_str,"    ",idntfITS1_end,"    mRNA
                                          product       internal transcribed spacer 1
                                          transl_table     1
                                          
",idntf5.8s_str,"    ",idntf5.8s_end,"    gene
                                          gene          5.8Sr
",idntf5.8s_str,"    ",idntf5.8s_end,"    mRNA
                                          product       5.8S ribosomal RNA
                                          
",idntfITS2_str,"    >",seqL,"     gene
                                   gene          its2
",idntfITS2_str,"    >",seqL,"     CDS
                                   product       internal transcribed spacer 2
",idntfITS2_str,"    >",seqL,"     mRNA
                                   product       internal transcribed spacer 2
                                   transl_table     1
                ")
  

                  #print(Ft_tx)
                  write(Ft_tx, 
                        file=fileConn,
                        append=TRUE)
                  
  }
close(fileConn)


# and rowbind the nested lists to a dataframe
df_lbl_d03 <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(sbt_al_dt03$label),
                                  "_")))
# modify columns names
colnames(df_lbl_d03) <- c("smplNm",
                        "loc01",
                        "yer")
#
hd_src_tbl <- c( "Sequence_ID",
                 "Organism",
                  "Collected_by",
                  "Collection_date",
                  "Country",
                  # "Isolation_source",
                  # "Isolate",
                  "Lat_Lon")
                  # "Specimen_voucher")
# get  population location names
collcloc03 <- strsplit(as.character(sbt_al_dt03$label), "_")
collcloc03 <- as.factor(sapply(collcloc03, "[[", 2))
cntr <- rep("Denmark",length(collcloc03))
Gidx <- which(grepl("^G",collcloc03))
cntr[Gidx] <- "Germany"
clctBy <- rep("Florian Luskow",length(collcloc03))

tib.bsq_Mnelei_18s_02 <- bioseq::as_tibble.DNAbin(bsq_Mnelei_18s)
# replace collectors name for Ballen samples
clctBy[grepl("Ballen",tib.bsq_Mnelei_18s_02$label)] <- "Steen Knudsen"

# get the latitude and the longitude
clo03.lat2rnd <- round(as.numeric(df_clo03$dec_lat2[match(collcloc03,df_clo03$locality2)]),digits = 2)
clo03.lon2rnd <- round(as.numeric(df_clo03$dec_lon2[match(collcloc03,df_clo03$locality2)]),digits = 2)
pos_Lat_lLon_for_NCBI <- paste0(clo03.lat2rnd," N ",clo03.lon2rnd," E")
# get the sampling year
collcyer03 <- strsplit(as.character(sbt_al_dt03$label), "_")
collcyer03 <- as.factor(sapply(collcyer03, "[[", 3))
# get all elements for the NCBI 'Sample Source Modifiers Table file'
Sequence_ID <- tib.bsq_Mnelei_18s$label
Collected_by  <- clctBy
Collection_year <- collcyer03 
Country <- cntr
Lat_Lon <- pos_Lat_lLon_for_NCBI
Organism <- "Mnemiopsis leidyi"
# combine as columns to get a NCBI 'Sample Source Modifiers Table file'
df_NCBI_src_tbl <- cbind(Sequence_ID,
                         Organism,
                          Collected_by,
                          Collection_year, 
                          Country,
                          Lat_Lon,
                         df_lbl_d03)
#make it a data frame
df_NCBI_src_tbl <- as.data.frame(df_NCBI_src_tbl)
df_NCBI_src_tbl$Collection_date <- as.numeric(as.character(collcyer03))
# also get the- collection month
df_NCBI_src_tbl$Collection_month <- df_lbl05$mnth[match(df_NCBI_src_tbl$smplNm,df_lbl05$smplNm)]
# make a collection date
df_NCBI_src_tbl$Collection_date <- paste0(df_NCBI_src_tbl$Collection_year,"-", 
df_NCBI_src_tbl$Collection_month,"-01") 
# only keep required columns
df_NCBI_src_tbl <- df_NCBI_src_tbl[hd_src_tbl]
# make a file name
outfile_NCBI_tab_src <- "Mnelei_sample_data_tab_src_delim_for_NCBI_v01.txt"
outfile_NCBI_tab_src <- paste0(wd00,"/",wd_NCBI,"/",outfile_NCBI_tab_src)

# write a tabe separated file for NCBI 'Sample Source Modifiers Table file'
write.table(df_NCBI_src_tbl,file=outfile_NCBI_tab_src, 
            quote = FALSE,
             row.names = FALSE,
            sep = "\t")
# read in file with annotations
infl02 <- "Mnemiopsis_18s_2024apr15_annotation_table.txt"
wd_NCBI.infl02 <- paste0(wd00,"/",
                         wd_NCBI,"/",
                         infl02)
df_annt01 <- read.table(wd_NCBI.infl02, header=T,sep=",")

#colnames(df_annt01)
# iterate ove characters to replace
elem_to_rplc <- c(">" ,"," ,"<")
for (e in elem_to_rplc)
{
  # use gsub for multiple colunms on the element to replace
  df_annt01[c("Minimum",
            "Maximum",
            "Length")] <- gsub(e,"", unlist(
                              df_annt01[c("Minimum",
                                          "Maximum",
                                          "Length")]))
# make the columns numeric
    df_annt01[c("Minimum",
              "Maximum",
              "Length")] <- sapply(df_annt01[c("Minimum",
                                     "Maximum",
                                     "Length")],as.numeric)
}


# file name for the text file with accession numbers
GBaccsF <- "GenBank_accession_numbers_for_the_Mnelei_submitted_sequences_2024apr17.txt"
# After having obtained the accession numbers from NCBI
# I got an email back with accession numbers
inf_accNos <- paste0(wd00,"/",
                     wd_NCBI,"/",
                     GBaccsF)
# read in the file using 'stringr'
library(stringr) # For str_trim 
# Read string data and split into data frame
txAcNdat <- readLines(inf_accNos)
# only keep rows that begin with 'SUB'
txAcNdat <- txAcNdat[grepl("^SUB",txAcNdat)]
# make it a data frame
df_accNos <- as.data.frame(do.call(rbind, strsplit(txAcNdat, 
                  split=" ")), stringsAsFactors=FALSE)
# only keep columns 1,2 and 5 
df_accNos <- df_accNos[,c(1,2,5)]
# change the column names to something more meaningful
colnames(df_accNos) <- c("SUBNo","MneleiAbbr","AccNo") 

#_______________________________________________________________________________
#_______________________________________________________________________________

# define a file to write to
fnm01_1 <- gsub("_fastafiles\\.","_",inpf01)
fnm01_2 <- gsub("fas\\.","trimmed\\.",fnm01_1)
fnm01_2 <- gsub(".fasta.fas","\\.trimmed\\.fas",fnm01_1)
fnm01_2 <- gsub("18s_16","18s_17",fnm01_2)
#fnm01_3 <- gsub("triplefin_out04_","triplefin_out06_",fnm01_2)
fnm02 <- gsub("\\.txt","",fnm01_2)
#check if row names are as they are supposed to be
#row.names(al_dt02)
fnm_towrite2 <- paste0(outdir,"/",fnm02)

#define input file as variable
pth_inpf03 <- fnm_towrite2
#Test the 'apex' package :
#read FASTA as multi.dna
mltdn_pip03 <- apex::read.multiFASTA(pth_inpf03)
#and read FASTA as DNAbin
dnb_pip03 <- ape::read.dna(pth_inpf03, format = "fasta")

# see the trimmed alignment
plot(mltdn_pip03, cex = 0.2)
# count the number sequences in the alignment
no_of_seq_pip3 <- mltdn_pip03@n.ind
# get the alignment
algn_pip03 <- apex::multidna2alignment(mltdn_pip03)
# count the characters in the sequence string
seq_lng_pip3 <- nchar(algn_pip03$seq[1])
df_apip3 <- as.data.frame(algn_pip03$seq)
str1_apip3 <- paste (df_apip3,collapse = "")
nofchars_apip3 <- nchar(str1_apip3)
no_char_in_mtrx_apip3 <- seq_lng_pip3*no_of_seq_pip3
indel_cnt_apip3 <- str_count(str1_apip3, pattern = "-")
perc_indel_apip3 <- indel_cnt_apip3/no_char_in_mtrx_apip3*100
# see the names of the sequences
#mltdn_pip03@labels
# rowbind the nested lists to a dataframe
df_pip03 <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(mltdn_pip03@labels),
                                  "_")))
#check object
colnames(df_pip03) <- c("smplNm","locality","year")
#df_pip03$locality3 <- paste0(df_pip03$locality,"_",df_pip03$locality2)
#make the multidna object a genind object
gi_pip03 <- apex::multidna2genind(mltdn_pip03)
#Make the genind object a data frame
df_geipip03 <- genind2df(gi_pip03)
#replace all NAs with "-"
df_geipip04 <- df_geipip03 %>% replace(is.na(.), "-")
df_geipip03 <- df_geipip04
#see column names of the data frame
#make genind and hierfstat objects by region
gi_pip_local <- adegenet::DNAbin2genind(dnb_pip03,pop=df_pip03$locality)
gi_pip_years <- adegenet::DNAbin2genind(dnb_pip03,pop=df_pip03$year) 
#make hierfstat objects - Note that it is required to set 'pop' to the 
# vector that matches the population to which the indivudal belongs
hfst_pip_lo <- hierfstat::genind2hierfstat(gi_pip_local,pop=df_pip03$locality)
hfst_pip_ye <- hierfstat::genind2hierfstat(gi_pip_years,pop=df_pip03$year)
hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
#replace all NAs with 0
hfst_pip_lo <- hfst_pip_lo %>% replace(is.na(.), 0)
hfst_pip_ye2 <- hfst_pip_ye %>% replace(is.na(.), 0)
#overwrite the previous dataframes
hfst_pip_lo3 <- hfst_pip_lo
hfst_pip_ye <- hfst_pip_ye2
#hfst_pip_ym <- hfst_pip_ym2
#pw_Nei_pip_ym <- hierfstat::pairwise.neifst(hfst_pip_ym, diploid = FALSE)
# # Replace NAs with zeros 
# pw_Nei_pip_un2 <- pw_Nei_pip_un %>% replace(is.na(.), 0)
# pw_Nei_pip_lo <- pw_Nei_pip_lo %>% replace(is.na(.), 0)
# pw_Nei_pip_ye <- pw_Nei_pip_ye %>% replace(is.na(.), 0)
# pw_Nei_pip_ym2 <- pw_Nei_pip_ym %>% replace(is.na(.), 0)
# #overwrite the previous dataframes
# pw_Nei_pip_un <- pw_Nei_pip_un2
# pw_Nei_pip_lo <- pw_Nei_pip_lo2
# pw_Nei_pip_ye <- pw_Nei_pip_ye2
# pw_Nei_pip_ym <- pw_Nei_pip_ym2

# Remove rows where the pop name is unique - 
# as the 'hierfstat::pairwise.WCfst'
# function cannot handle the unique pop rows only  represented by only a single individual
hfst_pip_lo <- hfst_pip_lo[hfst_pip_lo$pop %in% hfst_pip_lo$pop[duplicated(hfst_pip_lo$pop)],]
hfst_pip_ye2 <- hfst_pip_ye[hfst_pip_ye$pop %in% hfst_pip_ye$pop[duplicated(hfst_pip_ye$pop)],]
hfst_pip_lo <- hfst_pip_lo[order(hfst_pip_lo$pop),]
# plyr::count(hfst_pip_lo$pop)
# plyr::count(hfst_pip_ye2$pop)

#pairwise fst test Nei - this takes an 'hierfstat' object as input
#pw_Nei_pip_un <- hierfstat::pairwise.neifst(hfst_pip_un, diploid = FALSE)
pw_Nei_pip_lo <- hierfstat::pairwise.neifst(hfst_pip_lo, diploid = FALSE)
pw_Nei_pip_ye <- hierfstat::pairwise.neifst(hfst_pip_ye, diploid = FALSE)

#overwrite previous data frames
# The Pairwise fst test used in the heatmap colored table
# below makes use of the 'hierfstat' prepared object

#Pairwise fst test WC - this takes an 'hierfstat' object as input
 # pw_wc_pip_lo <- hierfstat::pairwise.WCfst(hfst_pip_lo, diploid = FALSE)
 # pw_wc_pip_ye <- hierfstat::pairwise.WCfst(hfst_pip_ye, diploid = FALSE)
#color the df with a heat map
#https://stackoverflow.com/questions/47733031/how-to-plot-dataframe-in-r-as-a-heatmap-grid
library(tidyverse)
# try it out on the pw_wc_pip_rg object prepared with
# the 'hierfstat' package
# make the data frame as a heat map
df_pw_pip <- as.data.frame(pw_Nei_pip_ye)
df_pw_pip2 <- df_pw_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#plot data frame as heat map
ggplot(df_pw_pip2, aes(x = rowname,
                          y = colname, 
                          fill = value)) +
  geom_tile()

#_________________________________________________________________
# make the data frame as a heat map
df_pw_pip <- as.data.frame(pw_Nei_pip_lo)
df_pw_pip2 <- df_pw_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#head(df_pw_pip2)
#plot data frame as heat map
ggplot(df_pw_pip2, aes(x = rowname,
                          y = colname, 
                          fill = value)) +
  geom_tile()
#_________________________________________________________________
# make the data frame as a heat map
df_pw_pip <- as.data.frame(pw_Nei_pip_ye )
df_pw_pip2 <- df_pw_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#head(df_pw_pip2)
#plot data frame as heat map
ggplot(df_pw_pip2, aes(x = rowname,
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
#plot data frame as heat map
ggplot(df_Nei_pip3, aes(x = rowname,
                        y = colname, 
                        fill = value)) +
  geom_tile()
#_________________________________________________________________
# make the data frame as a heat map
df_Nei_pip <- as.data.frame(pw_Nei_pip_lo)
df_Nei_pip3 <- df_Nei_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
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
#pw_Gst_pip_ym <- mmod::pairwise_Gst_Nei(gi_pip_yearm, linearized = FALSE)
pw_Gst_pip_lo <- mmod::pairwise_Gst_Nei(gi_pip_local, linearized = FALSE)
pw_Gst_pip_ye <- mmod::pairwise_Gst_Nei(gi_pip_years, linearized = FALSE)
# The Pairwise fst test used in the heatmap colored table
# below makes use of the 'hierfstat' prepared object
#Pairwise fst test WC - this takes an 'hierfstat' object as input
#pw_wc_pip_lo <- hierfstat::pairwise.WCfst(hfst_pip_lo3, diploid = FALSE)
#pw_wc_pip_ym <- hierfstat::pairwise.WCfst(hfst_pip_ym, diploid = FALSE)
#pw_wc_pip_ye <- hierfstat::pairwise.WCfst(hfst_pip_ye, diploid = FALSE)

# make the data frame as a heat map
df_pw_pip_lo <- as.data.frame(pw_Nei_pip_lo)
df_pw_pip_ye <- as.data.frame(pw_Nei_pip_ye)
#df_pw_pip_ye <- as.data.frame(pw_pip_ye)
#replace NAs with zeroes
df_pw_pip_lo[is.na(df_pw_pip_lo)] = 0
df_pw_pip_ye[is.na(df_pw_pip_ye)] = 0
#df_pw_pip_ye[is.na(df_pw_pip_ye)] = 0

#_____________________________________________________________
# Prepare min and max values for plot
#_____________________________________________________________
df_pw_pip <- df_pw_pip_ye
#df_pw_pip <- df_pw_pip_ye
#get the dplyr library to summarize
library(dplyr)
#get the max value per col
mx_pip <- df_pw_pip %>% summarise_if(is.numeric, max)
#get the min value per col
mn_pip <- df_pw_pip %>% summarise_if(is.numeric, min)
# get the max value among the max values obtained
mxx_pip <- max(mx_pip)
mnn_pip <- max(mn_pip)
#define upper and lower limits
upppl1 <- round(mxx_pip, 0)+1
lowmi1 <- round(mnn_pip, 0)-1
#rearrange dataframe to long 
df_pw_pip_2 <- df_pw_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
# see first header lines

#make a data frame for abbreviated sampling locations
locnm <- c("FynBogense","FynKerteminde","JyllandMariagerfjord",
           "NGermanyKielFjord","NGermanyMecklenburgerBuchtWismarBucht",
           "NJyllandLimfjord","NJyllandLimfjordLogstoer",
           "NWGermanyNSeaHelgolandRds","NWGermanyWaddenSeaBussumHaupstr",
           "SamsoeBallen","SEDenmarkMecklenburgerBucht","SjaellandSkovshoved")
locabb <- c("FBo","FKe","JMa","GKi","GWi","JLi","JLo","GHe","GBu","SBa","GWi","SSk")
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
#get the library 'pals' to be able to make use of colour gradients
library("pals")
#replace column names
colnames(df_pw_pip_2) <- c("rwnm_spcs","clnm_spcs","value")
#replace with abbreviated location names
locnm2 <- df_abloc$locabb[match(df_pw_pip_2$clnm_spcs, df_abloc$locnm)]
# assign back NCBI codes
locnm2[is.na(locnm2)]  <- df_pw_pip_2$clnm_spcs[is.na(match(df_pw_pip_2$clnm_spcs, df_abloc$locnm))]
#replace 'clnm_spcs' column
df_pw_pip_2$clnm_spcs <- locnm2

#prepare plot
h<-ggplot(df_pw_pip_2, aes(x = rwnm_spcs, 
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
df_clo03$locality13 <- df_clo03$locality6
df_clo03$locality13 <- gsub("SamsoeBallen","Ballen",df_clo03$locality13)
df_clo03$locality13 <- gsub("FynBogense","Bogense",df_clo03$locality13)
df_clo03$locality13 <- gsub("FynKerteminde","Kerteminde",df_clo03$locality13)
df_clo03$locality13 <- gsub("GermanyKielFjord","KielFjord",df_clo03$locality13)
df_clo03$locality13 <- gsub("NJyllandLimfjord","Loegstoer",df_clo03$locality13)
df_clo03$locality13 <- gsub("JyllandMariagerfjord","Mariagerfjord",df_clo03$locality13)
df_clo03$locality13 <- gsub("NGermanyMecklenburgerBuchtWismarBucht","MecklenburgerBucht",df_clo03$locality13)
df_clo03$locality13 <- gsub("Germany:Helgoland","NSeaHelgolandRds",df_clo03$locality13)
df_clo03$locality13 <- gsub("SjaellandSkovshoved","Skovshoved",df_clo03$locality13)
df_clo03$locality13 <- gsub("GermanyBusum","WaddenSeaBussumHaupstr",df_clo03$locality13)
# replace in column names and in row names
# colnames(df_pw_pip_lo) <- df_clo03$locality8[match(colnames(df_pw_pip_lo),df_clo03$locality13)]
# row.names(df_pw_pip_lo) <- df_clo03$locality8[match(row.names(df_pw_pip_lo),df_clo03$locality13)]

#df_pw_pip <- pw_Gst_pip_lo3
df_pw_pip <- df_pw_pip_lo
#get the dplyr library to summarize
library(dplyr)
#get the max value per col
mx_pip <- df_pw_pip %>% summarise_if(is.numeric, max)
#get the min value per col
mn_pip <- df_pw_pip %>% summarise_if(is.numeric, min)
# get the max value among the max values obtained
mxx_pip <- max(mx_pip)
mnn_pip <- max(mn_pip)
#define upper and lower limits
upppl1 <- round(mxx_pip, 0)+1
lowmi1 <- round(mnn_pip, 0)-1

#rearrange dataframe to long 
df_pw_pip_2 <- df_pw_pip %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
# see first header lines
#head(df_pw_pip_2)
#get the library 'pals' to be able to make use of colour gradients
library("pals")
#replace column names
colnames(df_pw_pip_2) <- c("rwnm_spcs","clnm_spcs","value")

#replace with abbreviated location names
locnm2 <- df_abloc$locabb[match(df_pw_pip_2$clnm_spcs, df_abloc$locnm)]
# assign back NCBI codes
locnm2[is.na(locnm2)]  <- df_pw_pip_2$clnm_spcs[is.na(match(df_pw_pip_2$clnm_spcs, df_abloc$locnm))]
#replace 'clnm_spcs' column
df_pw_pip_2$clnm_spcs <- locnm2
#replace with abbreviated location names
locnm2 <- df_abloc$locabb[match(df_pw_pip_2$rwnm_spcs, df_abloc$locnm)]
# assign back NCBI codes
locnm2[is.na(locnm2)]  <- df_pw_pip_2$rwnm_spcs[is.na(match(df_pw_pip_2$rwnm_spcs, df_abloc$locnm))]
#replace 'clnm_spcs' column
df_pw_pip_2$rwnm_spcs <- locnm2


#prepare plot
h2<-ggplot(df_pw_pip_2, aes(x = rwnm_spcs, 
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

#read FASTA as dna.bin
pip3 <- ape::read.dna(pth_inpf03, format = "fasta")
#create haplotypes from dna.bin
pipHaps3 <- pegas::haplotype(pip3)
#prepare hpt table
ind.hap3<-with(
  utils::stack(setNames(attr(pipHaps3, "index"), rownames(pipHaps3))),
  table(hap=ind, pop=rownames(pip3)[values]))
#make it a dataframe
df_ihpt01.3 <- as.data.frame(ind.hap3)
#limit to include only 'Freq' that equals 1
df_ihpt02.3 <- df_ihpt01.3[df_ihpt01.3$Freq == 1,]

liht2Nm <- strsplit(as.character(df_ihpt02.3$pop), "_")
locNm3 <- sapply(liht2Nm, "[[", 2)
#make it a table
hsll3 <- table(df_ihpt02$hap, locNm3)
ht3 <- pegas::haplotype(pip3)

tbl_hap_loc01.3 <- t(hsll3)
#make it a data frame instead
df_hap_loc01.3 <- as.data.frame(tbl_hap_loc01.3)
# reshape the data frame for long to wide
df_hap_loc02.3 <- reshape(data=df_hap_loc01.3,idvar="locNm3",
                          v.names = "Freq",
                          timevar = "Var2",
                          direction="wide")
# substitute to get matching location names
df_clo03$locality4 <- gsub("BalticSea","Baltic Sea",df_clo03$locality4)
df_clo03$locality4 <- gsub("CaspianSea","Caspian Sea",df_clo03$locality4)
df_clo03$locality4 <- gsub("CentralWAtlantic","CW Atlantic",df_clo03$locality4)
df_clo03$locality4 <- gsub("Germany:KielFjord","Germany, Kiel Fjord" ,df_clo03$locality4)
df_clo03$locality4 <- gsub("NEAtlantic" ,"NE Atlantic"  ,df_clo03$locality4)
df_clo03$locality4 <- gsub("AtlanticOcean:NWAtlantic"   ,"NW Atlantic"    ,df_clo03$locality4)

abbr8 <- c("B","CS","CWA","FBo","FKe","GBu","GHe","GHe","GKi","GKi","GKi","GWi","JLi","JMa","M","NEA","Nt","NWA","SBa","SSk","USA:G","USA:P","USA:W")
lNm8 <- c("BalticSea","CaspianSea","CentralWAtlantic","Bogense","Kerteminde","WaddenSeaBussumHaupstr","Germany:Helgoland","NSeaHelgolandRds","Germany:KielFjord","Germany:Maasholm","KielFjord","MecklenburgerBucht","Loegstoer","Mariagerfjord","Mediterranean","NEAtlantic","Netherlands","AtlanticOcean:NWAtlantic","Ballen","Skovshoved","USA:GalvestonBay","USA:Panacea","USA:WoodsHole")
df_lNm8<- as.data.frame(cbind(abbr8,lNm8))
df_clo03$locality9 <- df_lNm8$lNm8[match(df_clo03$locality7,df_lNm8$abbr8)]

# sum for duplicated values in row
# this is to add up the multiple entries row for the same 
# localities
# https://stackoverflow.com/questions/10180132/consolidate-duplicate-rows
library(plyr)
df_hap_loc03 <- plyr::ddply(df_hap_loc02.3,"locNm3",numcolwise(sum))
df_hap_loc03 <- df_hap_loc03[!is.na(df_hap_loc03[,1]),]
df_hap_loc02 <-   df_hap_loc03

#match between data frames
df_hap_loc02$dec_lat <- df_clo03$dec_lat2[match(df_hap_loc02$locNm3,df_clo03$locality2)]
df_hap_loc02$dec_lon <- df_clo03$dec_lon2[match(df_hap_loc02$locNm3,df_clo03$locality2)]
# limit the data frame to remove any rows that have NAs
df_hap_loc03 <-  df_hap_loc02[complete.cases(df_hap_loc02), ] 
#modify the colnames
colnames(df_hap_loc03) <- c(gsub("Freq\\.","",colnames(df_hap_loc03)))
# https://towardsdatascience.com/using-ggplot-to-plot-pie-charts-on-a-geographical-map-bb54d22d6e13
#count the number of columns, and subtract 2
enc1 <- ncol(df_hap_loc03)-2 
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
cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
colfunc <- colorRampPalette(cbbPalette2)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
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

#df_hap_loc03$rws2 <- df_hap_loc03$rws
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
                                  r = sqrt(rws/pi)*2.10), 
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
  
geom_scatterpie_legend(sqrt(df_hap_loc04$rws/pi)*2.10, x=-40, y=30) +
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
#p01

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
                                  r = sqrt(rws/pi)*0.24), 
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
geom_scatterpie_legend(sqrt(df_hap_loc04$rws/pi)*0.24, x=14, y=57) +
  
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
#install.packages("rnaturalearthdata", repos = "http://packages.ropensci.org", type = "source")
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

jitlvl <- 0.24
jitlvl <- 0.024
dev.off()
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
# cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
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
# add color range as a column
df_clo2$colfloc <- cl
# write out the table with colours
write.csv(df_clo2,file=paste0(wd00_wd05,"/df_clo2.csv"))

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
p04 <- p04 + theme(panel.background = element_rect(fill = 'white', color = 'white'),
                   panel.grid.major = element_line(color = 'white')) #, linetype = 'dotted'))#,
# add border around plot
p04 <- p04 + theme(panel.border = element_rect(color = "black",
                                               fill = NA,
                                               size = 1.0))
# change background of legend
p04 <- p04 + theme(legend.key = element_rect(fill = "white"))
# see the plot
#p04
#make filename to save plot to
figname06 <- paste0("map_samples_",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname06,sep="")
if(bSaveFigures==T){
  ggsave(p04,file=figname02,width=210,height=297*0.5,
         units="mm",dpi=300)
}


#make filename to save plot to
figname06 <- paste0("Fig01_map_samples_",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname06,sep="")
if(bSaveFigures==T){
  ggsave(p04,file=figname02,width=210,height=297*0.5,
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
                                  r = sqrt(rws/pi)*0.60), 
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
  
  geom_scatterpie_legend(sqrt(df_hap_loc04$rws/pi)*0.60, x=-10, y=47) +
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
#p05
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
                                  r = sqrt(rws/pi)*0.24), 
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
  geom_scatterpie_legend(sqrt(df_hap_loc04$rws/pi)*0.24, x=14, y=57) +
  
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  # set alpha values for color intensity of fill color in point
  #https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html
  #define limits of the plot
  ggplot2::coord_sf(xlim = c(4, 17.4),
                    ylim = c(53.8, 58.0),
                    expand = FALSE) +
  theme(aspect.ratio=2/7)
# see the plot
#p06
#dev.off()
#change labels on axis
p06 <- p06 + xlab("Longitude") + ylab("Latitude")
#change labels on legend
p06 <- p06 + labs(fill='haplotype')

#https://www.statology.org/ggplot-background-color/
p06 <- p06 + theme(panel.background = element_rect(fill = 'white', color = 'white')) #,
# panel.grid.major = element_line(color = 'red' , linetype = 'dotted'),
# panel.grid.minor = element_line(color = 'green', size = 2))
#p06
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

# copy column
df_clo03$locality17 <-  df_clo03$locality4
# substitute in column
df_clo03$locality17 <- gsub("Funen, Bogense","Funen, Bogense,",df_clo03$locality17)
# match back to get abbreviated location name
df_clo2$abbl<- df_clo03$locality8[match(df_clo2$locality,df_clo03$locality17)]
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
#cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
# This color range uses very dark red colors, and these might not
# be easy to see for someone that is red-green colour blind
cbbPalette1 <- c("firebrick4","firebrick2",#"orange",
                 "gold3")
# As an alternative use this palette of grey hues
cbbPalette1 <- c("gray40","gray65",#"orange",
                 "gray80", "cyan1","cyan3","cyan4","deepskyblue4")
#read FASTA as dna.bin
pip3 <- ape::read.dna(pth_inpf03, format = "fasta")
lNm4 <- strsplit(as.character(row.names(pip3)), "_")
lNm1 <- sapply(lNm4, "[[", 1)
lNm2 <- sapply(lNm4, "[[", 2)
rLnm <- unique(lNm2[grepl("[A-Z]{2}[0-9]{5}",lNm1)])
rLnm <- rLnm[order(rLnm)]
lrLnm <- length(rLnm)
colfunc2 <- colorRampPalette(cbbPalette2)
colfunc1 <- colorRampPalette(cbbPalette1)
cl <- cbbPalette
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
#for NCBI categories
colfpl1 <- colfunc1(lrLnm)
# bind to dataframes
df_Cl1 <- as.data.frame(cbind(rLnm,colfpl1))
df_Cl2 <- as.data.frame(cbind(df_clo2$abbl,df_clo2$colfloc))
# reorder the collection location names 
df_Cl2 <- df_Cl2[order(df_Cl2[,1]),]
# rename the column headers
colnames(df_Cl1) <- c("coll_loc","colfcol_loc")
colnames(df_Cl2) <- c("coll_loc","colfcol_loc")
# limit data frame to only include non matching remote locations
# to avoid having mulitple colors for locations that appear in both your own samples
# and the NCBI obtained samples
df_Cl1 <- df_Cl1[is.na(!match(df_Cl1$coll_loc,df_Cl2$coll_loc)),]
#bind the data frames together in to one data frame
df_cll <- rbind(df_Cl1,df_Cl2)
#copy DNAbin object again
dnb_pip4 <- pip3
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
  utils::stack(setNames(attr(ht4, "index"), rownames(ht4))),
  table(hap=ind, pop=rownames(dnb_pip4)[values]))
#make it a dataframe
df_ihpt04 <- as.data.frame(ind.hap4)
# make haplotype labels arabic numerals instead of roman numerals 
df_ihpt04$hap.ab <- as.numeric(as.roman(as.character(df_ihpt04$hap)))
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
# #make the plot
# plot(hN4, 
#      size = sqrt(attr(hN4,"freq")/pi), 
#      #size = log10(attr(hN4,"freq")), 
#      scale.ratio = 0.6, 
#      cex = 1.1, # set size of roman numerals on circles for haplotype ID
#      #bg= colfh, 
#      pie = hsl4, 
#      show.mutation = 2, threshold = 0, labels(T))
# 
# 
#
# write out table for colors to match catch locations
write.csv(df_cll, file=paste0(wd00_wd05,"/df_cll_colf_smpl_loc.csv"))
# make a color range
cbbPalette3 <- c("black","cyan","white")
cbbPalette3 <- c("black","firebrick2","brown","blue","cyan","white")
# make the color range a function
colfunc3 <- colorRampPalette(cbbPalette3)
# make a number of color steps based on the number unique sample years
clfy2 <- colfunc3(length(usmplyear3))
# reorder the list of sample years
smplyeord<- usmplyear3[order(usmplyear3)]
# bind together the list of years with the color steps
df_cfy3 <- as.data.frame(cbind(smplyeord,clfy2))
# Rename the column heasders in this data frame
colnames(df_cfy3) <- c("smplyr","colour")
# get colours to use for pies
colfh<- df_cll$colfcol_loc[match(colnames(hsl4),df_cll$coll_loc)]
colfh_y <- df_cfy3$colour[match(colnames(hsl5),df_cfy3$smplyr)]
# write the color table for sampling years
write.csv(df_cfy3,file=paste0(wd00_wd05,"/df_cfy3.csv"))
#_______________________________________________________________________________
#_______________________________________________________________________________
flnm <- c(paste("Fig03_v01_haplotype_network_",inp.f.fnm,"02.pdf",  sep = ""))
flnm <- c(paste("Fig03_v01_haplotype_network_",inp.f.fnm,"02.jpg",  sep = ""))
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
#make the plot
plot(hN4, 
     size = sqrt(attr(hN4,"freq")/pi), 
     #size = log10(attr(hN4,"freq")), 
     scale.ratio = 0.6, 
     cex = 1.1, # set size of roman numerals on circles for haplotype ID
     bg= colfh, 
     pie = hsl4, 
     show.mutation = 2, threshold = 0, labels(T))

#add a legend to the plot
legend("bottomleft",colnames(hsl4), 
       pt.bg=colfh,box.col=NA,
       #col=rainbow(ncol(new.hap.smplloc)), 
       pt.lwd=0.4,
       pch=21, ncol=1, cex=0.8)
title(main = "a",
      cex.main = 1.8,   font.main= 2, col.main= "black",
      adj = 0.01, line = 0.1)
#_______________________________________________________________________________
#make the plot
plot(hN4, 
     size = sqrt(attr(hN4,"freq")/pi), 
     #size = log10(attr(hN4,"freq")), 
     scale.ratio = 0.6, 
     cex = 1.1, # set size of roman numerals on circles for haplotype ID
     bg= colfh_y, 
     pie = hsl5, 
     show.mutation = 2, threshold = 0, labels(T))

#add a legend to the plot
legend("bottomleft",colnames(hsl5), 
       #col=rainbow(ncol(new.hap.smplye)), 
       pt.bg=colfh_y,box.col=NA,pt.lwd=0.4,
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
pth_inpf03 <- fnm_towrite2
# read in fasta
df_trimfas16 <- adegenet::genind2df(adegenet::DNAbin2genind(strataG::read.fasta(pth_inpf03 )))
# get the row names
rwnm01 <- rownames(df_trimfas16)
# split by delimeter and make a new data frame
df_r02 <- data.frame(do.call('rbind', strsplit(as.character(rwnm01),'_',fixed=TRUE)))
#change column names on data frame
colnames(df_r02) <- c("smplnm","collloc","collyear")
#replace in location names
df_clo2$locality <- gsub(",","",df_clo2$locality)

NCBIsmplN <-  df_r02$smplnm[grepl("[A-Z]{2}[0-9]{5}",df_r02$smplnm)]
ownsmplN  <-  df_r02$collloc[!grepl("[A-Z]{2}[0-9]{5}",df_r02$smplnm)]
#make empty columns
df_r02$declat2 <- NA
df_r02$declon2 <- NA
# match to get sampling latitude and longitude for NCBI samples
df_r02$declat2[grepl("[A-Z]{2}[0-9]{5}",df_r02$smplnm)] <- df_lN02$declat2[match(NCBIsmplN,df_lN02$accession_nmb)]
df_r02$declon2[grepl("[A-Z]{2}[0-9]{5}",df_r02$smplnm)] <- df_lN02$declon2[match(NCBIsmplN,df_lN02$accession_nmb)]
# match to get sampling latitude and longitude for own samples
df_r02$declat2[!grepl("[A-Z]{2}[0-9]{5}",df_r02$smplnm)] <- df_clo03$dec_lat2[match(ownsmplN,df_clo03$locality8)]
df_r02$declon2[!grepl("[A-Z]{2}[0-9]{5}",df_r02$smplnm)] <- df_clo03$dec_lon2[match(ownsmplN,df_clo03$locality8)]
#make numeric
df_r02$declat2 <- as.numeric(df_r02$declat2)
df_r02$declon2 <- as.numeric(df_r02$declon2)
# use the biogeo package to get degress minutes and seconds from decimal degrees
r02_declat2 <- biogeo::dd2dmslat(df_r02$declat2)
r02_declon2 <- biogeo::dd2dmslong(df_r02$declon2)

# get sampling position
df_r02$lat_deg <- r02_declat2$ydeg
df_r02$lat_min <- r02_declat2$ymin
df_r02$lat_sec <- r02_declat2$ysec
df_r02$lat_sph <- r02_declat2$NS
df_r02$lon_deg <- r02_declon2$xdeg
df_r02$lon_min <- r02_declon2$xmin
df_r02$lon_sec <- r02_declon2$xsec
df_r02$lon_sph <- r02_declon2$EW
# make an empty column for accession numbers
df_r02$NCBIaccsNo <- NA
#replace in this column with accession numbers
df_r02$NCBIaccsNo[grepl("[A-Z]{2}",df_r02$smplnm)] <- gsub("Mnelei","",df_r02$smplnm[grepl("[A-Z]{2}",df_r02$smplnm)])
# make an empty column for species names
df_r02$spcNm <- NA
df_r02$spcNm[grepl("Mnelei",df_r02$smplnm)] <- "Mnemiopsis leidyi"
df_r02$spcNm <- "Mnemiopsis leidyi"
#View(df_r02)
#replace the non Mnemiopsis leidyi samples
notMnelei <- df_r02$smplnm[!grepl("Mnelei",df_r02$smplnm)]
# if there are entries that are not Mnelei
# if(length(notMnelei)>1){
# notMnelei[grepl("^.*[A-Z]{2}.*$",notMnelei)] <- gsub("^(.*)([A-Z]{2}.*)$","\\1_\\2",notMnelei[grepl("^.*[A-Z]{2}.*$",notMnelei)])
# notMnelei[grepl("^.*\\.[A-Z]{1}.*$",notMnelei)] <- gsub("\\.","_",notMnelei[grepl("^.*\\.[A-Z]{1}.*$",notMnelei)])
# # split by delimeter and make a new data frame
# df_nMl <- data.frame(do.call('rbind', strsplit(as.character(notMnelei),'_',fixed=TRUE)))
# colnames(df_nMl) <- c("spcNm","NCBIaccsNos")
# #get non Mnemiopsis samples and modify
# df_r02$spcNm[!grepl("Mnelei",df_r02$smplnm)] <- df_nMl$spcNm
# df_r02$NCBIaccsNo[!grepl("Mnelei",df_r02$smplnm)] <- df_nMl$NCBIaccsNos}
# #replace accession numbers
# df_r02$NCBIaccsNosmplnm <- df_r02$NCBIaccsNo
# df_r02$NCBIaccsNosmplnm[is.na(df_r02$NCBIaccsNosmplnm)] <- df_r02$smplnm[is.na(df_r02$NCBIaccsNosmplnm)]
df_r02$NCBIaccsNosmplnm <- df_r02$smplnm
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
df_r02$collloc2 <- df_r02$collloc01
# round decimal degrees and paste together in a string
df_clo03$dec_lat2 <- round(as.numeric(df_clo03$dec_lat),3)
df_clo03$dec_lon2 <-  round(as.numeric(df_clo03$dec_lon),3)
# Keep 3 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
df_clo03$dec_lon2 <- sprintf("%.3f", df_clo03$dec_lon2)
df_clo03$dec_lat2 <- sprintf("%.3f", df_clo03$dec_lat2)
df_clo03$dec_latlon2 <- paste0(df_clo03$dec_lat2,"; ",df_clo03$dec_lon2)
#match to get sampling location 
df_r02$posloc <- df_clo03$dec_latlon2[match(df_r02$collloc,df_clo03$locality8)]
# get long name for sampling location
df_r02$lNmloc <- df_clo03$locality5[match(df_r02$collloc,df_clo03$locality8)]

# match with previously prepared data frame to get collectoin month
df_r02$collmnth <- df_ownsmpldt$smplmnt4[match(df_r02$smplnm,df_ownsmpldt$MnsmNo)]
#substitute to get accession number
df_r02$collyear2 <- paste0(df_r02$collyear,", ",df_r02$collmnth)
# substitute for the year and months pasted incorrectly
df_r02$collyear2 <- gsub(", NA","",df_r02$collyear2)
df_r02$collyear2 <- gsub("^([0-9]{4}), [0-9]{4}$","\\1",df_r02$collyear2)

# get the sample names that have Mnelei instead of the correct
# accession number
MneleiAbbrcodes <- df_r02$NCBIaccsNosmplnm[grepl("Mnelei",df_r02$NCBIaccsNosmplnm)]
# use match to get the accession numbers obtained from
# submitting the sequences to NCBI GenBank
Mnelei.Acc.Nos <- df_accNos$AccNo[match(MneleiAbbrcodes,df_accNos$MneleiAbbr)]
# add back the new accession numbers for the sequences deposited
df_r02$NCBIaccsNosmplnm[grepl("Mnelei",df_r02$NCBIaccsNosmplnm)] <- Mnelei.Acc.Nos
#df_r02$NCBIaccsNosmplnm


# define columns to keep
keeps <- c("spcNm",
           "NCBIaccsNosmplnm",
           "posloc",
           "collyear2",
           "collloc" )
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
lcnm8 <- gsub("Atlantic Ocean, ","",lcnm8)
abbr8 <- gsub("Ssk","SSk",abbr8)
#paste long nmae and abbreviatoin together to use for table header
ablo8 <- paste0(lcnm8," (",abbr8,")")
# exclude matches
ablo8<- ablo8[!grepl("GalvestonBay",ablo8)]
ablo8<- ablo8[!grepl("Panacea",ablo8)]
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
df_hap_loc01 <- as.data.frame(tbl_hap_loc01.3)
# reshape the data frame for long to wide
df_hap_loc02 <- reshape(data=df_hap_loc01,idvar="locNm3",
                        v.names = "Freq",
                        timevar = "Var2",
                        direction="wide")

df_hap_loc02$abNm2 <-   df_hap_loc02$locNm3  

#add a column with sampling locations to be able to 
#match between data frames
#df_hap_loc02$dec_loc3 <- df_clo$locality[match(df_hap_loc02$smplloca,df_clo2$locality)]
df_hap_loc02$dec_loc3 <- df_hap_loc02$abNm2

# sum for duplicated values in row
# this is to add up the multiple entries row for the same 
# localities
# https://stackoverflow.com/questions/10180132/consolidate-duplicate-rows
library(plyr)
df_hap_loc03 <- plyr::ddply(df_hap_loc02,"dec_loc3",numcolwise(sum))
#match between data frames
df_hap_loc03$dec_lat <- df_clo03$dec_lat2[match(df_hap_loc03$dec_loc3,df_clo03$locality8)]
df_hap_loc03$dec_lon <- df_clo03$dec_lon2[match(df_hap_loc03$dec_loc3,df_clo03$locality8)]
# make latitude longitude numeric
df_hap_loc03$dec_lat <- as.numeric(df_hap_loc03$dec_lat)
df_hap_loc03$dec_lon <- as.numeric(df_hap_loc03$dec_lon)
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


####
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
colfunc2 <- colorRampPalette(cbbPalette2)
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

#make a viridis colour range
cl03 <- pals::viridis(length(unique(df_hap_loc03[,c(2:enc)])))
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
                                  r = sqrt(rws*0.08/pi)), 
                              data = df_hap_loc04, 
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  scale_color_manual(values=c(rep("black",
                                  length(unique(df_hap_loc04[,c(2:enc)]))))) +
  scale_fill_manual(values=alpha(
    c(cl03),
    c(0.7)
  ))+
  #https://stackoverflow.com/questions/54078772/ggplot-scale-color-manual-with-breaks-does-not-match-expected-order
  geom_scatterpie_legend(sqrt((df_hap_loc04$rws*0.08)/pi), x=-10, y=47) +
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
#p07
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
#__________________
#
mtr_dd_pip <- ape::dist.dna(pip3)
mtr_dd_pip[is.na(mtr_dd_pip)] <- 0
tre_pip <- ape::njs(mtr_dd_pip)
#sort the branches in the tree
tre_pipr <- ape::ladderize(tre_pip, right = TRUE)

#https://joey711.github.io/phyloseq/plot_tree-examples.html
plot(tre_pipr, cex=0.4)
# Plot a neighbour join tree
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
# make the plot object
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
#p
# make a data frame of the tip label names
df_tiplb01 <- as.data.frame(cbind(c(tre_pipr$tip.label)))
sAb <- gsub("^(.*)_(.*)_(.*)$","\\2",df_tiplb01$V1)

df_tiplb01$seqNm <- df_olnwNm$smplolNm[match(sAb,df_olnwNm$sAbbrv)]
df_tiplb01$cat <- NA
colnames(df_tiplb01) <- c("newNm","seqNm", "cat")

df_tiplb01$cat[grepl("NCBI",df_tiplb01$seqNm)] <- "NCBI"
df_tiplb01$cat[grepl("[A-Z]{2}[0-9]{5}",df_tiplb01$seqNm)] <- "NCBI"
df_tiplb01$cat[grepl("Germany",df_tiplb01$seqNm)] <- "Germany"
df_tiplb01$cat[grepl("Jylland",df_tiplb01$seqNm)] <- "Jylland"
df_tiplb01$cat[grepl("_Fyn",df_tiplb01$seqNm)] <- "Fyn"
df_tiplb01$cat[grepl("Samsoe",df_tiplb01$seqNm)] <- "Samsoe"
df_tiplb01$cat[grepl("Sjaelland",df_tiplb01$seqNm)] <- "Sjaelland"
df_tiplb01$cat[grepl("Limfjord",df_tiplb01$seqNm)] <- "Jylland"
df_tiplb01$cat[grepl("Mecklenburger",df_tiplb01$seqNm)] <- "Germany"
tipcategories <- df_tiplb01
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")


#cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
colfunc <- colorRampPalette(cbbPalette2)
# copy the palette
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
df_clo2$locality <- gsub("NWGermanyNSeaHelgolandRds", "North Sea, Helgoland Roads", df_clo2$locality )
df_clo2$locality <- gsub("NWGermanyWaddenSeaBussumHaupstr","Germany, Büsum", df_clo2$locality )
df_clo2$locality <- gsub("SamsoeBallen","Samsøe, Ballen", df_clo2$locality )
df_clo2$locality <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", df_clo2$locality )
#identify unique localities
nloc <- length(unique(df_clo2$locality))
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# use the color function to make a color range across the number of locations
cl <- colfunc(nloc)
#make a dataframe that comprises colors and location names
dfcol01 <- as.data.frame(cbind(cl,unique(df_clo2$locality)))
colnames(dfcol01) <- c("col","loc")
dd = as.data.frame(tipcategories)
#Add another row to accomodate the sequence that only has 'Mecklenburger' in the sequence name
df_Meckl <- as.data.frame(c(dfcol01[grepl("Wismar",dfcol01$loc),][1],"Mecklenburger"))
colnames(df_Meckl) <- colnames(dfcol01)
dfcol01 <- rbind(dfcol01,df_Meckl)
df_clo03$locality19 <- NULL
df_clo03$locality19 <- df_clo03$locality17
ncm1 <- ncol(df_clo03)
df_tmp <- cbind(df_clo03[df_clo03$locality19=="Germany, Wismar Bight",][1:ncm1-1],"Mecklenburger")
colnames(df_tmp) <- colnames(df_clo03) 
df_clo04 <- rbind(df_tmp,df_clo03)
df_clo03 <- df_clo04
dfcol01$Abrloc <- df_clo03$locality8[match(dfcol01$loc,df_clo03$locality19)]
dfcol01$col2 <- df_cll$colfcol_loc[match(dfcol01$Abrloc,df_cll$coll_loc)]

#substitute to get second part of location name
loc2 <- gsub("(.*), (.*)","\\2",dfcol01$loc)
#substitute 
loc2 <- gsub(",","",loc2)
loc2 <- gsub(" ","",loc2)
loc2 <- gsub("Roads","Rds",loc2)
loc2 <- gsub("WismarBight","Wismar",loc2)
loc2 <- gsub("MariagerFjord","Mariagerfjord",loc2)
loc2 <- gsub("Büsum","Bussum",loc2)
# make new columns that are empty , you will fill them below
dd$cat2 <- NA
dd$loc2 <- NA
#count the number of rows
ncol2 <- nrow(dfcol01)
#for element in number rows, match the location, and assign the hex color
# to the new empty column
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
  "Netherlands")
# count the number of elements
noNCBIsmpl <- length(NCBIsmpl)
NCBIsmpl2<- gsub("Mnemiopsis_leidyi_","",dd$seqNm[dd$cat=="NCBI"])
dd$loc2[dd$cat=="NCBI"] <- df_lN02$location[match(NCBIsmpl2,df_lN02$accession_nmb)]
dd$loc2 <- gsub("Germany:Helgoland","NWGermanyNSeaHelgolandRds",dd$loc2)
dd$loc2 <- gsub("Germany:KielFjord","Germany:Maasholm",dd$loc2)
dd$loc2 <- gsub("USA:GalvestonBay","USA:WoodsHole",dd$loc2)
dd$loc2 <- gsub("USA:Panacea","USA:WoodsHole",dd$loc2)
#match to get location abbreviation on the df_lN02 data frame
df_lN02$loc.abbrv <- df_pip03$locality[match(df_lN02$accession_nmb,df_pip03$smplNm)]
# write the df_lN02 data frame to a csv file
write.csv(df_lN02,paste0(wd00_wd05,"/df_lN02.csv"))


# make a color range
cbbPalette3 <- c("brown4","brown3","brown2")
cbbPalette3 <- c("azure4","azure3","azure2")
cbbPalette3 <- c("cadetblue4","cadetblue3","cadetblue2")
cbbPalette3 <- cbbPalette1
colfunc2 <- colorRampPalette(cbbPalette3)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# Use the color function the make a color range across the number of samples
clNBCI <- colfunc2(noNCBIsmpl)
# and bind this into a data frame
dfNCBIcol <- as.data.frame(cbind(NCBIsmpl,clNBCI))
#count the number of rows
ncol3 <- nrow(dfNCBIcol)
#for element in number rows, match the location, and assign the hex color
# to the new column
for (e in seq(1:ncol3)){
  dd$col2[grepl(dfNCBIcol$NCBIsmpl[e],dd$loc2)] <- dfNCBIcol$clNBCI[e]
  dd$loc2[grepl(dfNCBIcol$NCBIsmpl[e],dd$loc2)] <- dfNCBIcol$NCBIsmpl[e]
}
#subset dataframe to only comprise column 4 and 5
dd2 <- dd[,c("loc2","col2")]
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
dd$loc3 <- dd$loc2
#Begin plotting the tree - assign the plot to an object
p01 <- p %<+% dd + 
  geom_tiplab(aes(fill =loc3, 
                  color=loc3
  ),
  #color = "black", # color for label font # if this is active, you cannot change the color of the the text manually in the section below
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
nct <- length(unique(factor(dd$loc3)))
# Make another color gradient
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
#cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
# make the color palette a ramp of colors
colfunc <- colorRampPalette(cbbPalette2)
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
cl <- colfunc(nct)
# Use transposed data frame with unique colors
dd4 <- dd2
dd4[1,] <- "black"
dd4[][grepl("Bogense",colnames(dd4))] <- "white"
dd4[][grepl("Kerteminde",colnames(dd4))] <- "white"
dd4[][grepl("Mariager",colnames(dd4))] <- "white"
dd4[][grepl("Kiel Fjord",colnames(dd4))] <- "white"
vl<- as.character(dd4[1,])
p01 <- p01 + scale_colour_manual(values=c(vl))
p01 <- p01 + scale_fill_manual(name = "col2", values = alpha(c(dd2),c(0.7) ))  
p01 <- p01 + theme(legend.position = "none")
# See the plot object
#p01
# make a data frame with colors
df_dd02 <- as.data.frame(t(as.data.frame(dd2)))
# assign row names to a column
df_dd02$smpl_loc <- rownames(df_dd02)
# write the csvdata frame as a csv file 
write.csv(df_dd02,file=paste0(wd00_wd05,"/df_dd02.csv"))
#dev.off()
#make filename to save plot to
figname01 <- paste0("Fig02_v02_NJtree_",inpf01,".png")

figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(p01,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}

# Make a color palette
cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
                 "yellow","white")
#cbbPalette2 <- c("white","yellow","orange","tomato","red","brown","black")
# make a color ramp from this palette
colfunc <- colorRampPalette(cbbPalette2)
#identify unique localities
nloc <- length(unique(df_clo2$locality))
#https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# use the color function to make an equal number of colors for the sampling localities
cl <- colfunc(nloc)
# ensure count is a numeric value
df_hap_loc03$rws <- as.numeric(df_hap_loc03$rws)
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p04 <- 
  ggplot(data = denm_map) +
  #ggplot(st_transform(norw_map, 9122)) +
  geom_sf(color = "black", fill = "azure3") +
  theme(aspect.ratio=3/7) +
  scatterpie::geom_scatterpie(data = df_hap_loc03,
                              aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  # make the area of the circle equal the sample size, by taking the square root to the count as the radius divided by pi
                                  r = sqrt(rws/pi)*0.3),
                              
                              cols = colnames(df_hap_loc04[,c(2:enc)])) +
  geom_scatterpie_legend((sqrt(df_hap_loc03$rws/pi))*0.3, x=14, y=57) +
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
#p04
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
# grep for Mediterranean sea to check position
#df_clo03[grepl("Mediterr",df_clo03$locality5),]
#replace NAs with zeros
df_hap_loc06 <- df_hap_loc05[!is.na(df_hap_loc05$rws),]
df_hap_loc06$dec_lat <- as.numeric(df_hap_loc06$dec_lat)
df_hap_loc06$dec_lon <- as.numeric(df_hap_loc06$dec_lon)
#df_hap_loc06[grepl("Hl",df_hap_loc06$dec_loc2),]
#df_hap_loc06[grepl("Nt",df_hap_loc06$dec_loc2),]
#https://guangchuangyu.github.io/2016/12/scatterpie-for-plotting-pies-on-ggplot/
#world3 <- ggplot2::map_data('world3')
jitlvl <- 0.017
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p08 <- 
  ggplot(data = world3) +
  #ggplot(st_transform(norw_map, 9122)) +
  geom_sf(color = "black", fill = "azure3") +
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = sqrt(rws*0.30/pi)), 
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
  geom_scatterpie_legend(sqrt(df_hap_loc06$rws/pi), x=-10, y=47, 
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
#p08
# also see : https://github.com/tidyverse/ggplot2/issues/2037
p09 <- ggplot(data = denm_map) +
  #ggplot(st_transform(norw_map, 9122)) +
  geom_sf(color = "black", fill = "azure3") +
  scatterpie::geom_scatterpie(aes(x=dec_lon, y=dec_lat, 
                                  #group = country, 
                                  r = sqrt(rws/pi)*0.2), 
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
  geom_scatterpie_legend(sqrt(df_hap_loc04$rws/pi)*0.2, x=13.4, y=57.0,labeller = function(ra) ra * 1/0.05) +
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
#p09
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
# paste path and file together
loc06fl <- paste(wd00_wd05,"/loc06.csv",sep="")
# write csv
write_csv(df_hap_loc08,file=loc06fl)

#_______________________________________________________________________________
# match between abbreviated location names to get long location names
#View(df_clo03)
df_clo03$locality9 <- df_clo03$locality5
df_clo03$locality9 <- gsub("BalticSea" ,"Baltic Sea" ,df_clo03$locality9)
df_clo03$locality9 <- gsub("CaspianSea" ,"Caspian Sea" ,df_clo03$locality9)
df_clo03$locality9 <- gsub("Germany:Helgoland","North Sea, Helgoland Roads" ,df_clo03$locality9)
df_clo03$locality9 <- gsub("CentralWAtlantic","CW Atlantic"  ,df_clo03$locality9)
df_clo03$locality9 <- gsub("Germany:KielFjord" ,"Germany, Kiel Fjord"   ,df_clo03$locality9)
df_clo03$locality9 <- gsub("NEAtlantic" ,"NE Atlantic" ,df_clo03$locality9)
df_clo03$locality9 <- gsub("AtlanticOcean:NWAtlantic" ,"NW Atlantic" ,df_clo03$locality9)
#match to get long name for locations
df_hap_loc08$lnNm <- df_clo03$locality9[match(df_hap_loc08$dec_loc3,df_clo03$locality2 )]
# add a new colunm to add overall geographic regions to
df_hap_loc08$ov.aL <- df_hap_loc08$lnNm
df_hap_loc09 <- df_hap_loc08
#df_hap_loc06 <- df_hap_loc07[complete.cases(df_hap_loc07),]
# paste path and file together
loc06fl <- paste(wd00_wd05,"/loc06.csv",sep="")
# write csv
write_csv(df_hap_loc09,file=loc06fl)
# start adding overall geographical regions based on matches
df_hap_loc08$ov.aL[grepl("Baltic",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("Germany",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("Jutland",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("Funen",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("North Sea",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("NEAtlantic",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("NE Atlantic",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("Samsoe",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("Sealand",df_hap_loc08$lnNm)] <- "NEurope"
df_hap_loc08$ov.aL[grepl("Netherlands",df_hap_loc08$lnNm)] <- "NEurope"

df_hap_loc08$ov.aL[grepl("NW",df_hap_loc08$lnNm)] <- "NWAtlantic"
df_hap_loc08$ov.aL[grepl("USA",df_hap_loc08$lnNm)] <- "NWAtlantic"

# sum up Hpt counts from different sampling locations within same geographic region
df_hap_loc07 <- aggregate(df_hap_loc08[,sapply(df_hap_loc08,
                                               is.numeric)],df_hap_loc08["ov.aL"],sum)
#replace the lat lon positions as they have been summed up
df_hap_loc07$dec_lat <- df_hap_loc08$dec_lat[match(df_hap_loc07$ov.aL,df_hap_loc08$ov.aL)]
df_hap_loc07$dec_lon <- df_hap_loc08$dec_lon[match(df_hap_loc07$ov.aL,df_hap_loc08$ov.aL)]

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
                                  r = sqrt(rws/pi)), 
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
  geom_scatterpie_legend(sqrt(df_hap_loc06$rws/pi), x=-40, y=40, 
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


#Make the dnabin a genind object
gei_pip3 <- adegenet::DNAbin2genind(pip3)
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
#h.pca.scores$pop
#get pop groups from row names
popyear <- gsub("(^.*)_(.*)_(.*)$","\\3",row.names(h.pca.scores))
poploc <- gsub("(^.*)_(.*)_(.*)$","\\2",row.names(h.pca.scores))
# grep for accession numbers, and only replace matching locations in accession numbers
df_clo03$locality10 <- NULL
df_clo03$locality10 <- df_clo03$locality6
df_clo03$locality10 <- gsub("FynBogense","Bogense",df_clo03$locality10)
df_clo03$locality10 <- gsub("SamsoeBallen","Ballen",df_clo03$locality10)
df_clo03$locality10 <- gsub("FynKerteminde","Kerteminde",df_clo03$locality10)
df_clo03$locality10 <- gsub("NJyllandLimfjord","Limfjord",df_clo03$locality10)
df_clo03$locality10 <- gsub("JyllandMariagerfjord","Mariagerfjord",df_clo03$locality10)
df_clo03$locality10 <- gsub("SjaellandSkovshoved","Skovshoved",df_clo03$locality10)
ncldfcl03 <- ncol(df_clo03)
ncm1 <- ncldfcl03-1
df_tmp <- cbind(df_clo03[df_clo03$locality10=="Limfjord",][1:ncm1],"Loegstoer")
colnames(df_tmp) <- colnames(df_clo03) 
df_clo04 <- rbind(df_tmp,df_clo03)

df_tmp <- cbind(df_clo04[df_clo04$locality10=="GermanyKielFjord",][1:ncm1],"Germany:Maasholm")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)

df_tmp <- cbind(df_clo04[df_clo04$locality10=="GermanyKielFjord",][1:ncm1],"KielFjord")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)

df_tmp <- cbind(df_clo04[df_clo04$locality10=="GermanyKielFjord",][1:ncm1],"Germany:KielFjord")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)

df_tmp <- cbind(df_clo04[df_clo04$locality10=="NWGermanyNSeaHelgolandRds",][1:ncm1],"NSeaHelgolandRds")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)

df_tmp <- cbind(df_clo04[df_clo04$locality10=="NGermanyMecklenburgerBuchtWismarBucht",][1:ncm1],"MecklenburgerBuchtWismarBucht")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)

df_tmp <- cbind(df_clo04[df_clo04$locality10=="NGermanyMecklenburgerBuchtWismarBucht",][1:ncm1],"MecklenburgerBuchtWismarBucht")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)

df_tmp <- cbind(df_clo04[df_clo04$locality10=="NGermanyMecklenburgerBuchtWismarBucht",][1:ncm1],"MecklenburgerBucht")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)

df_tmp <- cbind(df_clo04[df_clo04$locality10=="NGermanyMecklenburgerBuchtWismarBucht",][1:ncm1],"MecklenburgerBucht")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)


df_tmp <- cbind(df_clo04[df_clo04$locality10=="GermanyBusum",][1:ncm1],"WaddenSeaBussumHaupstr")
colnames(df_tmp) <- colnames(df_clo04) 
df_clo04 <- rbind(df_tmp,df_clo04)

# add a column with location names
h.pca.scores$pop <- poploc
# order the data frame by the column with location names
h.pca.scores <- h.pca.scores[order(h.pca.scores$pop),]
# get unique location names
ppll3 <- unique(h.pca.scores$pop)
df_clo04$locality11 <- df_clo04$locality9
df_clo04$locality11 <- gsub("Samsoe, Ballen" ,"Samsøe, Ballen", df_clo04$locality11)
df_clo04$locality11 <- gsub("Funen, Bogense," ,"Funen, Bogense", df_clo04$locality11)
df_clo04$locality11 <- gsub("Germany, Busum" ,"Germany, Büsum", df_clo04$locality11)
df_cll$coll_locAbb <- df_clo04$locality8[match(df_cll$coll_loc,df_clo04$locality11)]
ppll3 <- gsub("USA:P","USA:W",ppll3)
ppll3 <- gsub("USA:G","USA:W",ppll3)
dd$loc4 <- dd$loc3
dd$loc4 <-  gsub("Samsøe, Ballen","Ballen",dd$loc4)
dd$loc4 <-  gsub("Sealand, Skovshoved","Skovshoved",dd$loc4)
dd$loc4 <-  gsub("Funen, Bogense,","Bogense",dd$loc4)
dd$loc4 <-  gsub("Funen, Kerteminde","Kerteminde",dd$loc4)
dd$loc4 <-  gsub("Jutland, Limfjord","Limfjord",dd$loc4)
dd$loc4 <-  gsub("Jutland, Mariager Fjord","Mariagerfjord",dd$loc4)
dd$loc4 <-  gsub("Mecklenburger","MecklenburgerBuchtWismarBucht",dd$loc4)
dd$loc4 <-  gsub("North Sea, Helgoland Roads","NWGermanyNSeaHelgolandRds",dd$loc4)
dd$loc4 <-  gsub("Germany, Kiel Fjord","GermanyKielFjord",dd$loc4)
dd$loc4 <-  gsub("Germany, Büsum","GermanyBusum",dd$loc4)
dd$AbNm4 <- df_clo04$locality8[match(dd$loc4,df_clo04$locality10)]
# now match between data frames to get colors for locations
cl7 <- df_cll$colfcol_loc[match(ppll3,df_cll$coll_loc)]
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
#p01 <- p01 + scale_color_manual(values = cl7) 
p01 <- p01 + scale_fill_manual(values = cl7) 
# change the heading for the legend, this must be done for all 
#settings for the points
p01 <- p01 + labs(color='Location')
p01 <- p01 + labs(fill='Location')
p01 <- p01 + labs(shape='Location')
#p01

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
#
# set the color of the fill for lines and the fill for inside fill for the points
p02 <- p02 + scale_color_manual(values = colfh_y) 
p02 <- p02 + scale_fill_manual(values = colfh_y) 

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
figname01 <- paste0("Fig06_v01_pca_M_leyidy.png")
figname01 <- paste0("Fig06_v01_pca_M_leyidy.jpg")
figname02 <- paste(wd00_wd05,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}


#_______________________________________________________________________________
#_______________________________________________________________________________
figname01 <- "Fig03_v02_network.jpg"
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
gl_pip <- dartR::gi2gl(gei_pip3)
gei_pip2 <- gei_pip3
pop(gei_pip2) <- gsub("^(.*)_(.*)_(.*)$","\\2",rownames(gei_pip3$tab))
gl_pip <- dartR::gi2gl(gei_pip2)
poploc2 <- gsub("^(.*)_(.*)_(.*)$","\\2",rownames(gei_pip3$tab))
pop(gl_pip) <- poploc2
gl_pip$ind.names <- gsub("^(.*)_(.*)_(.*)$","\\1",rownames(gei_pip2$tab))
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
# make another PCA plot fig06_v02 - start
#_______________________________________________________________________________

# convert genind object to genlight object with the dartR package
gl_pip <- dartR::gi2gl(gei_pip2)
# get sample numbers from haplotype table
liht2Nm <- strsplit(as.character(gl_pip@ind.names), "_")
# get the first element of the split string
liht2Nm.1 <- sapply(liht2Nm, "[[", 1)
liht2Nm.2 <- sapply(liht2Nm, "[[", 2)
liht2Nm.3 <- sapply(liht2Nm, "[[", 3)
smplnt2 <- liht2Nm.2
popyear2 <- liht2Nm.3
# Add predetermined populations. by location
pop(gl_pip) <- poploc2
# make the dapc object
pnw.dapc.ml <- adegenet::dapc(gl_pip, n.pca = 3, n.da = 2) 
# Add predetermined populations. by year
pop(gl_pip) <- popyear2
# make the dapc object
pnw.dapc.my <- adegenet::dapc(gl_pip, n.pca = 3, n.da = 2) # warty comb jelly

# make filename to save the file to
figname01 <- paste0("Fig06_v02_dpca_location_M_leyidy.jpg")
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
# check if DAPC is similar to the PCA we can plot the data in a scatter plot.
scatter(pnw.dapc.ml, col = cl7, cex = 2.4, legend = TRUE, 
        clabel = T, posi.leg = "bottom", scree.pca = TRUE, 
        posi.pca = "bottomleft", posi.da="bottomright", cleg = 1.4)



# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 
#_______________________________________________________________________________
# make another PCA plot fig06_v02 - end
#_______________________________________________________________________________

# make filename to save the file to
figname01 <- paste0("Fig06_v03_dpca_year_M_leyidy.jpg")
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")

# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 5000, height = 4000, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 2), xpd=FALSE)
par(oma=c(1, 1, 0, 0))
#reset this parameter
par(mfrow = c(1, 1)) 
# check if DAPC is similar to the PCA we can plot the data in a scatter plot.
scatter(pnw.dapc.my, col = colfh_y, cex = 2.4, legend = TRUE, 
        clabel = T, posi.leg = "top", scree.pca = TRUE, 
        posi.pca = "topleft", posi.da="topright", cleg = 1.4)



# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1))




#_______________________________________________________________________________
# make tables with PhiST values

#_______________________________________________________________________________
#
df.p.lo <- df_pw_pip_lo
df.p.ye <- df_pw_pip_ye

df.p.lo <- round(df.p.lo,2)
df.p.ye <- round(df.p.ye,2)

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
"GWi",
"JLi",
"JLo",
"GHe",
"GBu",
"SBa",
"GWi",
"SSk")
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
df_clo03$locality8[df_clo03$locality8=="GMe"] <- "GWi"


# paste together location and abbreviation to use for the tables
v_locabb <- paste0(df_clo03$locality4," (",df_clo03$locality8,");")
# exclude by grepping
v_locabb <- v_locabb[!grepl("G:M",v_locabb)]
v_locabb <- v_locabb[!grepl("Panacea",v_locabb)]
v_locabb <- v_locabb[!grepl("GalvestonBay",v_locabb)]
v_locabb <- v_locabb[!grepl("Germany:Helgoland",v_locabb)]
v_locabb <- v_locabb[!grepl("Germany:Maasholm",v_locabb)]
v_locabb <- v_locabb[!grepl("Germany:KielFjord ",v_locabb)]
# use substitute to replace in names
v_locabb <- gsub("BalticSea","Baltic Sea",v_locabb)
v_locabb <- gsub("CaspianSea","Caspian Sea",v_locabb)
v_locabb <- gsub("CentralWAtlantic","Central W Atlantic",v_locabb)
v_locabb <- gsub("NEAtlantic","NE Atlantic",v_locabb)
v_locabb <- gsub("AtlanticOcean:NWAtlantic","NW Atlantic",v_locabb)
v_locabb <- gsub("Ssk","SSk",v_locabb)

#ensure only unique names are included
v_locabb <- unique(v_locabb)

# paste all elements in vector together in one string
locat_abb <- paste(v_locabb,collapse=" ")
# replace the last semi colon with a punctuation mark and a space
locat_abb <- gsub(";$",". ",locat_abb)
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
# substitue in column names and in row names
rownames(mdlo) <- gsub("Ssk","SSk",rownames(mdlo))
colnames(mdlo) <- gsub("Ssk","SSk",colnames(mdlo))
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
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)

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
capt_tbl01 <- paste0("Table 1. Comparison of population genetic fixation index (PhiST) obtained from sequences of nDNA ITS1-2 for various sampling locations. Low PhiST values suggest that individuals sampled are genetically very similar, and high PhiST values indicates genetic variation is high. Abbreviations for locations are: ",
                     locat_abb,
                     "The color gradient reflects the PhiST index with a yellow for low PhiST values, and dark blue for high PhiST values.")

#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.2f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.loc <- mdlo2 
# show the table
t.HTML01 <- mdlo2 %>%
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
d1 <- df.p.ye
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
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = 1)
# reverse the color scale: see: https://stackoverflow.com/questions/45868625/how-to-reverse-the-default-color-palette-for-ggplot2
colscale <- scale_fill_paletteer_c(palette="pals::ocean.haline", direction = -1)
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
capt_tbl02 <-        "Table 2. Comparison of PhiST for Mnemiopsis leidyi obtained from sequences of nDNA ITS1-2 for the sampling years. Sequences of nDNA ITS1-2 obtained from the National Center for Biotechnology Information (NCBI) GenBank were assigned a sampling year as inferred from the year the sequence was deposited on NCBI GenBank. A high PhiST index indicates no variation among the individuals sampled, a high PhiST index indicates there ishigh variation among sequences compared. The color gradient reflects this as a yellow for low PhiST values, and dark blue for high PhiST values."
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr <- function(c) sprintf("%.2f", c)
# apply the function to the matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mdlo2 <-  apply(mdlo, 2, fSpr)
# add back row names
rownames(mdlo2) <- rownames(mdlo)
mdlo2.yea <- mdlo2
# show the table
t.HTML02 <- mdlo2 %>%
  addHtmlTableStyle(align = "r") %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  htmlTable(caption = capt_tbl02)
t.HTML02

#_______________________________________________________________________________
# Get lower triangles from tables
#_______________________________________________________________________________

#http://www.sthda.com/french/wiki/ggplot2-heatmap-d-une-matrice-de-corr-lation-logiciel-r-et-visualisation-de-donn-es
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
upper_tri <- get_upper_tri(mdlo2.loc)


#_______________________________________________________________________________
# Make MDS plot
#_______________________________________________________________________________
#https://www.statmethods.net/advstats/mds.html

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

#count number of columns in data frame
nc09 <-ncol(df_hap_loc09)
#get columns from the second the 7th last column 
hpl10 <- df_hap_loc09[,2:(nc09-6)]
# match to get the other set of abbreviations for sampling locations
df_clo04$locality14 <- NULL
df_clo04$locality14 <- df_clo04$locality11
df_clo04$locality14 <- gsub("Funen, Bogense","Funen, Bogense,",df_clo04$locality14)
df_clo04$locality14 <- gsub("Samsøe, Ballen","Samsoe, Ballen",df_clo04$locality14)
df_clo04$locality14 <- gsub("Germany, Büsum","Germany, Busum",df_clo04$locality14)
#df_hap_loc09$dec_loc3 <- df_clo04$locality8[match(df_hap_loc09$lnNm,df_clo04$locality14)]
# assign row names to the data frame
rownames(hpl10) <- df_hap_loc09$dec_loc3
d <- dist(hpl10) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="n")
text(x, y, labels = row.names(hpl10), cex=.7) 

# Nonmetric MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
library(MASS)
d <- dist(hpl10) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
#fit # view results

# end plot
dev.off()
# define output file name for plot
figname01 <- "Fig07_v01_NMDS_plot.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 3600, height = 3600, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

# try this example code
# From this website: https://stackoverflow.com/questions/12302366/positioning-axes-labels/12302557#12302557
plot(1:100, cumsum(rnorm(100)), type="l", mgp=c(2.4,0.2,.5), las=1)
#and adjust the 'mgp=c(2.4,0.2,.5)'
# To see the effects

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(4, 4, 2, 2), xpd=FALSE)
par(oma=c(0, 0, 0, 0))
#reset this parameter
par(mfrow = c(1, 1)) 
# begin plot, with defined borders
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="NMDS1", ylab="NMDS2",
     #main="Nonmetric MDS",
     #pch=21,bg="red"
     type="n",
     mgp=c(2.4,1.1,0.00001), las=1
     )
text(x, y, labels = row.names(hpl10), cex=1.8) 
#axis(side = 1, col="black")
#axis(side = 2, col="black")
# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 


NMDS1 <- x
NMDS2 <- y
hlp10.lbls <- row.names(hpl10)
df_NMDS_01 <- as.data.frame(cbind(NMDS1,NMDS2,hlp10.lbls))
ggplot(df_NMDS_01) +
  geom_point(aes(x=NMDS1,
             y=NMDS2))

#_______________________________________________________________________________

#_______________________________________________________________________________
#https://www.r-statistics.com/2016/01/multidimensional-scaling-with-r-from-mastering-data-analysis-with-r/

dlo <- df_pw_pip_lo
dly <- df_pw_pip_ye
#read FASTA as dna.bin
pip3 <- ape::read.dna(pth_inpf03, format = "fasta")
#make the DNAbin object a genind object
geni_pip3 <- adegenet::DNAbin2genind(pip3)
#make the genind object a dataframe
df_pip3 <- adegenet::genind2df(geni_pip3)
#get the row names
orig_rwnm <- row.names(df_pip3)
# replace all NAs
df_pip3 <- df_pip3 %>% replace(is.na(.), "-")
# grep for "Bolinopsis" to exclude from the data frame
df_pip3 <- df_pip3[!grepl("Bolinopsis",rownames(df_pip3)),]
#make the date frame a matrix and a DNAbin object again
dnb_pip3 <- as.DNAbin(as.matrix(df_pip3))
# get abbreviated location name
ablo <- gsub("^(.*)_(.*)_(.*)$","\\2",row.names(dnb_pip3))
#rownames(pip) <- smplloca
d <- dlo
#d <- dpip
library(MASS)
#ape::dist.dna(dlo)
d <- dist(dlo) # euclidean distances between the rows

#d <- dist(hpl10) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim

# make the fitted NMDS a data frame
fmds <- as.data.frame(fit)

minstr <- min(fmds$stress)
maxstr <- max(fmds$stress)
strslvl <- ifelse(minstr==maxstr,maxstr, paste0(minstr,"-",maxstr) )
strslvl <- round(strslvl,2)
# make an empty column for overall location 
df_clo03$ov.loc <- NA
df_clo03$ov.loc[!is.na(df_clo03$locality)] <- "NE Atlantic"
df_clo03$ov.loc[grepl("USA",df_clo03$locality4)] <- "NW Atlantic"
df_clo03$ov.loc[grepl("NEAtlantic",df_clo03$locality4)] <- "NE Atlantic"
df_clo03$ov.loc[grepl("Mediterranean",df_clo03$locality4)] <- "Mediterranean"
df_clo03$ov.loc[grepl("Netherlands",df_clo03$locality4)] <- "NE Atlantic"
df_clo03$ov.loc[grepl("Baltic",df_clo03$locality4)] <- "NE Atlantic"
df_clo03$ov.loc[grepl("Caspian",df_clo03$locality4)] <- "Caspian Sea"
df_clo03$ov.loc[grepl("NWAtlantic",df_clo03$locality4)] <- "NW Atlantic"
df_clo03$ov.loc[grepl("Germany",df_clo03$locality4)] <- "NE Atlantic"
df_clo03$ov.loc[grepl("CentralWAtlantic",df_clo03$locality4)] <- "Central W Atlantic"
# replace long row names with abbreviated location names
#rownames(fmds) <- df_clo03$locality8[match(rownames(fmds),df_clo03$locality6)]
# see if there is a match in one column, and use this to add a string to another
# column, that specifies the overall sampling location
df_clo03$ov.loc[df_clo03$locality9=="CW Atlantic"] <- "CW Atlantic" 
df_clo03$ov.loc[df_clo03$locality9=="NE Atlantic"] <- "NE Atlantic" 
df_clo03$ov.loc[df_clo03$locality9=="NW Atlantic"] <- "NW Atlantic" 
# match to get overall location
fmds$ov.loc <- df_clo03$ov.loc[match(rownames(fmds),df_clo03$locality8)]
# define path and output file name
outfNm3 <- paste0(wd00_wd05,"/df_clo03.csv")
# write output filename
write.table(df_clo03,outfNm3,sep=";")

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
nrwsfmds <- length(rownames(fmds))
# plot it with ggplot
p14 <- ggplot(fmds, aes(points.1, -points.2, label = rownames(fmds))) +
  #geom_jitter() +
  geom_point(aes(shape= ov.loc, 
                 fill=ov.loc),size=3.0) +  
  geom_text(check_overlap = TRUE,
            #position = position_dodge(width = 1),
            vjust = -0.75) #+ theme_minimal() + xlab('') + ylab('') +
#scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL)

# add text to the plot to indicate stress level
p14 <- p14 +                               # Add text element to plot
  annotate("text", x = 0, y = -5.5, label = paste0("stress level: ",strslvl))

p14 <- p14 +  scale_shape_manual(values = c(as.integer(pchnmbs))) +
  #scale_color_manual(values = c(rep("#00AFBB", nrwsfmds))) +
  scale_fill_manual(values = c(colnmbs))

p14 <- p14 + theme(panel.background = element_rect(fill = 'white', 
                                                   color = 'white'))#,
#rect = element_line(color = 'black')) #, linetype = 'dotted'))#,
# change label for legend - Notice that you need to change for all 3 variables
# you called 'aes' in 'geom_jitter'
p14 <- p14 + labs(fill='Location')
p14 <- p14 + labs(color='Location')
p14 <- p14 + labs(shape='Location')
p14 <- p14 + xlab("NMDS1") + ylab("NMDS2")

# add border around plot
p14 <- p14 + theme(panel.border = element_rect(color = "black",
                                               fill = NA,
                                               size = 1.0))
# change background of legend
p14 <- p14 + theme(legend.key = element_rect(fill = "white"))

# https://stackoverflow.com/questions/66102618/ggplot-increasing-the-distance-between-axis-labels-and-axis-ticks
p14 <- p14 + theme(axis.text.x = element_text(color="black", 
                                              size=12.4,
                  margin = margin(t = 2, 
                                  unit = "mm")
                                              ),
                   
                   axis.text.y = element_text(color="black", 
                                              size=12.4,
                   margin = margin(l = 0, r=2,
                                   unit = "mm")))
#make filename to save plot to
figname14 <- paste0("Fig07_v02_smpl_location_NMDS_",inpf01,".jpg")

p14

figname02 <- paste(wd00_wd05,"/",figname14,sep="")
if(bSaveFigures==T){
  ggsave(p14,file=figname02,
         #width=210,height=297,
         width=210,height=(297*0.5),
         units="mm",dpi=300)
}
#_______________________________________________________________________________
#_______________________________________________________________________________

#_______________________________________________________________________________


dly <- df.p.ye
# exclude the years that are too variable
yrs_to_excl <- as.character(paste(c("unknown",2007),collapse = "|"))
indx.r.tk <- which(!grepl(yrs_to_excl,row.names(dly)))
indx.c.tk <- which(!grepl(yrs_to_excl,colnames(dly)))
dly <- dly[indx.r.tk,indx.c.tk]
library(MASS)
d <- dist(dly) # euclidean distances between the rows

fit <- isoMDS(d, k=2, p=80, maxit=100, tol=1e-3) # k is the number of dim
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
#colnmbs <- rep(c("white","white","white","black","black"),3)
#colnmbs <- rep(c("white","black"),6)
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
fmds


# plot it with ggplot
p15 <- ggplot(fmds, aes(points.1, -points.2, label = rownames(fmds))) +
  #geom_jitter() +
  geom_point(aes(shape= smplyear, 
                 fill=smplyear),size=3.0) +  
  # geom_text(check_overlap = TRUE,
  #           #position = position_dodge(width = 1),
  #           vjust = -0.75) #+ theme_minimal() + xlab('') + ylab('') +
  
  geom_text_repel()
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
p15 <- p15 + theme(axis.text.x = element_text(color="black", 
                                              size=12.4,
                                              margin = margin(t = 2, 
                                                              unit = "mm")
),

axis.text.y = element_text(color="black", 
                           size=12.4,
                           margin = margin(l = 0, r=2,
                                           unit = "mm")))
#make filename to save plot to
figname15 <- paste0("Fig07_v04_smpl_year_NMDS_",inpf01,".jpg")
# define pathe and file together in a string
figname02 <- paste(wd00_wd05,"/",figname15,sep="")
# make an if test to check whether the plot should be saved
# i.e. set to FALSE if you do not want the plot to be saved
if(bSaveFigures==T){
  ggsave(p15,file=figname02,
         #width=210,height=297,
         width=210,height=(297*0.5),
         units="mm",dpi=300)
}



#_______________________________________________________________________________
# Fig07_v03 - start
#_______________________________________________________________________________
#read FASTA as dna.bin
pip3 <- ape::read.dna(pth_inpf03, format = "fasta")
#make the DNAbin object a genind object
geni_pip3 <- adegenet::DNAbin2genind(pip3)
#make the genind object a dataframe
df_pip3 <- adegenet::genind2df(geni_pip3)
#get the row names
orig_rwnm <- row.names(df_pip3)
# replace all NAs
df_pip3 <- df_pip3 %>% replace(is.na(.), "-")
# grep for "Bolinopsis" to exclude from the data frame
df_pip3 <- df_pip3[!grepl("Bolinopsis",rownames(df_pip3)),]
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
# make the fitted NMDS a data frame

lcNm0 <- labels(dnb_pip3)
lcNm0 <- strsplit(lcNm0,"_")
lcNm1 <- sapply(lcNm0, "[[", 1)
lcNm2 <- sapply(lcNm0, "[[", 2)
lcNm3 <- sapply(lcNm0, "[[", 3)

shNm3 <- gsub("^(.*)_(.*)_(.*)$","\\2",labels(dnb_pip3))
#fmMDS$species <- lcNm4
fmMDS$species <- shNm3
#row.names(fmMDS$points) <- lcNm4
row.names(fmMDS$points) <- shNm3
strslvl2 <- fmMDS$stress
strslvl2 <- round(strslvl2,3)

# define output file name for plot
figname01 <- "Fig07_v03_smpl_location_NMDS_plot.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 3600, height = 3600, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(5, 5, 2, 2), xpd=FALSE)
par(oma=c(0, 0, 0, 0))
#reset this parameter
par(mfrow = c(1, 1)) 
#dev.off()
ordiplot(fmMDS,type="n",cex.axis=1.8,cex.lab = 1.8)
#ordiellipse(fmMDS,groups=lcNm4,draw="polygon",col="grey90",label=F)
text(x = 0.2, y = -0.4, cex=1.8,                # Add text element
     paste0("stress level: ",strslvl2))
#orditorp(fmMDS,display="species",col="red",air=0.01)
orditorp(fmMDS,display="sites",#col=c(rep("green",5),rep("blue",5)),
         air=0.2,cex=1.25)



# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 
#_______________________________________________________________________________
# Fig07_v03 - end
#_______________________________________________________________________________

#_______________________________________________________________________________
# Fig07_v05 - start
#_______________________________________________________________________________
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
# smplyown5 <- sapply(lcNm0[grepl("Mnelei",lcNm0)], "[[", 4)
# smplNown5<- sapply(lcNm0[grepl("Mnelei",lcNm0)], "[[", 1)
# #combine to a data frame
# df_smplown5 <- as.data.frame(cbind(smplNown5, smplyown5))
# #match to get sampling year
# df_smp4$smply[is.na(df_smp4$smply)] <- df_smplown5$smplyown5[match(df_smp4$lcNm1[is.na(df_smp4$smply)],df_smplown5$smplNown5)]
# 
# df_smp5 <- df_smp4
# colnames(df_smp5) <- c("smplyer","longNm")
fmMDS$species <- df_smp4$smply
lcNm5 <- df_smp4$smply
row.names(fmMDS$points) <- lcNm5
strslvl2 <- fmMDS$stress
strslvl2 <- round(strslvl2,3)


# define output file name for plot
figname01 <- "Fig07_v05_smpl_year_NMDS_plot.jpg"
pthfignm01 <- paste(wd00_wd05,"/",figname01,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm01,width = 3600, height = 3600, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#add extra space to the right of the plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
par(mar=c(5, 5, 2, 2), xpd=FALSE)
par(oma=c(0, 0, 0, 0))
#reset this parameter
par(mfrow = c(1, 1)) 
#dev.off()
ordiplot(fmMDS,type="n",cex.axis=1.8,cex.lab = 1.8)
#ordiellipse(fmMDS,groups=lcNm4,draw="polygon",col="grey90",label=F)
text(x = 0.0, y = -0.4, cex=1.8,                # Add text element
     paste0("stress level: ",strslvl2))
#orditorp(fmMDS,display="species",col="red",air=0.01)
orditorp(fmMDS,display="sites",#col=c(rep("green",5),rep("blue",5)),
         air=0.2,cex=1.25)



# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 
#_______________________________________________________________________________
# Fig07_v05 - end
#_______________________________________________________________________________

#_______________________________________________________________________________
# ANOSIM test on locations 01 - start
#_______________________________________________________________________________
#  read in the fasta file which makes it a DNAbin object
pip4 <- read.dna(file=pth_inpf03,format="fasta")
#make the DNAbin object a genind object
geni_pip4 <- adegenet::DNAbin2genind(pip4)
#make the genind object a dataframe
df_pip4 <- adegenet::genind2df(geni_pip4)
#get the row names
orig_rwnm <- row.names(df_pip4)
# replace all NAs
df_pip4 <- df_pip4 %>% replace(is.na(.), "-")
#make the date frame a matrix and a DNAbin object again
dnb_pip4 <- as.DNAbin(as.matrix(df_pip4))
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
pip5.ano <- anosim(m_comp5,grp5)
# turn off previous plots
dev.off()
# plot the test
# define name for outpu plot file
figname08 <- "Fig08_v01_ANOSIMplot_location_all.jpg"
pthfignm08 <- paste(wd00_wd05,"/",figname08,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm08,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#summary(pip5.ano)
par(mar=c(4, 5, 1, 1), xpd=FALSE,
    oma=c(0,0,0,0),
    mfrow = c(1, 1))

plot(pip5.ano,
     xlab="Location",
     ylab="Dissimilarity ranks between \n and within classes")
# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
# ANOSIM test on locations 01 - end
#_______________________________________________________________________________

#_______________________________________________________________________________
# ANOSIM test on locations 02 - start
#_______________________________________________________________________________
#  read in the fasta file which makes it a DNAbin object
pip4 <- read.dna(file=pth_inpf03,format="fasta")
#make the DNAbin object a genind object
geni_pip4 <- adegenet::DNAbin2genind(pip4)
#make the genind object a dataframe
df_pip4 <- adegenet::genind2df(geni_pip4)
#get the row names
orig_rwnm <- row.names(df_pip4)
# split the string with the row name
# and rowbind the nested lists to a dataframe
# to get the different elements
df_pp4Nm <- data.frame(do.call
                       ('rbind', 
                         strsplit(as.character(orig_rwnm),
                                  "_")))
# modify the column names to something more meaningful
colnames(df_pp4Nm) <- c("smplNm",
                        "smplloc",
                        "smplyear")
# identify the unique column names, and use dput to limit the vector
# manually afterwards to only include the NE Atlantic samples
locNm5 <- unique(df_pp4Nm$smplloc)
locNm5 <- locNm5[order(locNm5)]
#dput(locNm5)
# make a vector with sample location names to keep
smplloc.tk <- c("B", "FBo", "FKe", "GBu", "GHe", "GKi", "GWi", 
                "JLi", "JMa", "NEA", "SBa", "SSk")
# use this vector to only keep the rows that has the NE Atlantic samples
df_pip4 <- df_pip4[(df_pp4Nm$smplloc %in% smplloc.tk),]
# replace all NAs
df_pip4 <- df_pip4 %>% replace(is.na(.), "-")
#make the date frame a matrix and a DNAbin object again
dnb_pip4 <- as.DNAbin(as.matrix(df_pip4))
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
pip5.ano <- anosim(m_comp5,grp5)
# turn off previous plots
dev.off()
# plot the test
# define name for outpu plot file
figname08 <- "Fig08_v02_ANOSIMplot_location_NEAtl_smpls.jpg"
pthfignm08 <- paste(wd00_wd05,"/",figname08,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm08,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#summary(pip5.ano)
par(mar=c(4, 5, 1, 1), xpd=FALSE,
    oma=c(0,0,0,0),
    mfrow = c(1, 1))

plot(pip5.ano,
     xlab="Location",
     ylab="Dissimilarity ranks between \n and within classes")
# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
# ANOSIM test on locations 02 - end
#_______________________________________________________________________________


#_______________________________________________________________________________
# ANOSIM test on year for only NEA samples - start
#_______________________________________________________________________________

# make the data frame a  matrix and then a DNAbin object  
dnb_pip4 <- as.DNAbin(as.matrix(df_pip4))
# make a distance matrix 
dst_pip4 <- ape::dist.dna(dnb_pip4, model= "raw")
# make it a matrix
mtx_pip4 <- as.matrix(dst_pip4)
# make it a data frame
df_pip6 <- as.data.frame(mtx_pip4)
# split the row name string by a character 
lpip6 <- strsplit(as.character(row.names(mtx_pip4)), "_")
# get the 3 rd element per vector in the list - this holds the sample year
grp.yea6 <- sapply(lpip6, "[[", 3)
# order years
grp.y6 <- unique(grp.yea6)[order(unique(grp.yea6))]
# make a sequence of numbers to use as index numbers for the years
nof.l6 <- seq(1,length(grp.y6))
# bind the columns in a data frame that has index numbers per year
df_nfl6 <- as.data.frame(cbind(nof.l6,grp.y6))
# attach the group location back to the data frame
df_pip6$grp.yea <- grp.yea6
# match the location name to get an index number, and ensure this number is numeric
df_pip6$grp.yea <- as.numeric(df_nfl6$nof.l6[match(df_pip6$grp.yea,df_nfl6$grp.y6)])
#make community matrix
m_comp6 <- sapply(df_pip6[,1:ncol(df_pip6)],as.numeric)
# ensure the row names in the matrix are as in the data frame
row.names(m_comp6) <- row.names(df_pip6)
# group by site
df_comp6 <- as.data.frame(m_comp6)
# remove all non duplicated rows - the ANOSIM test requires that there are
# multiple representations per group
df_comp6 <- df_comp6[duplicated(df_comp6[c('grp.yea')]), ]


#row.names(df_comp6)
slpip6 <- strsplit(as.character(row.names(df_comp6)), "_")
# get third elements
grp6 <- as.factor(sapply(slpip6, "[[", 3))

# make it a matrix
m_comp6 <- as.matrix(df_comp6)
# make the ANOSIM test
pip6.ano <- anosim(m_comp6,grp6)

dev.off()
# plot the test
# define name for outpu plot file
figname08 <- "Fig08_v03_ANOSIMplot_year_for_only_NEA_samples.jpg"
pthfignm08 <- paste(wd00_wd05,"/",figname08,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm08,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#summary(pip5.ano)
par(mar=c(4, 5, 1, 1), xpd=FALSE,
    oma=c(0,0,0,0),
    mfrow = c(1, 1))

plot(pip6.ano,
     xlab="Sampling year",
     ylab="Dissimilarity ranks between \n and within classes")
# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
# ANOSIM test on year for only NEA samples - end
#_______________________________________________________________________________


#_______________________________________________________________________________
# ANOSIM test on locations 01 - start
#_______________________________________________________________________________
#  read in the fasta file which makes it a DNAbin object
pip4 <- read.dna(file=pth_inpf03,format="fasta")
#make the DNAbin object a genind object
geni_pip4 <- adegenet::DNAbin2genind(pip4)
#make the genind object a dataframe
df_pip4 <- adegenet::genind2df(geni_pip4)
#get the row names
orig_rwnm <- row.names(df_pip4)
# replace all NAs
df_pip4 <- df_pip4 %>% replace(is.na(.), "-")
#make the date frame a matrix and a DNAbin object again
dnb_pip4 <- as.DNAbin(as.matrix(df_pip4))
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
grp6 <- as.factor(sapply(slpip5, "[[", 3))
# make it a matrix
m_comp5 <- as.matrix(df_comp5)
# make the ANOSIM test
pip5.ano <- anosim(m_comp5,grp6)
# turn off previous plots
dev.off()
# plot the test
# define name for output plot file
figname08 <- "Fig08_v04_ANOSIMplot_year_for_all_locations.jpg"
pthfignm08 <- paste(wd00_wd05,"/",figname08,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm08,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#summary(pip5.ano)
par(mar=c(4, 5, 1, 1), xpd=FALSE,
    oma=c(0,0,0,0),
    mfrow = c(1, 1))

plot(pip5.ano,
     xlab="Sampling year",
     ylab="Dissimilarity ranks between \n and within classes")
# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
# ANOSIM test on locations 01 - end
#_______________________________________________________________________________


#_______________________________________________________________________________
# ANOSIM test on year for only NEA samples for only 2017-2019 samples - start
#_______________________________________________________________________________

# make the data frame a  matrix and then a DNAbin object  
dnb_pip4 <- as.DNAbin(as.matrix(df_pip4))
# make a distance matrix 
dst_pip4 <- ape::dist.dna(dnb_pip4, model= "raw")
# make it a matrix
mtx_pip4 <- as.matrix(dst_pip4)
# make it a data frame
df_pip6 <- as.data.frame(mtx_pip4)
# split the row name string by a character 
lpip6 <- strsplit(as.character(row.names(mtx_pip4)), "_")
# get the 3 rd element per vector in the list - this holds the sample year
grp.yea6 <- sapply(lpip6, "[[", 3)
# order years
grp.y6 <- unique(grp.yea6)[order(unique(grp.yea6))]
# make a sequence of numbers to use as index numbers for the years
nof.l6 <- seq(1,length(grp.y6))
# bind the columns in a data frame that has index numbers per year
df_nfl6 <- as.data.frame(cbind(nof.l6,grp.y6))
# attach the group location back to the data frame
df_pip6$grp.yea <- grp.yea6
# match the location name to get an index number, and ensure this number is numeric
df_pip6$grp.yea <- as.numeric(df_nfl6$nof.l6[match(df_pip6$grp.yea,df_nfl6$grp.y6)])
#make community matrix
m_comp6 <- sapply(df_pip6[,1:ncol(df_pip6)],as.numeric)
# ensure the row names in the matrix are as in the data frame
row.names(m_comp6) <- row.names(df_pip6)
# group by site
df_comp6 <- as.data.frame(m_comp6)
# remove all non duplicated rows - the ANOSIM test requires that there are
# multiple representations per group
df_comp6 <- df_comp6[duplicated(df_comp6[c('grp.yea')]), ]
# make a sequence of years, and use this to subset the data frame to only 
# the recent sampling
sq_yrs <- seq(2017,2019,1)
sq_yrs_ts <- paste(sq_yrs,collapse = "|")
df_comp6 <- df_comp6[grepl(sq_yrs_ts,row.names(df_comp6)),]
#row.names(df_comp6)
slpip6 <- strsplit(as.character(row.names(df_comp6)), "_")
# get third elements
grp6 <- as.factor(sapply(slpip6, "[[", 3))
# make it a matrix
m_comp6 <- as.matrix(df_comp6)

# make the ANOSIM test
pip6.ano <- anosim(m_comp6,grp6)

dev.off()
# plot the test
# define name for outpu plot file
figname08 <- "Fig08_v05_ANOSIMplot_year_for_only_2017_2019_NEA_samples.jpg"
pthfignm08 <- paste(wd00_wd05,"/",figname08,sep="")
# set to save plot as pdf file with dimensions 8.26 to 2.9
# 8.26 inches and 2.9 inhes equals 210 mm and 74.25 mm
# and 210 mm and 74.25 mm matches 1/4 of a A4 page
#pdf(pthfignm01,width=(1.6*2.9),height=(0.8*8.26))
jpeg(pthfignm08,width = 3600, height = 2400, res =300)#,width=(1.6*2.9),height=(0.8*8.26))

#summary(pip5.ano)
par(mar=c(4, 5, 1, 1), xpd=FALSE,
    oma=c(0,0,0,0),
    mfrow = c(1, 1))

plot(pip6.ano,
     xlab="Sampling year",
     ylab="Dissimilarity ranks between \n and within classes")
# end plot
dev.off()
#reset this parameter
par(mfrow = c(1, 1)) 

#_______________________________________________________________________________
# ANOSIM test on year for only NEA samples for only 2017-2019 samples - end
#_______________________________________________________________________________



# make the dnabin object a phydata object
phyd_pip4 <- as.phyDat(dnb_pip4)
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
dnb_pip5 <- as.DNAbin(phyd_pip5)
# make a distance matrix 
dst_pip4 <- ape::dist.dna(dnb_pip4, model= "raw")
dst_pip5 <- ape::dist.dna(dnb_pip5, model= "raw")
# make it a matrix
mtx_pip4 <- as.matrix(dst_pip4)
mtx_pip5 <- as.matrix(dst_pip5)



# copy data frame
dnb_pip5.1 <- dnb_pip5
ovlcNm5.2 <- gsub("(.*)_(.*)_(.*)","\\2",labels(dnb_pip5))
dnb_pip5.2 <- updateLabel(dnb_pip5.1, labels(dnb_pip5.1), ovlcNm5.2)
#make the DNAbin object a genind object
geni_pip5.2 <- adegenet::DNAbin2genind(dnb_pip5.2)
#make the genind object a dataframe
df_pip5.2 <- adegenet::genind2df(geni_pip5.2)
# replace all NAs
df_pip5.2 <- df_pip5.2 %>% replace(is.na(.), "-")

hfst_pip_lo5.2 <- hierfstat::genind2hierfstat(geni_pip5.2,pop=ovlcNm5.2)
hfst_pip_lo5.2 <- hfst_pip_lo5.2[order(hfst_pip_lo5.2$pop),]
#replace all NAs with 0
hfst_pip_lo5.2 <- hfst_pip_lo5.2 %>% replace(is.na(.), 0)
# as the 'hierfstat::pairwise.WCfst'
# function cannot handle the unique pop rows only  represented by only a single individual
hfst_pip_lo5.2 <- hfst_pip_lo5.2[hfst_pip_lo5.2$pop %in% hfst_pip_lo5.2$pop[duplicated(hfst_pip_lo5.2$pop)],]
#pairwise fst test Nei - this takes an 'hierfstat' object as input
#pw_Nei_pip_un <- hierfstat::pairwise.neifst(hfst_pip_un, diploid = FALSE)
pw_Nei_pip_lo5.2 <- hierfstat::pairwise.neifst(hfst_pip_lo5.2, diploid = FALSE)
df_pw_pip5.2 <- as.data.frame(pw_Nei_pip_lo5.2)
#replace NAs with zeroes
df_pw_pip5.2[is.na(df_pw_pip5.2)] = 0
library(MASS)
#ape::dist.dna(dlo)
dst5.2 <- dist(df_pw_pip5.2) # euclidean distances between the rows
library("ecodist")
#
pip5.2.nmds <- ecodist::nmds(dst5.2, nits=20, mindim=1, maxdim=4)
# get the mean stress level
m.stress <- mean(pip5.2.nmds$stress)
strslvl2 <- round(m.stress,3)

locNm5.2 <- colnames(df_pw_pip5.2)
# choose the best two-dimensional solution to work with
pip5.2.nmin <- min(pip5.2.nmds, dims=2)
pip5.2.vf <- vf(pip5.2.nmin, dst5.2, nperm=1000)



mxlcNm5.2 <- max(unique(as.numeric(factor(locNm5.2))))
rmxlcNm5.2<- ceiling((mxlcNm5.2)/5)
pchsL5.2 <- rep(c(21,22,23,24,25),rmxlcNm5.2)
colL5.2 <- c(rep("white",5),
             rep("black",5),
             rep("blue",5),
             rep("seagreen",5))
pchsL5.2 <- pchsL5.2[1:mxlcNm5.2]
colL5.2 <- colL5.2[1:mxlcNm5.2]
df_pip5.2cl <- as.data.frame(cbind(unique(as.numeric(factor(locNm5.2))),pchsL5.2,colL5.2,unique(locNm5.2)))
colnames(df_pip5.2cl) <- c("locNo","pchsym","col5.2","locNm")

pip5.2lcNm <- cbind(pip5.2.nmin,locNm5.2)
df_pip5.2cl$colfcol_loc <- df_cll$colfcol_loc[match(df_pip5.2cl$locNm,df_cll$coll_loc)]
pchsL5.2 <- df_pip5.2cl$pchsym[match(locNm5.2,df_pip5.2cl$locNm)]
colfL5.2 <- df_pip5.2cl$colfcol_loc[match(locNm5.2,df_pip5.2cl$locNm)]
# get min and max from NMDS data frame
# to be able to plot the stress level on ggplot plot
mny <- min(pip5.2.nmin$X2)
mnx <- min(pip5.2.nmin$X1)
mxy <- max(pip5.2.nmin$X2)
mxx <- max(pip5.2.nmin$X1)
xdf <- mxx-mnx
ydf <- mxy-mny
spot.x.f <- -0.2*xdf
spot.y.f <- 0.7*ydf

library(ggrepel)
# use ggplot to make a nmds plot with all samples
p.nmds.p5.2 <- ggplot(pip5.2.nmin, aes(x=X1,y=X2,
                                       label=locNm5.2,
                                       shape=locNm5.2,
                                       fill=locNm5.2)) +
  theme_classic() +
  #geom_point() +
  geom_jitter(width = 0.0003, height = 0.0001, size = 3) +
  xlab("NMDS1") + ylab("NMDS2")
p.nmds.p5.2 <- p.nmds.p5.2 + 
  #geom_text() +
  # use : https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
  # also see : https://stackoverflow.com/questions/6996538/dynamic-position-for-ggplot2-objects-especially-geom-text
  ggrepel::geom_text_repel(min.segment.length = 21.1,
                           max.overlaps=2,
                           box.padding = 0.3, na.rm = T) +
  ggplot2::scale_shape_manual(values = as.numeric(pchsL5.2) ) +
  ggplot2::scale_fill_manual(values = colfL5.2 ) +
  #http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
  annotate(geom="text", 
                 x=spot.x.f, y=spot.y.f, size= 4,
                 label=paste0("stress level: ",strslvl2),
                 color="black")
# see the nmds plot
p.nmds.p5.2



# p.nmds.p5.2 <- p.nmds.p5.2 + theme(panel.background = element_rect(fill = 'white', 
#                                                                    color = 'white'))#,
p.nmds.p5.2 <- p.nmds.p5.2 + labs(fill='Location')
p.nmds.p5.2 <- p.nmds.p5.2 + labs(color='Location')
p.nmds.p5.2 <- p.nmds.p5.2 + labs(shape='Location')

# https://stackoverflow.com/questions/66102618/ggplot-increasing-the-distance-between-axis-labels-and-axis-ticks
p.nmds.p5.2 <- p.nmds.p5.2 + theme(axis.text.x = element_text(color="black", 
                                                              size=12.4,
                                                              margin = margin(t = 2, 
                                                                              unit = "mm")
),

axis.text.y = element_text(color="black", 
                           size=12.4,
                           margin = margin(l = 0, r=2,
                                           unit = "mm")))

#make filename to save plot to
figname15 <- paste0("Fig07_v06_smpl_loc_NMDS_",inpf01,".png")
# define pathe and file together in a string
figname02 <- paste(wd00_wd05,"/",figname15,sep="")
# make an if test to check whether the plot should be saved
# i.e. set to FALSE if you do not want the plot to be saved
if(bSaveFigures==T){
  ggsave(p.nmds.p5.2,file=figname02,
         #width=210,height=297,
         width=210,height=(297*0.5),
         units="mm",dpi=300)
}



# make a distance matrix 
dst_pip4 <- ape::dist.dna(dnb_pip4, model= "raw")
dst_pip5 <- ape::dist.dna(dnb_pip5, model= "raw")
# make it a matrix
mtx_pip4 <- as.matrix(dst_pip4)
mtx_pip5 <- as.matrix(dst_pip5)



# copy data frame
dnb_pip4.1 <- dnb_pip4
ovlcNm4.2 <- gsub("(.*)_(.*)_(.*)","\\2",labels(dnb_pip4))
dnb_pip4.2 <- updateLabel(dnb_pip4.1, labels(dnb_pip4.1), ovlcNm4.2)
#make the DNAbin object a genind object
geni_pip4.2 <- adegenet::DNAbin2genind(dnb_pip4.2)
#make the genind object a dataframe
df_pip4.2 <- adegenet::genind2df(geni_pip4.2)
# replace all NAs
df_pip4.2 <- df_pip4.2 %>% replace(is.na(.), "-")

hfst_pip_lo4.2 <- hierfstat::genind2hierfstat(geni_pip4.2,pop=ovlcNm4.2)
hfst_pip_lo4.2 <- hfst_pip_lo4.2[order(hfst_pip_lo4.2$pop),]
#replace all NAs with 0
hfst_pip_lo4.2 <- hfst_pip_lo4.2 %>% replace(is.na(.), 0)
# as the 'hierfstat::pairwise.WCfst'
# function cannot handle the unique pop rows only  represented by only a single individual
hfst_pip_lo4.2 <- hfst_pip_lo4.2[hfst_pip_lo4.2$pop %in% hfst_pip_lo4.2$pop[duplicated(hfst_pip_lo4.2$pop)],]
#pairwise fst test Nei - this takes an 'hierfstat' object as input
#pw_Nei_pip_un <- hierfstat::pairwise.neifst(hfst_pip_un, diploid = FALSE)
pw_Nei_pip_lo4.2 <- hierfstat::pairwise.neifst(hfst_pip_lo4.2, diploid = FALSE)
df_pw_pip4.2 <- as.data.frame(pw_Nei_pip_lo4.2)
#replace NAs with zeroes
df_pw_pip4.2[is.na(df_pw_pip4.2)] = 0
library(MASS)
#ape::dist.dna(dlo)
dst4.2 <- dist(df_pw_pip4.2) # euclidean distances between the rows
#dst4.3 <- ape::dist.dna(dnb_pip4) # euclidean distances between the rows
#dst4.3[is.na(dst4.3)] <- 0
library("ecodist")
#
pip4.2.nmds <- ecodist::nmds(dst4.2, nits=20, mindim=1, maxdim=2)
#pip4.3.nmds <- ecodist::nmds(dst4.3, nits=20, mindim=1, maxdim=2)
locNm4.2 <- colnames(df_pw_pip4.2)
#locNm4.3 <- labels(dnb_pip4)
#row.names(df_comp6)
#locNm4.3 <- strsplit(as.character(locNm4.3), "_")
# get third elements
#locNm4.3 <- as.factor(sapply(locNm4.3, "[[", 2))
# choose the best two-dimensional solution to work with
pip4.2.nmin <- min(pip4.2.nmds, dims=2)
#pip4.3.nmin <- min(pip4.3.nmds, dims=2)
pip4.2.vf <- ecodist::vf(pip4.2.nmin, dst4.2, nperm=1000)
#pip4.3.vf <- ecodist::vf(pip4.3.nmin, dst4.3, nperm=1000)


mxlcNm4.2 <- max(unique(as.numeric(factor(locNm4.2))))
#mxlcNm4.3 <- max(unique(as.numeric(factor(locNm4.3))))
rmxlcNm4.2 <- ceiling((mxlcNm4.2)/5)
#rmxlcNm4.3 <- ceiling((mxlcNm4.3)/5)
pchsL4.2 <- rep(c(21,22,23,24,25),rmxlcNm4.2)
#pchsL4.3 <- rep(c(21,22,23,24,25),rmxlcNm4.3)
colL4.0 <- c(rep("white",5),
             rep("black",5),
             rep("blue",5),
             rep("seagreen",5))
pchsL4.2 <- pchsL4.2[1:mxlcNm4.2]
#pchsL4.3 <- pchsL4.3[1:mxlcNm4.3]
colL4.2 <- colL4.0[1:mxlcNm4.2]
#colL4.3 <- colL4.0[1:mxlcNm4.3]
df_pip4.2cl <- as.data.frame(cbind(unique(as.numeric(factor(locNm4.2))),pchsL4.2,colL4.2,unique(locNm4.2)))

# df_pip4.3cl <- as.data.frame(
#   cbind(
#     unique(factor(locNm4.3)),
#     pchsL4.3,
#     colL4.3,
#     as.character(unique(locNm4.3))
#   ))
colnames(df_pip4.2cl) <- c("locNo","pchsym","col4.2","locNm")
#colnames(df_pip4.3cl) <- c("locNo","pchsym","col4.3","locNm")
pip4.2lcNm <- cbind(pip4.2.nmin,locNm4.2)
#pip4.3lcNm <- cbind(pip4.3.nmin,locNm4.3)
df_pip4.2cl$colfcol_loc <- df_cll$colfcol_loc[match(df_pip4.2cl$locNm,df_cll$coll_loc)]
#df_pip4.3cl$colfcol_loc <- df_cll$colfcol_loc[match(df_pip4.3cl$locNm,df_cll$coll_loc)]
pchsL4.2 <- df_pip4.2cl$pchsym[match(locNm4.2,df_pip4.2cl$locNm)]
#pchsL4.3 <- df_pip4.3cl$pchsym[match(locNm4.3,df_pip4.3cl$locNm)]
colfL4.2 <- df_pip4.2cl$colfcol_loc[match(locNm4.2,df_pip4.2cl$locNm)]
#colfL4.3 <- df_pip4.3cl$colfcol_loc[match(locNm4.3,df_pip4.3cl$locNm)]


library(ggrepel)
# use ggplot to make a nmds plot with all samples
p.nmds.p4.2 <- ggplot(pip4.2.nmin, aes(x=X1,y=X2,
                                       label=locNm4.2,
                                       shape=locNm4.2,
                                       fill=locNm4.2)) +
  theme_classic() +
  #geom_point() +
  geom_jitter(width = 0.0003, height = 0.0001, size = 3) +
  xlab("NMDS1") + ylab("NMDS2")
p.nmds.p4.2 <- p.nmds.p4.2 + 
  #geom_text() +
  # use : https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
  # also see : https://stackoverflow.com/questions/6996538/dynamic-position-for-ggplot2-objects-especially-geom-text
  ggrepel::geom_text_repel(min.segment.length = 21.1,
                           max.overlaps=2,
                           box.padding = 0.3, na.rm = T) +
  ggplot2::scale_shape_manual(values = as.numeric(pchsL4.2) ) +
  ggplot2::scale_fill_manual(values = colfL4.2 )
# see the nmds plot
p.nmds.p4.2



# p.nmds.p4.2 <- p.nmds.p4.2 + theme(panel.background = element_rect(fill = 'white', 
#                                                                    color = 'white'))#,
p.nmds.p4.2 <- p.nmds.p4.2 + labs(fill='Location')
p.nmds.p4.2 <- p.nmds.p4.2 + labs(color='Location')
p.nmds.p4.2 <- p.nmds.p4.2 + labs(shape='Location')


# https://stackoverflow.com/questions/66102618/ggplot-increasing-the-distance-between-axis-labels-and-axis-ticks
p.nmds.p4.2 <- p.nmds.p4.2 + theme(axis.text.x = element_text(color="black", 
                                                              size=12.4,
                                                              margin = margin(t = 2, 
                                                                              unit = "mm")
),

axis.text.y = element_text(color="black", 
                           size=12.4,
                           margin = margin(l = 0, r=2,
                                           unit = "mm")))

#make filename to save plot to
figname15 <- paste0("Fig07_v07_smpl_loc_NMDS_",inpf01,".jpg")
# define pathe and file together in a string
figname02 <- paste(wd00_wd05,"/",figname15,sep="")
# make an if test to check whether the plot should be saved
# i.e. set to FALSE if you do not want the plot to be saved
if(bSaveFigures==T){
  ggsave(p.nmds.p4.2,file=figname02,
         #width=210,height=297,
         width=210,height=(297*0.5),
         units="mm",dpi=300)
}

# df_pip4.3cl <- df_pip4.3cl[order(df_pip4.3cl$locNo),]
# clf4.3 <- df_pip4.3cl$colfcol_loc
# pchsL4.3 <- df_pip4.3cl$pchsym
# p.nmds.p4.3 <- ggplot(pip4.3.nmin, aes(x=X1,y=X2,
#                                        label=locNm4.3,
#                                        shape=locNm4.3,
#                                        fill=locNm4.3)) +
#   theme_classic() +
#   #geom_point() +
#   geom_jitter(width = 0.0003, height = 0.0001, size = 3) +
#   xlab("NMDS1") + ylab("NMDS2")
# p.nmds.p4.3 <- p.nmds.p4.3 + 
#   #geom_text() +
#   # use : https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
#   # also see : https://stackoverflow.com/questions/6996538/dynamic-position-for-ggplot2-objects-especially-geom-text
#   ggrepel::geom_text_repel(min.segment.length = 1.41,
#                            max.overlaps=522,
#                            box.padding = 0.103, na.rm = T) +
#   ggplot2::scale_shape_manual(values = as.numeric(pchsL4.3) ) +
#   ggplot2::scale_fill_manual(values = clf4.3 )
# # see the nmds plot
# p.nmds.p4.3
# 
# 
# 
# dst5.3 <- ape::dist.dna(dnb_pip5) # euclidean distances between the rows
# #as.matrix(dnb_pip5)
# dst5.3[is.na(dst5.3)] <- 0
# library("ecodist")
# #head(dst5.3)
# mtxds5.3 <- as.matrix(dst5.3)!=0
# ldst5.3_0 <- length(as.matrix(dst5.3)[as.matrix(dst5.3)==0])
# ldst5.3_0==dim(mtxds5.3)[1]*dim(mtxds5.3)[2]

#_______________________________________________________________________________
# make nucleotide diversity tables - start
#_______________________________________________________________________________
#https://www.rdocumentation.org/packages/strataG/versions/2.4.905
if(!require(strataG)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  # install from GitHub
  devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
}
# get the library for the strataG package
library(strataG)
# copy the data frame
df_pip5 <- df_pip4
# make the data frame a  matrix and then a DNAbin object  
dnb_pip4 <- ape::as.DNAbin(as.matrix(df_pip5))
# split the row name string by a character 
rNmp6 <- strsplit(as.character(row.names(df_pip5)), "_")
# get the 2nd element per vector in the list - this holds the location
lcNm9 <- sapply(rNmp6, "[[", 2)
# count up frequency of samples
df_popinf01 <- as.data.frame(table(lcNm9))
# make the column a character column
df_popinf01$lcNm9 <- as.character(df_popinf01$lcNm9)
# make a vector with total counts
vTtcnt <- c("Total",sum(df_popinf01$Freq))
#  append a row with total counts
df_popinf01 <- rbind(df_popinf01, vTtcnt)
# change column names
colnames(df_popinf01) <- c("popNm","smplSz")
# get unique sampling locations
popNms <- unique(df_popinf01$popNm)
# make empty lists to collect nucleotide diverstity values
# collect haplotype diversity values
# and collect Tajimas D values and probabilities
lst_nd    <- list()
lst_hd    <- list()
lst_nd.sd <- list()
lst_hd.sd <- list()
lst_tD <- list()
lst_tD.pvln <- list()
lst_tD.pvlb <- list()
lst_FusF <- list()
# iterate over sampling location
for (popl in popNms)
{
   if (popl=="Total") {
    popl2 <- "*"
  } else {   popl2 <- popl }
  # subset the data frame to only comprise the current pop name 
  df_psb09 <- df_pip5[grepl(popl2,rownames(df_pip5)),]
  # make the data frame a matrix and make it a DNAbin object
  dnb_ps09 <- ape::as.DNAbin(as.matrix(df_psb09))  
  # get the nucleotide diversity and the haplotype diversity
  nc.dv      <- pegas::nuc.div(dnb_ps09,variance = F)  
  # setting 'variance = T' returns the standard deviation for the subsamples 
  nc.dv.sd   <- pegas::nuc.div(dnb_ps09,variance = T)[2]  
  hc.dv      <- pegas::hap.div(dnb_ps09,variance = F)
  hc.dv.sd   <- pegas::hap.div(dnb_ps09,variance = T)[2]
  # make the data frame a matrix
  mtx_ps09 <- as.matrix(df_psb09)
  # split matrix using row names
  lstv_ps09 <- split(mtx_ps09, f=row.names(mtx_ps09))
  # check if the rows equals the sequence contents for the first sequence
  #unlist(df_psb09[1,]) == unlist(lstv_ps09[1])
  # determine Fus F
  FusF.val <- strataG::fusFs(lstv_ps09)
  # get tajimas D value
  #https://stackoverflow.com/questions/33301632/tajimas-d-for-sequences-of-different-length
  tD         <- pegas::tajima.test(dnb_ps09)[[1]]
  tD.pvln    <- pegas::tajima.test(dnb_ps09)[[2]]
  tD.pvlb    <- pegas::tajima.test(dnb_ps09)[[3]]
  # add results to the empty lists
  lst_nd[[popl]]         <- nc.dv
  lst_nd.sd[[popl]]      <- nc.dv.sd
  lst_hd[[popl]]         <- hc.dv
  lst_hd.sd[[popl]]      <- hc.dv.sd
  lst_tD[[popl]]         <- tD
  lst_tD.pvln[[popl]]    <- tD.pvln
  lst_tD.pvlb[[popl]]    <- tD.pvlb
  lst_FusF[[popl]]    <- FusF.val$gene1
  # end iteration over subsets
}
# make the lists data frames instead of lists
df_ndv <- as.data.frame(do.call(rbind,lst_nd))
df_hdv <- as.data.frame(do.call(rbind,lst_hd))
df_ndv.sd <- as.data.frame(do.call(rbind,lst_nd.sd))
df_hdv.sd <- as.data.frame(do.call(rbind,lst_hd.sd))
df_tD <- as.data.frame(do.call(rbind,lst_tD))
df_tD.pvln <- as.data.frame(do.call(rbind,lst_tD.pvln))
df_tD.pvlb <- as.data.frame(do.call(rbind,lst_tD.pvlb))
df_FusF <- as.data.frame(do.call(rbind,lst_FusF))
# add row names that holds the population location name
df_ndv$poplcNm <- row.names(df_ndv)
df_hdv$poplcNm <- row.names(df_hdv)
df_ndv.sd$poplcNm <- row.names(df_ndv.sd)
df_hdv.sd$poplcNm <- row.names(df_hdv.sd)
df_tD$poplcNm <- row.names(df_tD)
df_tD.pvln$poplcNm <- row.names(df_tD.pvln)
df_tD.pvlb$poplcNm <- row.names(df_tD.pvlb)
df_FusF$poplcNm <- row.names(df_FusF)
# change column names
colnames(df_ndv)       <- c("ndv","poplcNm")
colnames(df_hdv)       <- c("hdv","poplcNm")
colnames(df_ndv.sd)    <- c("ndv.sd","poplcNm")
colnames(df_hdv.sd)    <- c("hdv.sd","poplcNm")
colnames(df_tD)    <- c("tD","poplcNm")
colnames(df_tD.pvln)    <- c("tD.pvln","poplcNm")
colnames(df_tD.pvlb)    <- c("tD.pvlb","poplcNm")
colnames(df_FusF)    <- c("FusF","poplcNm")
#combine data frames in a list
lst_df.nhdv <- list(df_ndv,df_hdv,df_ndv.sd,df_hdv.sd,df_tD,df_tD.pvln,df_tD.pvlb,df_FusF)
#merge all data frames in list
df_nhdv <- lst_df.nhdv %>% purrr::reduce(full_join, by='poplcNm')
# use the DNAbin object to get the frequency of haplotypes per population
mtx_HptFrq <- pegas::haploFreq(dnb_pip4,split="_", what=2)
# set all hapltype counts to 1 , to be able to sum up unique haplotypes
mtx_HptFrq[mtx_HptFrq>=1] <- 1
# get the haplotype frequency - i.e the number of haplotypes per location
hptF.p.loc <- colSums(mtx_HptFrq)
# get the haplotype frequency - i.e the number of haplotypes per haplotype
hptF.p.hpt <- rowSums(mtx_HptFrq)
# only count occurences of haplotypes
hptF.p.hpt[hptF.p.hpt>=1] <- 1
# count the total number of haplotypes
tot.hpt.cnt <- sum(hptF.p.hpt)
# make it a data frame
df_hptF.p.loc <- as.data.frame(hptF.p.loc)
# get the rownames to make a column
df_hptF.p.loc$popNm  <- row.names(df_hptF.p.loc)
# also add the total count of haplotypes
df_hptF.p.loc <- rbind(c(tot.hpt.cnt,"Total"),df_hptF.p.loc)
# match back to get the number of Haplotypes per population
df_popinf01$Nh <- df_hptF.p.loc$hptF.p.loc[match(df_popinf01$popNm,df_hptF.p.loc$popNm)]
# match to get haplotype diverstiy and nucleotide diversity and standard deviations
df_popinf01$hdv       <- df_nhdv$hdv[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf01$hdv.sd    <- df_nhdv$hdv.sd[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf01$ndv       <- df_nhdv$ndv[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf01$ndv.sd    <- df_nhdv$ndv.sd[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf01$tD    <- df_nhdv$tD[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf01$tD.pvln    <- df_nhdv$tD.pvln[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf01$tD.pvlb    <- df_nhdv$tD.pvlb[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf01$tD.pvlb    <- df_nhdv$tD.pvlb[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf01$FusF    <- df_nhdv$FusF[match(df_popinf01$popNm,df_nhdv$poplcNm)]
df_popinf02 <- df_popinf01
#exclude the beta probability
df_popinf02$tD.pvlb <- NULL
# store the populations names in a vector
poplocNms <- df_popinf01$popNm
# add back to row names
rownames(df_popinf02) <- poplocNms
df_popinf02$popNm <- NULL
#df_popinf02 <- as.data.frame(lapply(df_popinf02, function(x) as.numeric(as.character(x))))
rownames(df_popinf02) <- poplocNms
df_popinf02 <- as.data.frame(lapply(df_popinf02, function(x) as.numeric(as.character(x))))
df_popinf02[is.na(df_popinf02)] <- 0
# substitute in the list of abbreviated locations
ablo9 <- gsub("Netherlands \\(Nt\\); ","",ablo8)
ablo9 <- gsub("; USA, WoodsHole \\(USA:W\\)","",ablo9)
# make a table caption
capt_tbl02 <-        paste0("Table 6. Estimates for diversity and neutrality for the populations in the analysed regions of samples for Mnemiopsis leyidi. Abbreviations above columns are: Number of haplotypes per population (Nh), haplotypic diversity (h), and standard deviation for haplotypic diversity (h.sd), nucleotide diversity (nd), and standard deviation for nucleotide diversity (nd.sd), probability of Tajimas D (p.TD). Sampled locations are abbreviated: ",
                      ablo9,". Additional sequences were obtained from M. leidyi from NBCI GenBank as indicated by accession numbers.")
#make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
fSpr2 <- function(c) sprintf("%.2f", c)
fSpr0 <- function(c) sprintf("%.0f", c)
# apply the function to the data frame -  this makes it a matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
mtx_popinf02 <-  apply(df_popinf02, 2,fSpr2)
mtx_popinf02[,1:2] <- as.character(round(as.numeric(as.character(mtx_popinf02[,1:2])),digits=2))

# match to get coordinates for sampling positions
smplcoord <- df_clo04$dec_latlon2[match(poplocNms,df_clo04$locality8)]
smplcoord[is.na(smplcoord)] <- ""
# order data frame by location name and then by individual sample name
df_pip03.2 <- df_pip03[order(df_pip03$locality,df_pip03$smplNm),]
# get lettercode in sample and get number code in samples
df_pip03.2$ltc.smpl <- gsub("^([A-Za-z]{+})([0-9]{+})$","\\1",df_pip03.2$smplNm)
df_pip03.2$noc.smpl <- gsub("^([A-Za-z]{+})([0-9]{+})$","\\2",df_pip03.2$smplNm)
# make the number numeric, to be albe to find the maximum and the minimum
df_pip03.2$noc.smpl <- as.numeric(df_pip03.2$noc.smpl)
# use dplyr to summarise upper and lower value of sample per group
# to use later on in a table that holds all samples
df_pip03.3 <- df_pip03.2 %>%
  # determine the columns to work on with dplyr
  dplyr::select(locality,ltc.smpl,noc.smpl) %>%
  # specify the columns to group by
  dplyr::group_by(locality,ltc.smpl) %>%
  # get the max value of the noc.smpl column per groups defined above, and 
  # put in a new column called 'mxnb', and also do this for the minimum value
  dplyr::summarise(mxnb = max(noc.smpl), minb = min(noc.smpl))
#pad with zeros to three characters for own Mnelei samples
#see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
df_pip03.3$mxnb <- ifelse(nchar(df_pip03.3$mxnb)<3,stringr::str_pad(df_pip03.3$mxnb, 3, pad = "0"),df_pip03.3$mxnb)
df_pip03.3$minb <- ifelse(nchar(df_pip03.3$minb)<3,stringr::str_pad(df_pip03.3$minb, 3, pad = "0"),df_pip03.3$minb)
# paste together to get range of sample numbers per sample location
df_pip03.3$smplrng <- paste0(df_pip03.3$ltc.smpl,df_pip03.3$minb,"-",df_pip03.3$ltc.smpl,df_pip03.3$mxnb)
# concatenate strings per group - see this example: https://datacornering.com/how-to-concatenate-text-by-group-in-r/
# notice I added'dplyr::' to the dplyr functions, as I also have 'plyr' up and running 
require(dplyr)
df_pip03.3 <- df_pip03.2 %>%
  dplyr::select(smplNm, locality) %>% 
  dplyr::group_by(locality) %>%
  dplyr::mutate(all_spmls = paste(smplNm, collapse = ", "))  %>%
  dplyr::distinct(locality,all_spmls)
# copy the column to a vector
all_spmlsNos <- c(df_pip03.3$all_spmls,"")
# get sample month for samples
ncbismplmnth <- df_lN02$smplyear[match(df_pip03$smplNm[grepl("[A-Z]{2}[0-9]{5}",df_pip03$smplNm)],df_lN02$accession_nmb)]
smplmnth.ncbi <- cbind(df_pip03$smplNm[grepl("[A-Z]{2}[0-9]{5}",df_pip03$smplNm)],ncbismplmnth)
df_smplmnth.ncbi <- as.data.frame(smplmnth.ncbi)
colnames(df_smplmnth.ncbi) <- c("smplNo","mnth")
ownsmplmnth <- df_ownsmpldt$smplmnt4[match(df_pip03$smp[grepl("Mnelei",df_pip03$smplNm)],df_ownsmpldt$MnsmNo)]
smplmnth.Mnelei <- cbind(df_pip03$smp[grepl("Mnelei",df_pip03$smplNm)],ownsmplmnth)
df_smplmnth.Mnelei <- as.data.frame(smplmnth.Mnelei)
colnames(df_smplmnth.Mnelei) <- c("smplNo","mnth")
#put all data frames into list
df_list <- list(df_smplmnth.ncbi, df_smplmnth.Mnelei)      
#merge all data frames together
df_smpl.mnth <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)  
# match back to pip dataframe
df_pip03$month <- df_smpl.mnth$mnth[match(df_pip03$smplNm,df_smpl.mnth$smplNo)]
# make an empty column to fill
df_pip03$month2 <- NA
# evaluate if year and month column are the same, and then use the month or nothing
df_pip03$month2[!df_pip03$month==df_pip03$year] <- df_pip03$month[!df_pip03$month==df_pip03$year]
df_pip03$month2[df_pip03$month==df_pip03$year] <- ""
# assign to vectors
syear <- df_pip03$year[match(poplocNms,df_pip03$locality)]
syear[is.na(syear)] <- ""
smnth <- df_pip03$month2[match(poplocNms,df_pip03$locality)]
smnth[is.na(smnth)] <- ""

# bind the population location names back on to the matrix array
mtx_popinf02 <- cbind(poplocNms,smplcoord,mtx_popinf02,syear,smnth,all_spmlsNos)
# add back row names
rownames(mtx_popinf02) <- poplocNms
# change the column names
colnames(mtx_popinf02) <- c("Population","Coordinates","Sample size","Nh","h","h.sd","nd","nd.sd", 
                            "Tajima's D", "p.TD","Fu's F","year","month","sample No")
# count the number of columns in the matrix
nclmtx <- ncol(mtx_popinf02)
#make a vector for text adjustment
vftxtadj <- c("l",rep("r",nclmtx-1))

# show the table
t.HTML05 <- mtx_popinf02 %>%
  addHtmlTableStyle(align = vftxtadj) %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  #addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02, rnames = FALSE)
t.HTML05


#_______________________________________________________________________________
# make nucleotide diversity tables - end
#_______________________________________________________________________________


# use the DNAbin object to get the frequency of haplotypes per population
mtx_HptFrq <- pegas::haploFreq(dnb_pip4,split="_", what=2)

hgstHptNmb <- max(as.numeric(as.roman(as.character(row.names(mtx_HptFrq)))))
row.names(mtx_HptFrq) <- as.numeric(as.roman(as.character(row.names(mtx_HptFrq))))

df_HptFrq2 <- as.data.frame(mtx_HptFrq)
totcnt.Hpt <- colSums(df_HptFrq2)
df_totcnt.Hpt <- as.data.frame(t(totcnt.Hpt))
# also add the total count of haplotypes
df_HptFrq3 <- rbind(df_HptFrq2,df_totcnt.Hpt)

df_HptFrq3$HptNmb <- rownames(df_HptFrq3)

df_HptFrq3$HptNmb[as.numeric(df_HptFrq3$HptNmb)>hgstHptNmb] <- "Total"
lastcol <- ncol(df_HptFrq3)
df_HptFrq3 <- df_HptFrq3[,c(lastcol,seq(1,lastcol-1))]
colnames(df_HptFrq3)[1] <- c("Haplotype") 
df_HptFrq3[df_HptFrq3==0] <- ""

# make a table caption
capt_tbl02 <-        "Table 7. Haplotype frequency per population. Abbreviations for sampling locations are explained in table 1."

# show the table
t.HTML06 <- df_HptFrq3 %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  #addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02, rnames = FALSE)
t.HTML06

#and to export in a file
#tableHTML::write_tableHTML(t.HTML06, file = paste(wd00_wd05,"/Table07_v01_haplotypefrq_",inf01,".html", sep=""))


#_______________________________________________________________________________
#install.packages("pegas")
#packageVersion("pegas")
#R.Version()

# Or the development version from GitHub:
# install.packages("devtools")
#devtools::install_github("r-lib/devtools")

if(!require(pegas)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  # install from GitHub
  devtools::install_github("emmanuelparadis/pegas/pegas")
}
#
packageVersion("pegas")

