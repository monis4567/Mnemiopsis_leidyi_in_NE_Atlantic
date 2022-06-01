
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
#as.matrix(pipHaps)
#prepare hpt table
ind.hap<-with(
  stack(setNames(attr(pipHaps, "index"), rownames(pipHaps))),
  table(hap=ind, pop=rownames(pip)[values]))
#make it a dataframe
df_ihpt01 <- as.data.frame(ind.hap)
#limit to include only 'Freq' that equals 1 or more
df_ihpt02 <- df_ihpt01[df_ihpt01$Freq >= 1,]	
# split string and get lists nested in a list
locations <- strsplit(as.character(df_ihpt02$pop), "_")
smplloca <- sapply(locations, "[[", 2)
#Modify location names
smplloca <- gsub("FynBogense","Funen, Bogense", smplloca )
smplloca <- gsub("FynKerteminde","Funen, Kerteminde", smplloca )
smplloca <- gsub("JyllandMariagerfjord","Jutland, Mariager Fjord", smplloca )
smplloca <- gsub("NGermanyKielFjord","Germany, Kiel Fjord", smplloca )
smplloca <- gsub("NGermanyMecklenburgerBuchtWismarBucht","Germany, Wismar Bight", smplloca )
smplloca <- gsub("NJyllandLimfjord","Jutland, Limfjord", smplloca )
smplloca <- gsub("NWGermanyNSeaHelgolandRds", "North Sea, Helgoland Roads", smplloca )
smplloca <- gsub("NWGermanyWaddenSeaBussumHaupstr","Germany, Büsum", smplloca )
smplloca <- gsub("SamsoeBallen","Samsøe, Ballen", smplloca )
smplloca <- gsub("SjaellandSkovshoved","Sealand, Skovshoved", smplloca )
smplloca <- gsub("Germany:Maasholm","Germany, Kiel Fjord" ,smplloca)
smplloca <- gsub("AtlanticOcean:NWAtlantic","NW Atlantic" ,smplloca)
smplloca <- gsub("BalticSea","Baltic Sea" ,smplloca)
smplloca <- gsub("NEAtlantic" ,"NE Atlantic",smplloca)
smplloca <- gsub("CentralWAtlantic" ,"CW Atlantic",smplloca)
smplloca <- gsub("CaspianSea" ,"Caspian Sea",smplloca)

#make it a table
new.hap.smplloc <- table(df_ihpt02$hap, smplloca)

#make a haplotype network
pipNet <- pegas::haploNet(pipHaps)
#plot the network
plot(pipNet, size = attr(pipNet,"freq"), 
     scale.ratio = 1.2, cex = 0.9, pie = new.hap.smplloc, 
     show.mutation = 0, threshold = 0, labels(FALSE))
#add a legend to the plot
# legend("topright",colnames(new.hap.smplloc), 
#        col=rainbow(ncol(new.hap.smplloc)), 
#        pch=19, ncol=1)