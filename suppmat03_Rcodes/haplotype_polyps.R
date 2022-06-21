#This following code is used to plot the haplotype and add the haplotype frequencies for the scyphozoan polyp COI data.
#As imput files are needed:
#--an alignment of sequences in fasta format
#--a .csv file with, for each sequence, the different grouping labels that you want to be applied to the sequence. The sequences should be in exactly the same order as the sequences in the .fas file (if you want to add pie charts of haplotype frequencies)

#the following steps are taken:
#1 read in the data
#2 extract all unique haplotypes and attach custom labels
#3 generate a plain network without fhaplotype frequencies.
#4 optionally the haplotype frequencies can be added as coloured pie charts:
#5 read in label names (location data)
#add names in a single string, separated by an underscore (_)
#6 plot the haplotypes with frequencies
#7 add legend 
#8 resize plot window to the desired size
#9 save plot as pdf
#10 edit network (overlapping nodes, labels etc) in a vector editor such as inkscape or adobe illustrator


#install.packages("pegas")
User <- "Lodewijk"

#set working directory
if (User == "Lodewijk") { setwd ("C:/Users/Desktop_LvW/OneDrive/PhD/publications/polyps")}
if (User == "Henk") { setwd ("H:/Calc")}

wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
wd01 <- "suppmat01_inp_files"
wd05 <- "suppmat05_out_files"
setwd(wd00)
#define paths for in and output directories
wd00_wd01 <- paste0(wd00,"/",wd01)
wd00_wd05 <- paste0(wd00,"/",wd05)
#load packages
library(pegas) #for analyses and haplotype network
library(magrittr)
library(dplyr)

#Important: all these functions rely on the order of the sequences 
#in the fasta file being exactly the same as the order of the labels!!! 
#No error message is given if this is not the case.


#read in the data
#x <- read.dna(file="aligment_alle_Aurelias_definitief.fas",format="fasta")
inf01 <- "algn_Mnelei_18s_10.aligned.fasta.fas"
wd00_wd01_inf01 <- paste0(wd00_wd01,"/",inf01)
x <- read.dna(file=wd00_wd01_inf01,format="fasta")
#extract all unique haplotypes and attach custom labels
h <- haplotype(x)
hlab <- length(labels(h))
h <- haplotype(x,labels=c(1:hlab))#default labels are roman numerals. Vector of names has to be of equal length as the number of haplotypes.
#h_selection<-haplotype(x_selection)
#which polyps belong to which haplotype?
h_indices<-attr(h, "index")
#assign labels
names(h_indices)<-c(1:hlab)
str(h_indices)
#generate plain network 
net <- haploNet(h)
#haploNet
#net <- haploNet(h_selection)
#plot
plot(net,
     size=attr(net, "freq"),
     scale.ratio = 3,
     cex = 0.8,
     show.mutation=1,
     threshold=)
#extract all names
names<-data.frame(labels(x))


#
str(names)

#write.table(names,file="names_def.csv",sep=",")
haploFreq(x,haplo=h)
#get list of all sequence labels with haplotype numbers
haplotypes<-stack(setNames(attr(h, "index"), rownames(h)))

#Grep among the list of files, to get the file that holds the different names
fnm_clo3 <- list.files(wd00_wd05)[grepl("clo03",list.files(wd00_wd05))]
fnm_df_lN02 <- list.files(wd00_wd05)[grepl("df_lN02",list.files(wd00_wd05))]
#read locations file
df_clo3<-read.csv(file=paste0(wd00_wd05,"/",fnm_clo3),header=TRUE,sep=";")
df_lN02<-read.csv(file=paste0(wd00_wd05,"/",fnm_df_lN02),header=TRUE,sep=",")

df_clo3$locality8[match(df_lN02$location,df_clo3$locality4)]
labels(x)
# split string and get lists nested in a list
lbls01 <- strsplit(as.character(labels(x)), "_")
# get only NCBI labels
lbls02 <- sapply(lbls01, "[[", 3)
# grep only elements that includes numbers
ncbilbls <- lbls02[grepl("[0-9]",lbls02)]
# get location names for NCBI accesion numbers
lnbls<- df_lN02$location[match(ncbilbls,df_lN02$accession_nmb)]
# replace the labels that has numbers in them with the location name
lbls02[grepl("[0-9]",lbls02)] <- lnbls

lbls03 <- lbls02
df_lblnms <- as.data.frame(cbind(lbls02,lbls03)) 
colnames(df_lblnms) <- c("longNm","AbrNm")
df_lblnms$AbrNm <- df_clo3$locality8[match(df_lblnms$longNm,df_clo3$locality5)]
llnNm1 <- df_lblnms$longNm[is.na(df_lblnms$AbrNm)]
# make vector with long names appearing in the haplotype labels
Nm_inLonNm <- c("USA:WoodsHole","NEAtlantic","Mediterranean","CaspianSea","CentralWAtlantic","USA:GalvestonBay","USA:Panacea","Germany:KielFjord","Germany:Maasholm","Germany:Helgoland","BalticSea","AtlanticOcean:NWAtlantic","Netherlands","Loegstoer","Mariagerfjord","Kerteminde","NSeaHelgolandRds","WaddenSeaBussumHaupstr","Skovshoved","KielFjord","Bogense","MecklenburgerBuchtWismarBucht","MecklenburgerBucht","Ballen")
# make vector with long names appearing in the df_clo3 df
Nm_inclo03 <- c("USA:WoodsHole","NEAtlantic","Mediterranean","CaspianSea","CentralWAtlantic","USA:GalvestonBay","USA:Panacea","NGermanyKielFjord","NGermanyKielFjord","NWGermanyNSeaHelgolandRds","BalticSea","AtlanticOcean:NWAtlantic","Netherlands","NJyllandLimfjord","JyllandMariagerfjord","FynKerteminde","NWGermanyNSeaHelgolandRds","GermanyBusum","SjaellandSkovshoved","NGermanyKielFjord","FynBogense","NGermanyMecklenburgerBuchtWismarBucht","NGermanyMecklenburgerBuchtWismarBucht","SamsoeBallen")
# combine to a df
df_rplNm <- as.data.frame(cbind(Nm_inLonNm,Nm_inclo03))
# match first the longNm from the Haplotype labels to the df with both the haplotype labels and the df_clo03 labels
df_lblnms$longNm <- df_rplNm$Nm_inclo03[match(df_lblnms$longNm,df_rplNm$Nm_inLonNm)]
# then match the new long name with th df_clo03 table to get the abbreviation
df_lblnms$AbrNm <- df_clo3$locality8[match(df_lblnms$longNm,df_clo3$locality6)]



#add rownames for haploFreq function
rownames(x) <- df_lblnms$AbrNm
#The following code is used to plot the haplotype and add the haplotype frequencies. 

#create haplotypes from dna.bin
#x
px <- pegas::haplotype(x)
#view haplotype 
#h
#as.matrix(pipHaps)
#prepare hpt table
ih<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(x)[values]))
#make it a dataframe
ih01 <- as.data.frame(ih)

#limit to include only 'Freq' that equals 1
ih02 <- ih01[ih01$Freq >= 1,]	

hpfrq <- haploFreq(x,rownames(x),haplo=h)
#hpfrq <- as.data.frame(hpfrq)
#make it a table
hsl <- table(ih02$hap, ih02$pop)
#per area THIS IS THE ONE USED IN THE PAPER
#freq_area<-haploFreq(x,split="_",what=2,haplo=h)
#plot the new haplotype
plot(net,
     size=attr(net, "freq"),
     scale.ratio = 2.7,
     cex = 0.8,
     pie=hsl,
     show.mutation=2,# mutations are shown with small dots on the links
     threshold=0#turn off display of alternative links
     )

legend(30,20, colnames(hsl), col=rainbow(ncol(hsl)), pch=20,bty="n")

#https://cran.r-project.org/web/packages/pegas/vignettes/PlotHaploNet.pdf
#_______________________________________________________________________________
# Advice from Cody Aylward:  - https://doi.org/10.1007/s10592-018-1130-3

# Hello Steen,
# 
# Thanks for your comments on my paper. I looked back at my old hard drive and unfortunately cannot find the exact R script used for that figure. However, I have a few suggestions. I read through this vignette to refresh my memory on constructing haplotype networks in Pegas (https://cran.r-project.org/web/packages/pegas/vignettes/PlotHaploNet.pdf), I found it quite helpful.
# 
# My first suggestion is to re-size of your circles. 
#It looks like you have one haplotype that is substantially more common than the rest, 
#such that it is covering other circles completely. 
#I am assuming you made this figure using the default setting, 
#where the diameter of each circle is equal to the frequency of each haplotype. 
#An alternative is to set the diameter of each circle equal to the square root
#of the frequency of each haplotype divided by pi 
#(such that the area of each circle is equal to the frequency of 
#each haplotype). This is explained on pages 11-12 of the vignette.
#Another alternative is to manually size the circles. 
#This was easier for my paper since I only had 12 haplotypes, 
#and might be a bit more of a task for your larger number of haplotypes. 
#However, if the square root approach doesn't shrink your largest circle 
#sufficiently, you might consider manually creating a vector 
#for your circle sizes. If you use this approach, 
#I would group haplotype frequencies into different size categories. 
#For example, haplotype frequencies of 1-5 would be size 1, haplotype 
#frequencies 6-20 are size 2, etc. 
#That may be a bit more work, but it gives you more direct manipulation 
#of the network. Page 6 of the vignette shows how to order the size
#vector to match the order of your haplotypes in the network.
# 
# I think once the circles are resized your network will be much cleaner. 
#If that isn't sufficient, there may be a way to manually adjust 
#the coordinates of each circle. I did a very quick search and didn't 
#find an answer right away. If you still need to modify the network after 
#re-sizing the circles I would be happy to take a longer look at 
#manually adjusting coordinates.
# 
# Best of luck with your study,
# 
# Cody Aylward, M.S.
#______________________________________________________________________________

# vignette try out

library(pegas) # loads also ape
data(woodmouse)