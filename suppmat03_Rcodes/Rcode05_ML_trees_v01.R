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
# get the ggtree package
if(!require(ggtree)){
  BiocManager::install("ggtree")
}
if(!require(Biostrings)){
  BiocManager::install("Biostrings")
}
library(Biostrings)
# read in package libraries
library(dplyr)
library(ggtree)
library(pegas)
library(ape)
library(phangorn)
library(Biostrings)
#define overall working directory
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
#define input working directories
wd01 <- "suppmat01_inp_files"
wd05 <- "suppmat05_out_files"

wd_NCBI <- "NCBI_seq_submission_for_Mnemiopsis"
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

# READ in accession numbers

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


# then read in the fasta file which makes it a DNAbinxn object
dnb_pip4 <- read.dna(file=wd00_wd05_inf01,format="fasta")

#make the DNAbin object a genind object
geni_pip4 <- adegenet::DNAbin2genind(dnb_pip4)
#make the genind object a dataframe
df_pip4 <- adegenet::genind2df(geni_pip4)

# To be able to replace the labels in the DNAbin object
# the replacement of labels names has to be done on the data frame version
# of the DNAbin object, afterwards the  data frame version can then be 
# turned back in to a DNAbin object
lbN <- row.names(df_pip4)
# then grep for any names that has 'Mnelei' in the name
lb.Ml <- lbN[grepl("Mnelei",lbN)]
# prepare the name as a data frame that can be accessed
dlb.Ml <- data.frame(do.call
                     ('rbind', 
                       strsplit(as.character(lb.Ml),
                                "_")))
# change the column names of this data frame to something meaningful
colnames(dlb.Ml) <- c("MneleiAbbrNo","locAbbr","smplYr")
# use match to get the corresponding NCBI accession number
dlb.Ml$AccNo <- df_accNos$AccNo[match(dlb.Ml$MneleiAbbrNo,
                                      df_accNos$MneleiAbbr)]
# use the columns in the data frame to paste together a new name
lb.Nw <- paste(dlb.Ml$AccNo,
               dlb.Ml$locAbbr,
               dlb.Ml$smplYr,sep="_")
# then get the names that match 'Mnelei, and replace these names
lbN[grepl("Mnelei",lbN)] <- lb.Nw
# then overwrite the old names with this new vector of new names
row.names(df_pip4) <- lbN

# ::::::::::: IMPORTANT !! START ::::::::::: 
# In order to get the 'pml' function from phangorn  working  
# all cases of NA in the sequences must be replaced, with '-'
# otherwise R will crash
df_pip4[is.na(df_pip4)] <- "-"
# 
# ::::::::::: IMPORTANT !! END ::::::::::: 
#View(df_pip4)
#make the date frame a matrix and a DNAbin object again
dnb_pip4 <- as.DNAbin(as.matrix(df_pip4))
#labels(dnb_pip4)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#https://www.molecularecologist.com/2016/02/26/quick-and-dirty-tree-building-in-r/

if(!require(ape)){
  install.packages("ape")
}
if(!require(phangorn)){
  install.packages("phangorn")
}
if(!require(seqinr)){
  install.packages("seqinr")
}
library(ape)
library(phangorn)
library(seqinr)
library(treeio)
library(ggtree)
library(ggplot2)
# https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html
library(phangorn)
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

# replace the Mnelei abbr numbers with NCBI accession numbers obtained
# after depositing sequences on NCBI
# then get the names
lbN <- df_pp4Nm$smplNm
# then grep for any names that has 'Mnelei' in the name
lb.Ml <- lbN[grepl("Mnelei",lbN)]
# use match to get the corresponding NCBI accession number
lb.AccNo <- df_accNos$AccNo[match(lb.Ml,
                                      df_accNos$MneleiAbbr)]

lbN[grepl("Mnelei",lbN)] <- lb.AccNo
# then overwrite the old names with this new vector of new names
df_pp4Nm$smplNm <- lbN



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
# make the alignment a data matrix
df_phdp4 <- as.data.frame(as.matrix(phyd_pip4))
# count rows and columns
nrd4 <- nrow(df_phdp4)
ncd4 <- ncol(df_phdp4)

library(phangorn)

class(phyd_pip4)
# calculate distances

dm4 <- phangorn::dist.p(phyd_pip4)
dm5 <- phangorn::dist.p(phyd_pip5)

dm4 <- phangorn::dist.ml(phyd_pip4, "F81")
dm5 <- phangorn::dist.ml(phyd_pip5, "F81")
tree_NJ4 <- NJ(dm4)
tree_NJ5 <- NJ(dm5)

# as alternative for a starting tree:
# tree4 <- pratchet(phyd_pip4)          # parsimony tree
# tree5 <- pratchet(phyd_pip5)          # parsimony tree
# tree4 <- nnls.phylo(tree4, dm4)   # need edge weights
# tree5 <- nnls.phylo(tree5, dm5)   # need edge weights

# 1. alternative: quick and dirty: GTR + G
#(fit <- pml_bb(phyd_pip4, model="GTR+G"))
#(fit4 <- pml_bb(phyd_pip4, model="GTR+G"))

# ___________________________________________________
# (fit5 <- pml_bb(phyd_pip4, model="GTR+G"))
# # 2. alternative: choose with modelTest
# # mt4 <- modelTest(phyd_pip4, multicore=TRUE)
# # optmdl4 <- mt4[order(mt4$BIC),][1,1]
# # chooses best model from the table according to BIC (default)
# #fit <- pml_bb(mt)
# #fit_pip4 <- pml(tree_NJ4, phyd_pip4)
# fit_pip5 <- pml(tree_NJ5, phyd_pip5)
# # perform bootstrap with maximum likelihood
# # bs_pip4 <- bootstrap.pml(fit_pip4, bs=100, optNni=TRUE, multicore=TRUE)
# # bs_pip5 <- bootstrap.pml(fit_pip5, bs=100, optNni=TRUE, multicore=TRUE)
# 
# #bs_pip4 <- bootstrap.pml(fit4, bs=100, optNni=TRUE, multicore=TRUE)
# bs_pip5 <- bootstrap.pml(fit5, bs=100, optNni=TRUE, multicore=TRUE)
# # make starshaped layout trees
# #tree_ml_pip4 <- plotBS(fit_pip4$tree, bs_pip4)
# tree_ml_pip5 <- plotBS(fit_pip5$tree, bs_pip5)
# 
# # https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
# # also see : https://stackoverflow.com/questions/52956181/how-do-i-annotate-a-ggtree-object-with-bootstrap-values
# # gt04 <- ggtree::ggtree(tree_ml_pip4, right =T, #layout="slanted",
# #                        options(ignore.negative.edge=TRUE)) +
# #   ggtree::geom_tiplab(align=F, size=1.82) +
# #   ggtree::geom_nodelab(aes(label = label),color="red",hjust=-.3, size=4.2) +
# #   ggtree::geom_treescale() # adds the scale
# # gt04
# 
# # try again with only NEA samples
# gt05 <- ggtree::ggtree(tree_ml_pip5, right =T, #layout="slanted",
#                        options(ignore.negative.edge=TRUE)) +
#   ggtree::geom_tiplab(align=F, size=1.82) +
#   ggtree::geom_text2(aes(label=round(as.numeric(label)),
#                          subset=as.numeric(label)> 90,
#                          x=branch), vjust=0, color="red",hjust=0.3, size=3.2) +
#   ggtree::geom_treescale() # adds the scale
# gt05
# 
# # split label names by underscore to get a data frame
# df_pp5Nm <- data.frame(do.call
#                        ('rbind', 
#                          strsplit(as.character(names(phyd_pip5)),
#                                   "_")))
# # modify columns names
# colnames(df_pp5Nm) <- c("smplNm","loc","yer")
# # substitute to get location name
# NEAloc <-  gsub("(.*)_(.*)_(.*)","\\2",tree_ml_pip5$tip.label)
# # match with data frame that has colours for locations 
# cr <- df_cll$colfcol_loc[match(NEAloc,df_cll$coll_loc)]
# # specify the clade number for the long clade
# tr2 <- tidytree::groupClade(tree_ml_pip5, c(154))
# # plot a tree with colors on the long clade
# p5.01 <- ggtree(tr2, right =T, aes(color=group)) + 
#   ggtree::geom_tiplab(align=F, size=1.82) +
#   geom_tippoint(size = 2, color = factor(cr)) +
#   ggtree::geom_text2(aes(label=round(as.numeric(label)),
#                          subset=as.numeric(label)> 90, 
#                          x=branch), vjust=0, color="red",hjust=0.3, size=3.2) +
#   
#   scale_color_manual(values=c("black", "firebrick", "steelblue"))
# p5.01
# # plot a tree without colors on the long clade
# p5.02 <- ggtree(tr2, right =T, ladderize = T) + 
#   #theme_tree2() + 
#   ggtree::geom_tiplab(align=T, size=2.2, linesize=.5, hjust=-0.08) +
#   geom_tippoint(size = 2, color = factor(cr)) +
#   #ggtree::geom_text2(aes(label=node,x=branch), size=2.2) +
#   ggtree::geom_text2(aes(label=round(as.numeric(label)),
#                          subset=as.numeric(label)> 50, 
#                          x=branch), vjust=0, color="black",hjust=0.8, size=3.2) +
#   scale_color_manual(values=c("black", "firebrick", "steelblue"))
# #p5.02
# # subset the phydat object by the full names that match bein NEA samples
# df_pip5 <- df_pip4[match(NmslocNEA,rownames(df_pip4)),]
# # replace all NAs with '-'
# df_pip5.1 <- df_pip5 %>% replace(is.na(.), "-")
# #https://stackoverflow.com/questions/68925906/set-x-axis-on-ggtree-heatmap-in-r
# # make the data frame a matrix and then make it a DNAbin object
# dnb_pip5.02 <- as.DNAbin(as.matrix(df_pip5.1))
# #write the DNAbin object as a fasta file
# ape::write.FASTA(dnb_pip5.02,paste0(wd00_wd05,"/pip5_dnaseq.fasta"))
# # plot alignment together with tree, the offset defines the space between aligment and tree
# p5.03 <- ggtree::msaplot(p5.02, offset=0.04, fasta=paste0(wd00_wd05,"/pip5_dnaseq.fasta"))
# p5.03 <- p5.03 + scale_fill_manual(values=c("gray70","violetred3", "blue","orange", "green1"))
# p5.03 <- p5.03 + ggplot2::labs(fill='nt')
# #p5.03
# bSaveFigures=T
# # make an if test to check whether the plot should be saved
# # i.e. set to FALSE if you do not want the plot to be saved
# if(bSaveFigures==T){
#   ggsave(p5.03,file=paste0(wd00_wd05,"/Fig02_v03_MLtree_w_align.png"),
#          #width=210,height=297,
#          width=210,height=(297*0.5),
#          units="mm",dpi=300)
# }
# 
# #remove unneeded elements
# rm(bs_pip5,df_pip5,fit5,gt05,phyd_pip5)
#_______________________________________________________________________________
# start  -  make ML tree with alignment included for all the samples, 
# including the NCBI Genbank obtained samples
#_______________________________________________________________________________

#_______________________________________________________________________________
#_______________________________________________________________________________

(fit4 <- pml_bb(phyd_pip4, model="GTR+G"))
#(fit4 <- pml_bb(phyd_pip4, model="F81"))
#dna_dist <- dist.ml(mammals10, model="JC69")

class(tree_NJ4)
class(phyd_pip4)
fit_pip4 <- phangorn::pml(tree_NJ4, phyd_pip4)
bs_pip4 <- bootstrap.pml(fit4, bs=100, optNni=TRUE, multicore=TRUE)
tree_ml_pip4 <- plotBS(fit_pip4$tree, bs_pip4)

# https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
# also see : https://stackoverflow.com/questions/52956181/how-do-i-annotate-a-ggtree-object-with-bootstrap-values
gt04 <- ggtree::ggtree(tree_ml_pip4, right =T, #layout="slanted",
                       options(ignore.negative.edge=TRUE)) +
  ggtree::geom_tiplab(align=F, size=1.82) +
  ggtree::geom_nodelab(aes(label = label),color="red",hjust=-.3, size=4.2) +
  ggtree::geom_treescale() # adds the scale
gt04

# substitute to get location name
Nmsloc <-  gsub("(.*)_(.*)_(.*)","\\2",tree_ml_pip4$tip.label)
# match with data frame that has colours for locations 
cr <- df_cll$colfcol_loc[match(Nmsloc,df_cll$coll_loc)]

# specify the clade number for the long clade
tr4 <- tidytree::groupClade(tree_ml_pip4, c(154))

# I tried 'rescale' but this did not improve on differences in topology 
# branch lengths
# https://www.rdocumentation.org/packages/geiger/versions/2.0.10/topics/rescale
# plot a tree with colors on the long clade
p4.01 <- ggtree(tr4, right =T, aes(color=group)) + 
  ggtree::geom_tiplab(align=F, size=1.82) +
  geom_tippoint(size = 2, color = factor(cr)) +
  ggtree::geom_text2(aes(label=round(as.numeric(label)),
                         subset=as.numeric(label)> 90, 
                         x=branch), vjust=0, color="red",hjust=0.3, size=3.2) +
  
  scale_color_manual(values=c("black", "firebrick", "steelblue"))
p4.01
# plot a tree without colors on the long clade
p4.02 <- ggtree(tr4, right =T, ladderize = T) + 
  #theme_tree2() + 
  ggtree::geom_tiplab(align=T, size=1.2, linesize=.5, hjust=-0.08) +
  geom_tippoint(size = 2, color = factor(cr)) +
  #ggtree::geom_text2(aes(label=node,x=branch), size=2.2) +
  ggtree::geom_text2(aes(label=round(as.numeric(label)),
                         subset=as.numeric(label)> 50, 
                         x=branch), vjust=0, color="black",hjust=0.8, size=3.2) +
  scale_color_manual(values=c("black", "firebrick", "steelblue"))
#p4.02
# replace all NAs with '-'
df_pip4.1 <- df_pip4 %>% replace(is.na(.), "-")
#https://stackoverflow.com/questions/68925906/set-x-axis-on-ggtree-heatmap-in-r
# make the data frame a matrix and then make it a DNAbin object
dnb_pip4.02 <- as.DNAbin(as.matrix(df_pip4.1))
#write the DNAbin object as a fasta file
ape::write.FASTA(dnb_pip4.02,paste0(wd00_wd05,"/pip4_dnaseq.fasta"))
# plot alignment together with tree, the offset defines the space between aligment and tree
p4.03 <- ggtree::msaplot(p4.02, offset=0.05, fasta=paste0(wd00_wd05,"/pip4_dnaseq.fasta"))
p4.03 <- p4.03 + scale_fill_manual(values=c("gray70","violetred3", "blue","orange", "green1"))
p4.03 <- p4.03 + ggplot2::labs(fill='nt')

bSaveFigures=T
# make an if test to check whether the plot should be saved
# i.e. set to FALSE if you do not want the plot to be saved
if(bSaveFigures==T){
  ggsave(p4.03,file=paste0(wd00_wd05,"/Fig02_v04_MLtree_w_align.png"),
         #width=210,height=297,
         width=210,height=(297*0.5),
         units="mm",dpi=300)
}

#_______________________________________________________________________________
# end  -  make ML tree with alignment included for all the samples, 
# including the NCBI Genbank obtained samples
#_______________________________________________________________________________




# Other websites for reference
# # https://www.molecularecologist.com/2016/02/26/quick-and-dirty-tree-building-in-r/
# # https://wiki.duke.edu/display/AnthroTree/2.2+Simple+Parsimony+Analysis+in+R
# # https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html
## 2.4 Inferring Phylogeny using Maximum Likelihood in R (phangorn)
## https://wiki.duke.edu/pages/viewpage.action?pageId=131172124

# read in file with accession numbers to use for replacing Mnelei sample numbers
infl02 <- "list_of_Mnelei_numbers_in_table_to_replace_w_acc_nos.txt"
wd_NCBI.infl02 <- paste0(wd00,"/",
                         wd_NCBI,"/",
                         infl02)
df_acct01 <- read.csv(wd_NCBI.infl02, header=F,sep=",")
# rearrange the data frame to a longer format.
# use tidyr::pivot_longer
df_acct01 <- df_acct01 %>% tibble::rownames_to_column(var = "rowid") %>%
  tidyr::pivot_longer(-rowid, names_to = "colVNm",
                      values_to = "AbbrNm.smpl") %>%
  dplyr::group_by(AbbrNm.smpl)
# substitute all ' ' with nothing
df_acct01$AbbrNm.smpl <- gsub(" ","",df_acct01$AbbrNm.smpl)
# pick the Mnelei samples
Mnelei.smpls<- df_acct01$AbbrNm.smpl[grepl("Mnelei",df_acct01$AbbrNm.smpl)]
# make it a data frame
df_Mnelei.smpls <- as.data.frame(do.call(rbind, strsplit(Mnelei.smpls, 
                    split="-")), stringsAsFactors=FALSE)
# apply gsub to the entire data frame to get only numbers
# and make them numeric
mtx_smplNos <- apply(df_Mnelei.smpls, 2, function(y) as.numeric(gsub("Mnelei", "", y)))
df_smplNos <- as.data.frame(mtx_smplNos)
# make an empty column to add to
df_smplNos$Nw.spml.rng <- NA
df_smplNos$accno.rng <- NA
# iterrate over sample ranges
for (rn in (1:nrow(df_smplNos)))
{
  print(rn)
  rn <- 12
  # make a sequence of numbers to cover the range in the samples listed
  sqrng <- seq(df_smplNos$V1[rn],df_smplNos$V2[rn],1)
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  sqrng <- ifelse(nchar(sqrng)<3,stringr::str_pad(sqrng, 3, pad = "0"),
                  sqrng)
  # paste the species abbreviation onto the numbers
  smpl.rng <- paste0("Mnelei",sqrng)  
  # get the new accesion numbers
  Nw.accsNos<- df_accNos$AccNo[match(smpl.rng,df_accNos$MneleiAbbr)]
  Nw.accsNos <- Nw.accsNos[order(Nw.accsNos)]
  Nw.accsNos.1line <- paste(Nw.accsNos, collapse = ", ")
  smpl.rng.1line  <- paste(smpl.rng, collapse = ", ")
  df_smplNos$Nw.spml.rng[rn] <- smpl.rng.1line
  df_smplNos$accno.rng[rn] <- Nw.accsNos.1line
}
# combine to a new data frame
df_Mneleispml.w.accs <- cbind(df_Mnelei.smpls,df_smplNos,Mnelei.smpls)
# paste together a path and a filename
wd00_NCBI.outf <- paste0(wd00,"/",
                         wd_NCBI,
                         "/table_w_new_acc_nos_matched_w_Mnelei_smpl_nos.csv")
# write a csv file that has the data frame
write.table(df_Mneleispml.w.accs,
            file=wd00_NCBI.outf,
            col.names = T,
          sep=";")


