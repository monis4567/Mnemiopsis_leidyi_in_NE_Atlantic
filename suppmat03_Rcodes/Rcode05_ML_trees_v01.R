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


# calculate distances
dm4 <- phangorn::dist.ml(phyd_pip4, "F81")
dm5 <- phangorn::dist.ml(phyd_pip5, "F81")
tree_NJ4 <- NJ(dm4)
tree_NJ5 <- NJ(dm5)
# as alternative for a starting tree:
tree4 <- pratchet(phyd_pip4)          # parsimony tree
tree5 <- pratchet(phyd_pip5)          # parsimony tree
tree4 <- nnls.phylo(tree4, dm4)   # need edge weights
tree5 <- nnls.phylo(tree5, dm5)   # need edge weights
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
fit_pip4 <- pml(tree_NJ4, phyd_pip4)
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

