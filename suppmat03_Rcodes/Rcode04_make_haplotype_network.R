#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#remove everything in the working environment, without a warning!!
#rm(list=ls())


#remove everything in the working environment, without a warning!!
#rm(list=ls())
#define overall working directory
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
#define input working directories
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
#read in the data
inf01 <- "Mnelei_tmp01.fas"
inp.f.fnm <- "Mnelei_tmp01.fas"
wd00_wd05_inf01 <- paste0(wd00_wd05,"/",inf01)
dnb_fasf1 <- read.dna(file=wd00_wd05_inf01,format="fasta")
# split to get individual names
lbf <- strsplit(as.character(labels(dnb_fasf1)), "_")
#ge the different  sample numbers
smplno1 <- sapply(lbf, "[[", 1)
smplno2 <- sapply(lbf, "[[", 2)
smplno3 <- sapply(lbf, "[[", 3)
# combine to a table
tb_splNm1 <- table(smplno1, smplno2)
# extract haplotypes
hpt3 <- pegas::haplotype(dnb_fasf1)
#which polyps belong to which haplotype?
hpt3_indices<-attr(hpt3, "index")
# get number of labels 
hlab3 <- length(labels(hpt3))
#assign labels
names(hpt3_indices)<-c(1:hlab3)
# make roman numerals arabic numerals instead
labs.hpt3 <- as.numeric(as.roman(labels(hpt3)))
# use arabian numerals instead of roman numerals for haplotypes
hpt3 <- haplotype(dnb_fasf1,labels=c(labs.hpt3))
#make a haplonetwork
hpt3net <- pegas::haploNet(hpt3)
#mjn3net <- pegas::mjn(as.matrix(tb_splNm1))
#prepare hpt table
ind.hap3<-with(
  stack(setNames(attr(hpt3, "index"), rownames(hpt3))),
  table(hap=ind, pop=rownames(dnb_fasf1)[values]))
#make it a dataframe
df_ihpt03 <- as.data.frame(ind.hap3)
# make haplotype labels arabic numerals instead of roman numerals 
df_ihpt03$hap.arab <- as.numeric(as.roman(as.character(df_ihpt03$hap)))
# split the string with sequence name which holds the location name
lbf3  <- strsplit(as.character(df_ihpt03$pop), "_")
# get the second element from this string split
df_ihpt03$pop.loc <- sapply(lbf3, "[[", 2)
#remove any zero occurences
df_ihpt03 <- df_ihpt03[df_ihpt03$Freq >= 1,]	
# make a table of locations and haplotype numbers
hsl3 <- table(df_ihpt03$hap.arab, df_ihpt03$pop.loc)

#plot the network
plot(hpt3net, size = (sqrt(attr(hpt3net,"freq")/pi)), 
     scale.ratio = 0.2, cex = 0.5, pie = hsl3, 
     show.mutation = 0, threshold = 0) #, labels(TRUE))

# 


flnm <- c(paste("Fig02_v03_haplotype_network_",inp.f.fnm,"02.jpg",  sep = ""))
#paste output directory and filename together in a string
outflnm <- paste(wd00_wd05,"/",flnm,sep="")
# Exporting PFD files via postscript()           
# pdf(outflnm,
#     width=(1*1.0*8.2677),height=(4*1.0*2.9232))
jpeg(outflnm,
     width=(3200),height=(4800),res=300)
#define lpot arrangement
tbt.par <- par(mfrow=c(1, 1),
               oma=c(0,0,0,0), #define outer margins
               mai=c(0,0,0,0), #define inner margins
               mar=c(0,0,2,0))

#plot the network
plot(hpt3net, size = (sqrt(attr(hpt3net,"freq")/pi)), 
     scale.ratio = 0.2, cex = 0.5, pie = hsl3, 
     show.mutation = 0, threshold = 0) #, labels(TRUE))
#add a legend to the plot
legend("bottomleft",unique(df_ihpt03$pop.loc), 
       #pt.bg=colfh,
       box.col=NA,
       #col=rainbow(ncol(new.hap.smplloc)), 
       pt.lwd=0.4,
       pch=21, ncol=1, cex=0.8)
title(main = "a",
      cex.main = 1.8,   font.main= 2, col.main= "black",
      adj = 0.01, line = 0.1)

#_______________________________________________________________________________
par(tbt.par)
# end svg file to save as
dev.off()  
#reset this parameter
par(mfrow = c(1, 1)) 






#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
# 
# 
# 
# 
# 
# 
# # make a haplotype object
# h2 <- haplotype(dnb_x2)
# hlab2 <- length(labels(h2))
# h2 <- haplotype(dnb_x2,labels=c(1:hlab2))#default labels are roman numerals. Vector of names has to be of equal length as the number of haplotypes.
# #h_selection<-haplotype(x_selection)
# #which polyps belong to which haplotype?
# h2_indices<-attr(h2, "index")
# #assign labels
# names(h2_indices)<-c(1:hlab2)
# #str(h_indices)
# #generate plain network 
# net2 <- haploNet(h2)
# #get list of all sequence labels with haplotype numbers
# haplotypes2<-stack(setNames(attr(h2, "index"), rownames(h2)))
# 
# #haploNet
# 
# 
# plot(net2,
#      size=(sqrt((attr(net2, "freq"))/pi)),
#      scale.ratio = 0.3,
#      cex = 1.2,
#      bg= colfh, 
#      pie = new.hap.smplloc2 , 
#      show.mutation=2,
#      threshold=0)
# 
# 
# rownames(dnb_x2)
# #prepare hpt table
# ind.hap2<-with(
#   stack(setNames(attr(h2, "index"), rownames(h2))),
#   table(hap=ind, pop=rownames(dnb_x2)[values]))
# #make it a dataframe
# df_ihpt02 <- as.data.frame(ind.hap2)
# #limit to include only 'Freq' that equals 1
# df_ihpt02 <- df_ihpt01[df_ihpt01$Freq == 1,]	
# 
# 
# 
# 
# 
# #names(getHaploNetOptions())
# #extract all unique haplotypes and attach custom labels
# h <- haplotype(x)
# hlab <- length(labels(h))
# h <- haplotype(x,labels=c(1:hlab))#default labels are roman numerals. Vector of names has to be of equal length as the number of haplotypes.
# #h_selection<-haplotype(x_selection)
# #which polyps belong to which haplotype?
# h_indices<-attr(h, "index")
# #assign labels
# names(h_indices)<-c(1:hlab)
# #str(h_indices)
# #generate plain network 
# net <- haploNet(h)
# #haploNet
# #net <- haploNet(h_selection)
# #plot
# plot(net,
#      size=attr(net, "freq"),
#      scale.ratio = 3,
#      cex = 0.8,
#      show.mutation=2,
#      threshold=0)
# 
# plot(net,
#      size=(sqrt((attr(net, "freq"))/pi)),
#      scale.ratio = 0.85,
#      cex = 1.1,
#      show.mutation=2,
#      threshold=0)
# 
# #extract all names
# names<-data.frame(labels(x))
# 
# #write.table(names,file="names_def.csv",sep=",")
# haploFreq(x,haplo=h)
# #get list of all sequence labels with haplotype numbers
# haplotypes<-stack(setNames(attr(h, "index"), rownames(h)))
# #Grep among the list of files, to get the file that holds the different names
# fnm_clo3 <- list.files(wd00_wd05)[grepl("clo03",list.files(wd00_wd05))]
# fnm_df_lN02 <- list.files(wd00_wd05)[grepl("df_lN02",list.files(wd00_wd05))]
# #read locations file
# df_clo3<-read.csv(file=paste0(wd00_wd05,"/",fnm_clo3),header=TRUE,sep=";")
# df_lN02<-read.csv(file=paste0(wd00_wd05,"/",fnm_df_lN02),header=TRUE,sep=",")
# # split string and get lists nested in a list
# lbls01 <- strsplit(as.character(labels(x)), "_")
# # get only NCBI labels
# lbls02 <- sapply(lbls01, "[[", 3)
# # grep only elements that includes numbers
# ncbilbls <- lbls02[grepl("[0-9]",lbls02)]
# # get location names for NCBI accesion numbers
# lnbls<- df_lN02$location[match(ncbilbls,df_lN02$accession_nmb)]
# # replace the labels that has numbers in them with the location name
# lbls02[grepl("[0-9]",lbls02)] <- lnbls
# # copy the vector with labels
# lbls03 <- lbls02
# df_lblnms <- as.data.frame(cbind(lbls02,lbls03)) 
# colnames(df_lblnms) <- c("longNm","AbrNm")
# df_lblnms$AbrNm <- df_clo3$locality8[match(df_lblnms$longNm,df_clo3$locality5)]
# llnNm1 <- df_lblnms$longNm[is.na(df_lblnms$AbrNm)]
# # make vector with long names appearing in the haplotype labels
# Nm_inLonNm <- c("USA:WoodsHole","NEAtlantic","Mediterranean","CaspianSea","CentralWAtlantic","USA:GalvestonBay","USA:Panacea","Germany:KielFjord","Germany:Maasholm","Germany:Helgoland","BalticSea","AtlanticOcean:NWAtlantic","Netherlands","Loegstoer","Mariagerfjord","Kerteminde","NSeaHelgolandRds","WaddenSeaBussumHaupstr","Skovshoved","KielFjord","Bogense","MecklenburgerBuchtWismarBucht","MecklenburgerBucht","Ballen")
# # make vector with long names appearing in the df_clo3 df
# Nm_inclo03 <- c("USA:WoodsHole","NEAtlantic","Mediterranean","CaspianSea","CentralWAtlantic","USA:GalvestonBay","USA:Panacea","NGermanyKielFjord","NGermanyKielFjord","NWGermanyNSeaHelgolandRds","BalticSea","AtlanticOcean:NWAtlantic","Netherlands","NJyllandLimfjord","JyllandMariagerfjord","FynKerteminde","NWGermanyNSeaHelgolandRds","GermanyBusum","SjaellandSkovshoved","NGermanyKielFjord","FynBogense","NGermanyMecklenburgerBuchtWismarBucht","NGermanyMecklenburgerBuchtWismarBucht","SamsoeBallen")
# # combine to a df
# df_rplNm <- as.data.frame(cbind(Nm_inLonNm,Nm_inclo03))
# # match first the longNm from the Haplotype labels to the df with both the haplotype labels and the df_clo03 labels
# df_lblnms$longNm <- df_rplNm$Nm_inclo03[match(df_lblnms$longNm,df_rplNm$Nm_inLonNm)]
# # then match the new long name with th df_clo03 table to get the abbreviation
# df_lblnms$AbrNm <- df_clo3$locality8[match(df_lblnms$longNm,df_clo3$locality6)]
# #add rownames for haploFreq function
# rownames(x2) <- df_lblnms$AbrNm
# #The following code is used to plot the haplotype and add the haplotype frequencies. 
# #create haplotypes from dna.bin
# #x
# px <- pegas::haplotype(x2)
# #prepare hpt table
# ih<-with(
#   stack(setNames(attr(h, "index"), rownames(h))),
#   table(hap=ind, pop=rownames(x2)[values]))
# #make it a dataframe
# ih01 <- as.data.frame(ih)
# #limit to include only 'Freq' that equals 1
# ih02 <- ih01[ih01$Freq >= 1,]	
# hpfrq <- haploFreq(x2,rownames(x2),haplo=h)
# #make it a table
# hsl <- table(ih02$hap, ih02$pop)
# #per area THIS IS THE ONE USED IN THE PAPER
# #freq_area<-haploFreq(x,split="_",what=2,haplo=h)
# #plot the new haplotype
# 
# plot(net,
#      #size=(sqrt(attr(net, "freq")/pi)),
#      size=attr(net, "freq"),
#      scale.ratio = 2.7,
#      cex = 2.8,
#      pie=ih, #hsl
#      #pie=ih02,
#      show.mutation=2,# mutations are shown with small dots on the links
#      threshold=0#turn off display of alternative links
#      )
# 
# legend(30,20, colnames(hsl), col=rainbow(ncol(hsl)), pch=20,bty="n")
# 
# 
# 
# plot(net,
#      size=(sqrt((attr(net, "freq"))/pi)),
#      scale.ratio = 0.85,
#      cex = 1.1,
#      #pie=ih,
#      show.mutation=2,
#      threshold=0)
# #https://cran.r-project.org/web/packages/pegas/vignettes/PlotHaploNet.pdf
# #_______________________________________________________________________________
# # Advice from Cody Aylward:  - https://doi.org/10.1007/s10592-018-1130-3
# 
# # Hello Steen,
# # 
# # Thanks for your comments on my paper. 
# #I looked back at my old hard drive and unfortunately cannot 
# #  find the exact R script used for that figure. However, I have a few suggestions. 
# #I read through this vignette to refresh my memory on constructing haplotype 
# #networks in Pegas (https://cran.r-project.org/web/packages/pegas/vignettes/PlotHaploNet.pdf), 
# #I found it quite helpful.
# # 
# # My first suggestion is to re-size of your circles. 
# #It looks like you have one haplotype that is substantially more common than the rest, 
# #such that it is covering other circles completely. 
# #I am assuming you made this figure using the default setting, 
# #where the diameter of each circle is equal to the frequency of each haplotype. 
# #An alternative is to set the diameter of each circle equal to the square root
# #of the frequency of each haplotype divided by pi 
# #(such that the area of each circle is equal to the frequency of 
# #each haplotype). This is explained on pages 11-12 of the vignette.
# #Another alternative is to manually size the circles. 
# #This was easier for my paper since I only had 12 haplotypes, 
# #and might be a bit more of a task for your larger number of haplotypes. 
# #However, if the square root approach doesn't shrink your largest circle 
# #sufficiently, you might consider manually creating a vector 
# #for your circle sizes. If you use this approach, 
# #I would group haplotype frequencies into different size categories. 
# #For example, haplotype frequencies of 1-5 would be size 1, haplotype 
# #frequencies 6-20 are size 2, etc. 
# #That may be a bit more work, but it gives you more direct manipulation 
# #of the network. Page 6 of the vignette shows how to order the size
# #vector to match the order of your haplotypes in the network.
# # 
# # I think once the circles are resized your network will be much cleaner. 
# #If that isn't sufficient, there may be a way to manually adjust 
# #the coordinates of each circle. I did a very quick search and didn't 
# #find an answer right away. If you still need to modify the network after 
# #re-sizing the circles I would be happy to take a longer look at 
# #manually adjusting coordinates.
# # 
# # Best of luck with your study,
# # 
# # Cody Aylward, M.S.
# #______________________________________________________________________________
# 
# # vignette try out
# 
# library(pegas) # loads also ape
# library(ape) 
# #detach("package:genetics", unload=TRUE)
# #detach("package:ape", unload=TRUE)
# #install.packages("ape")
# #install.packages("pegas")
# data(woodmouse)
# 
# set.seed(10)
# x <- woodmouse[sample.int(nrow(woodmouse), 80, TRUE), ]
# region <- rep(c("regA", "regB"), each = 40)
# pop <- rep(paste0("pop", 1:4), each = 20)
# table(region, pop)
# # extract the haplotypes which are used to reconstruct the RMST after computing the
# # pairwise Hamming distances:
# h <- pegas::haplotype(x)
# #hDb <- as.DNAbin(h)
# d <- ape::dist.dna(h, "N")
# 
# d <- ape::dist.dna(x)
# nt <- rmst(d, quiet = TRUE)
# nt
# plot(nt)
# plot(nt, fast = TRUE)
# # By default, not all links are drawn. This is controlled with the option threshold which
# # takes two values in order to set the lower and upper bounds of the number of mutations for
# # a link to be drawn:
# plot(nt, threshold = c(1, 14))
# # The visual aspect of the links is arbitrary: the links of the backbone MST are shown with
# # continuous segments, while “alternative” links are shown with dashed segments.
# args(pegas:::plot.haploNet)
# # Like for most plot methods, the first argument (x) is the object to be plotted. Until
# # pegas 0.14, all other arguments were defined with default values. In recent versions, as
# # shown above, only size and shape are defined with default values; the other options, if not
# # modified in the call to plot, are taken from a set of parameters which can be modified as
# # explained in Section 3.4.
# # The motivation for this new function definition is that in most cases users need to modify
# # size and shape with their own data, such as haplotype frequencies or else, and these might
# # be changed repeatedly (e.g., with different data sets or subsets). On the other hand, the
# # other options are more likely to be used to modify the visual aspect of the graph, so it could
# # be more useful to change them once during a session as explained later in this document.
# # The size of the haplotype symbols can be used to display haplotype frequencies. The
# # function summary can extract these frequencies from the "haplotype" object:
# 
# (sz <- summary(h))
# # It is likely that these values are not ordered in the same way than haplotypes are ordered
# # in the network
# (nt.labs <- attr(nt, "labels"))
# #It is simple to reorder the frequencies before using them into plot
# sz <- sz[nt.labs]
# plot(nt, size = sz)
# (R <- haploFreq(x, fac = region, haplo = h))
# # A  similar mechanism can be used to show variables such as region or pop. The function
# # haploFreq is useful here because it computes the frequencies of haplotypes for each region
# # or population:
# (P <- haploFreq(x, fac = pop, haplo = h))
# # Like with size, we have to reorder these matrices so that their rows are in the same order
# # than in the network:
# R <- R[nt.labs, ]
# P <- P[nt.labs, ]
# # We may now plot the network with either information on haplotype frequencies by just
# # changing the argument pie:
# plot(nt, size = sz, pie = R, legend = c(-25, 30))
# plot(nt, size = sz, pie = P, legend = c(-25, 30))
# 
# # The option legend can be:
# # FALSE (the default): no legend is shown;
# # TRUE: the user is asked to click where the legend should be printed;
# # a vector of two values with the coordinates where the print the legend (for non-
# # interactive use like in this vignette
# # New Features in pegas 1.0
# # This section details some of the improvements made to haplotype network drawing after
# # pegas 0.14.
# # 3.1 Improved ‘Replotting’
# # The graphical display of networks is a notoriously difficult problem, especially when there is
# # an undefined number of links (or edges). The occurrence of reticulations makes line crossings
# # almost inevitable. The packages igraph and network have algorithms to optimise the layouts
# # of nodes and edges when plotting such networks.
# #
# # The function replot (introduced in pegas 0.7, March 2015) lets the user modify the
# # layout of nodes interactively by clicking on the graphical window where a network has been
# # plotted beforehand. replot—which cannot be used in this non-interactive vignette—has
# # been improved substantially:
# # The explanations printed when the function is called are more detailed and the node
# # to be moved is visually identified after clicking.
# # The final coordinates, for instance saved with xy <- replot(), can be used directly
# # into plot(nt, xy = xy). This also makes possible to input coordinates calculated
# # with another software.
# # In previous versions, the limits of the plot tended to drift when increasing the number
# # of node moves. This has been fixed, and the network is correctly displayed whatever
# # the number of moves done
# 
# # Haplotype Symbol Shapes
# # Haplotypes can be represented with three different shapes: circles, squares, or diamonds.
# # The argument shape of plot.haploNet is used in the same way than size as explained
# # above (including the evental need to reorder the values). Some details are given below
# # about how these symbols are scaled.
# # There are two ways to display a quantitative variable using the size of a circle: either
# # with its radius (r) or with the area of the disc defined by the circle. This area is πr2, so if we
# # want the area of the symbols to be proportional to size, we should square-root these last
# # values. However, in practice this masks variation if most values in size are not very different
# # (see below). In pegas, the diameters of the circles (2r) are equal to the values given by size.
# # If these are very heterogeneous, they could be transformed with size = sqrt(.... keeping
# # in mind that the legend will be relative to this new scale.
# # The next figure shows both ways of scaling the size of the circles: the top one is the
# # scaling used in pegas.
# 
# par(xpd = TRUE)
# size <- c(1, 3, 5, 10)
# x <- c(0, 5, 10, 20)
# plot(0, 0, type="n", xlim=c(-2, 30), asp=1, bty="n", ann=FALSE)
# other.args <- list(y = -5, inches = FALSE, add = TRUE,bg = rgb(1, 1, 0, .3))
# o <- mapply(symbols, x = x, circles = sqrt(size / pi),MoreArgs = other.args)
# other.args$y <- 5
# o <- mapply(symbols, x = x, circles = size / 2,MoreArgs = other.args)
# text(x, -1, paste("size =", size), font = 2, col = "blue")
# text(30, -5, expression("circles = "*sqrt(size / pi)))
# text(30, 5, "circles = size / 2")
# #
# # For squares and diamonds (shape = "s" and shape = "d", respectively), they are scaled
# # so that their areas are equal to the disc areas for the same values given to size. The figure
# # below shows these three symbol shapes superposed for several values of this parameter. Note
# # that a diamond is a square rotated 45° around its center
# 
# x  <- c(0, 6, 13, 25)
# plot(0, 0, type="n", xlim=c(-2, 30), asp=1, bty="n", ann=FALSE)
# other.args$y <- 0
# o <- mapply(symbols, x = x, circles = size/2, MoreArgs = other.args)
# other.args$col <- "black"
# other.args$add <- other.args$inches <- NULL
# o <- mapply(pegas:::square, x = x, size = size, MoreArgs = other.args)
# o <- mapply(pegas:::diamond, x = x, size = size, MoreArgs = other.args)
# text(x, -7, paste("size =", size), font = 2, col = "blue")
# # The Function mutations
# # mutations() is a low-level plotting function which displays information about the mutations
# # related to a particular link of the network. This function can be used interactively. For
# # instance, the following is copied from an interactive R session:
# #mutations(nt)
# 
# plot(nt)
# mutations(nt, 18, x = -8.9, y = 16.3, data = h)
# #Like any low-level plotting function, mutations() can be called as many times
# #as needed to display similar information on other links.
# #The option style takes the value "table"
# # (the default) or "sequence". In the second, the positions of the mutations are drawn on a
# # horizontal segment representing the sequence:
# plot(nt)
# mutations(nt, 18, x = -8.9, y = 16.3, data = h)
# mutations(nt, 18, x = 10, y = 17, data = h, style = "s")
# #
# # Getting and Setting Options
# # The new version of pegas has two ways to change some of the parameters of the plot:
# # either by changing the appropriate option(s) in one of the above functions, or by
# # setting these values with the function setHaploNetOptions, in which case all subsequent
# # plots will be affected.2 The list of the option values currently in use can be printed
# # with getHaploNetOptions. There is a relatively large number of options that affect
# # either plot.haploNet() or mutations(). Their names are quite explicit so that the user
# # should find which one(s) to modify easily:
# # names(getHaploNetOptions())
# # # We see here several examples with the command plot(nt, size = 2) which is repeated
# # # after calling setHaploNetOptions
# # plot(nt, size = 2)
# # 
# # setHaploNetOptions(haplotype.inner.color = "#CCCC4D",
# #                    haplotype.outer.color = "#CCCC4D",
# #                    show.mutation = 3, labels = FALSE)
# # plot(nt, size = 2 )
# # setHaploNetOptions(haplotype.inner.color = "blue",
# #                    haplotype.outer.color = "blue", show.mutation = 1)
# # par(bg = "yellow3")
# # plot(nt, size = 2)
# # 
# # setHaploNetOptions(haplotype.inner.color = "navy",
# #                    haplotype.outer.color = "navy")
# # par(bg = "lightblue")
# # plot(nt, size = 2)
# # par(bg = "white")
# #______________________________________________________________________________
# # end vignette
# #______________________________________________________________________________
# 
# # get summary  for haplotype object
# (szp <- summary(h))
# 
# # It is likely that these values are not ordered in the same way than haplotypes are ordered
# # in the network
# (nt.labsp <- attr(net, "labels"))
# #It is simple to reorder the frequencies before using them into plot
# szp <- szp[as.numeric(nt.labsp)]
# plot(net, size = sqrt(szp / pi))
# regions2 <- rownames(x2)
# (Rp <- haploFreq(x2, fac = regions2, haplo = h))
# # A  similar mechanism can be used to show variables such as region or pop. The function
# # haploFreq is useful here because it computes the frequencies of haplotypes for each region
# # or population:
# #(P <- haploFreq(x, fac = pop, haplo = h))
# # Like with size, we have to reorder these matrices so that their rows are in the same order
# # than in the network:
# Rp <- Rp[as.numeric(nt.labsp), ]
# 
# # We may now plot the network with either information on haplotype frequencies by just
# # changing the argument pie:
# plot(net, size = sqrt(szp / pi), 
#      scale.ratio = 1.8,
#      cex = 1.8,
#      threshold=0,#turn off display of alternative links
#      show.mutation=2,# mutations are shown with small dots on the links
#      pie = Rp, 
#      legend = c(-25, 30))
# 
# 
# hsl2 <- sqrt(attr(net, "freq") / pi)
# 
# plot(net,
#      size=sqrt(attr(net, "freq") / pi),
#      scale.ratio = 2,
#      cex = 0.8,
#      pie=hsl,
#      show.mutation=2,# mutations are shown with small dots on the links
#      threshold=0#turn off display of alternative links
# )
# 
# 
# plot(net, size = sqrt(szp / pi), 
#      scale.ratio = 1.8,
#      cex = 1.8,
#      threshold=0,#turn off display of alternative links
#      show.mutation=2,# mutations are shown with small dots on the links
#      pie = Rp, 
#      legend = c(-25, 30))
