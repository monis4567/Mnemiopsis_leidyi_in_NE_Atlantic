#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-


#remove everything in the working environment, without a warning!!
#rm(list=ls())

#_______________________________________________________________________________
# make nucleotide diversity tables - start
#_______________________________________________________________________________
# you will need 'gls' see:
# https://stackoverflow.com/questions/69583193/installing-r-gsl-package-on-ubuntu
# In the terminal use these commands
# $ sudo add-apt-repository ppa:dns/gnu
# $ sudo apt-get update
# $ sudo apt install libgsl-dev
# # Then install 'gsl'
# install.packages('gsl')
# # and then install the other bits required for 'strataG's
# install.packages('copula')
# # to get 'sprex' installed you need "rjags"
# # in the termnial install 'rjags' like this:
# # $ sudo apt install jags
# #Then, in R
# install.packages("rjags")
# install.packages('swfscMisc')
# install.packages('Matrix')
# install.packages('units')
library('units')
library('Matrix')
library('sprex')
#
# #https://www.rdocumentation.org/packages/strataG/versions/2.4.905
# if(!require(strataG)){
#   # make sure you have Rtools installed
#   if (!require('devtools')) install.packages('devtools')
#   # install sprex development version
#   devtools::install_github("ericarcher/sprex", force = T)
#   # install from GitHub
#   devtools::install_github('ericarcher/strataG', build_vignettes = TRUE, force = T)
#   #install.packages('strataG')
# }
# get the library for the strataG package
library(strataG)
library(bioseq)
library(ape)
library(dplyr)
library(pegas)
library(htmlTable)
# define the directory with the NCBI accession numbers
wd_ncbi <- "NCBI_seq_submission_for_Mnemiopsis"
# define the file with the NCBI accession numbers
inf_ncbi_accNo <- "GenBank_accession_numbers_for_the_Mnelei_submitted_sequences_2024apr17.txt"
# read the file with the accession numbers
pthinf_ncbi_accNo <- paste(wd00,wd_ncbi,inf_ncbi_accNo, sep="/")
dlm_fl_ncbi  <- read.delim(pthinf_ncbi_accNo)
# grep for "Mnelei" in the 1 st column
df_fl_ncbi <- dlm_fl_ncbi[grepl("Mnelei",dlm_fl_ncbi[,1]),]
# substitute a character that occurs more than 1 time
# https://stackoverflow.com/questions/9562535/gsub-reduce-all-repeating-characters-to-one-instance
df_fl_ncbi <- gsub("([ ])\\1+"," ",df_fl_ncbi)
#split by delimiter
#https://stackoverflow.com/questions/7069076/split-column-at-delimiter-in-data-frame
df_no_ncbi <- data.frame(do.call('rbind', strsplit(as.character(df_fl_ncbi),' ',fixed=TRUE)))
# alter the column names
colnames(df_no_ncbi) <- c("SUBNo",
                          "MneleiNo",
                          "NCBIAccNo") 



#combine the alignments in a list
lst_aM.ITS <- list(alvrs.M.ITS1, alvrs.M.ITS2)
# also make a vector with the gene names
gene.nms <- c("ITS1", "ITS2")
# count the number of elements in the list with alignments
n.alvrs.M <- length(lst_alvrs.M.ITS)
# make an empty list to collect frequency tables in
lst.frq.tbls <- list()
lst.nucdiv.tbls <- list()
#make a sequence of numbers that can reflect the genes for alignmetns
n.f.alvrs.M <- seq(1, n.alvrs.M, 1)
#iterate over the  sequence of numbers that  reflect the genes
for (ng in n.f.alvrs.M)
{
  print(ng)
  
  #}
  # get the corresponding element from the list of alignments
  alvrsM <- lst_alvrs.M.ITS[[ng]]
  #get the gene name
  GnNm <- gene.nms[ng]
  # Now the alignment has been picked from the list,
  # Then prepare the NJ trees, and the label categories
  # and a color range to assign the labels in the trees
  # also get the labels on the sequences
  LalvrsM <- names(alvrsM)
  alvrsM <- bioseq::as_DNAbin(alvrsM)
  names(alvrsM) <- LalvrsM
  #
  dnb_pip4 <-  alvrsM
  # The DNAbin object MUST have rownames in order to be able to use the
  # 'pegas::haploFreq' function
  # Transforming the DNAbin object
  # To a genind object  and then to a data frame and then to a matrix,
  # and then back into a a DNAbin object allows for ensuring the rownames
  # are present in the DNAbin object
  #make the DNAbin object a genind object
  geni_pip <- adegenet::DNAbin2genind(dnb_pip4)
  #make the genind object a dataframe
  df_pip5 <- adegenet::genind2df(geni_pip)
  
  #make the df object a matrix and make this a dnabin object
  dnb_pip05 <- ape::as.DNAbin(as.matrix(df_pip5))
  # checj there are rownames present in the DNAbin object
  is.vector(rownames(dnb_pip4))
  is.vector(rownames(dnb_pip05))
  # split the row name string by a character
  rNmp6 <- strsplit(as.character(LalvrsM), "_")
  # get the 2nd element per vector in the list - this holds the location
  lcNm9 <- sapply(rNmp6, "[[", 2)
  # count up frequency of samples
  df_popinf01 <- as.data.frame(table(lcNm9))
  # make the column a character column
  df_popinf01$lcNm9 <- as.character(df_popinf01$lcNm9)
  # make a vector with total counts
  vTtcnt <- c("Total", sum(df_popinf01$Freq))
  #  append a row with total counts
  df_popinf01 <- rbind(df_popinf01, vTtcnt)
  # change column names
  colnames(df_popinf01) <- c("popNm", "smplSz")
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
    if (popl == "Total") {
      popl2 <- "*"
    } else {
      popl2 <- popl
    }
    # subset the data frame to only comprise the current pop name
    df_psb09 <- df_pip5[grepl(popl2, LalvrsM), ]
    # make the data frame a matrix and make it a DNAbin object
    dnb_ps09 <- ape::as.DNAbin(as.matrix(df_psb09))
    # get the nucleotide diversity and the haplotype diversity
    nc.dv      <- pegas::nuc.div(dnb_ps09, variance = F)
    # setting 'variance = T' returns the standard deviation for the subsamples
    nc.dv.sd   <- pegas::nuc.div(dnb_ps09, variance = T)[2]
    hc.dv      <- pegas::hap.div(dnb_ps09, variance = F)
    hc.dv.sd   <- pegas::hap.div(dnb_ps09, variance = T)[2]
    # make the data frame a matrix
    mtx_ps09 <- as.matrix(df_psb09)
    # split matrix using row names
    lstv_ps09 <- split(mtx_ps09, f = row.names(mtx_ps09))
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
  df_ndv <- as.data.frame(do.call(rbind, lst_nd))
  df_hdv <- as.data.frame(do.call(rbind, lst_hd))
  df_ndv.sd <- as.data.frame(do.call(rbind, lst_nd.sd))
  df_hdv.sd <- as.data.frame(do.call(rbind, lst_hd.sd))
  df_tD <- as.data.frame(do.call(rbind, lst_tD))
  df_tD.pvln <- as.data.frame(do.call(rbind, lst_tD.pvln))
  df_tD.pvlb <- as.data.frame(do.call(rbind, lst_tD.pvlb))
  df_FusF <- as.data.frame(do.call(rbind, lst_FusF))
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
  colnames(df_ndv)       <- c("ndv", "poplcNm")
  colnames(df_hdv)       <- c("hdv", "poplcNm")
  colnames(df_ndv.sd)    <- c("ndv.sd", "poplcNm")
  colnames(df_hdv.sd)    <- c("hdv.sd", "poplcNm")
  colnames(df_tD)    <- c("tD", "poplcNm")
  colnames(df_tD.pvln)    <- c("tD.pvln", "poplcNm")
  colnames(df_tD.pvlb)    <- c("tD.pvlb", "poplcNm")
  colnames(df_FusF)    <- c("FusF", "poplcNm")
  #combine data frames in a list
  lst_df.nhdv <- list(df_ndv,
                      df_hdv,
                      df_ndv.sd,
                      df_hdv.sd,
                      df_tD,
                      df_tD.pvln,
                      df_tD.pvlb,
                      df_FusF)
  #merge all data frames in list
  df_nhdv <- lst_df.nhdv %>% purrr::reduce(full_join, by = 'poplcNm')
  #3
  hpt.pip <- pegas::haplotype(dnb_pip05)
  # use the DNAbin object to get the frequency of haplotypes per population
  # this fucntion must have the rownames in the DNAbin object, otherwise,
  # it cannot work.
  mtx_HptFrq <- pegas::haploFreq(dnb_pip05,
                                 split = "_",
                                 what = 2,
                                 haplo = hpt.pip)
  # set all hapltype counts to 1 , to be able to sum up unique haplotypes
  mtx_HptFrq[mtx_HptFrq >= 1] <- 1
  # get the haplotype frequency - i.e the number of haplotypes per location
  hptF.p.loc <- colSums(mtx_HptFrq)
  # get the haplotype frequency - i.e the number of haplotypes per haplotype
  hptF.p.hpt <- rowSums(mtx_HptFrq)
  # only count occurences of haplotypes
  hptF.p.hpt[hptF.p.hpt >= 1] <- 1
  # count the total number of haplotypes
  tot.hpt.cnt <- sum(hptF.p.hpt)
  # make it a data frame
  df_hptF.p.loc <- as.data.frame(hptF.p.loc)
  # get the rownames to make a column
  df_hptF.p.loc$popNm  <- row.names(df_hptF.p.loc)
  # also add the total count of haplotypes
  df_hptF.p.loc <- rbind(c(tot.hpt.cnt, "Total"), df_hptF.p.loc)
  # match back to get the number of Haplotypes per population
  df_popinf01$Nh <- df_hptF.p.loc$hptF.p.loc[match(df_popinf01$popNm, df_hptF.p.loc$popNm)]
  # match to get haplotype diverstiy and nucleotide diversity and standard deviations
  df_popinf01$hdv       <- df_nhdv$hdv[match(df_popinf01$popNm, df_nhdv$poplcNm)]
  df_popinf01$hdv.sd    <- df_nhdv$hdv.sd[match(df_popinf01$popNm, df_nhdv$poplcNm)]
  df_popinf01$ndv       <- df_nhdv$ndv[match(df_popinf01$popNm, df_nhdv$poplcNm)]
  df_popinf01$ndv.sd    <- df_nhdv$ndv.sd[match(df_popinf01$popNm, df_nhdv$poplcNm)]
  df_popinf01$tD    <- df_nhdv$tD[match(df_popinf01$popNm, df_nhdv$poplcNm)]
  df_popinf01$tD.pvln    <- df_nhdv$tD.pvln[match(df_popinf01$popNm, df_nhdv$poplcNm)]
  df_popinf01$tD.pvlb    <- df_nhdv$tD.pvlb[match(df_popinf01$popNm, df_nhdv$poplcNm)]
  df_popinf01$tD.pvlb    <- df_nhdv$tD.pvlb[match(df_popinf01$popNm, df_nhdv$poplcNm)]
  df_popinf01$FusF    <- df_nhdv$FusF[match(df_popinf01$popNm, df_nhdv$poplcNm)]
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
  df_popinf02 <- as.data.frame(lapply(df_popinf02, function(x)
    as.numeric(as.character(x))))
  df_popinf02[is.na(df_popinf02)] <- 0
  
  # use a previously defined data frame with long names and four letter
  # abbreviations for locations
  df_4LtabbrNm <- df_4LtabbrNm[order(df_4LtabbrNm$llabbr), ]
  # define a vector with the long names locations you wnat to retain
  lngNmtk <- c(
    "BalticSea",
    "CaspianSea",
    "CentralWAtlantic",
    "Funen, Bogense",
    "Funen, Kerteminde",
    "Germany, Büsum",
    "North Sea, Helgoland Roads",
    "Germany, Kiel Fjord",
    
    "Germany, Wismar Bight",
    "Jutland, Limfjord",
    "Jutland, Mariager Fjord",
    "Mediterranean",
    "NEAtlantic",
    "Netherlands",
    "AtlanticOcean:NWAtlantic",
    "Samsøe, Ballen",
    "Sealand, Skovshoved",
    "USA:GalvestonBay",
    "USA:Panacea",
    "USA:WoodsHole"
  )
  lngNmtk <- paste(lngNmtk, collapse = "|")
  df_4LtabbrNm.01 <-  df_4LtabbrNm[grepl(lngNmtk, df_4LtabbrNm$llongNm), ]
  # substitute in the list of  locations
  df_4LtabbrNm.01$llongNm <- gsub("BalticSea", "Baltic Sea", df_4LtabbrNm.01$llongNm)
  df_4LtabbrNm.01$llongNm <- gsub("CaspianSea", "Caspian Sea", df_4LtabbrNm.01$llongNm)
  df_4LtabbrNm.01$llongNm <- gsub("CentralWAtlantic",
                                  "Central W Atlantic",
                                  df_4LtabbrNm.01$llongNm)
  df_4LtabbrNm.01$llongNm <- gsub("NEAtlantic", "NE Atlantic", df_4LtabbrNm.01$llongNm)
  df_4LtabbrNm.01$llongNm <- gsub("AtlanticOcean:NWAtlantic",
                                  "NW Atlantic",
                                  df_4LtabbrNm.01$llongNm)
  df_4LtabbrNm.01$llongNm <- gsub("USA:", "NW Atlantic, USA:", df_4LtabbrNm.01$llongNm)
  df_4LtabbrNm.01$llongNm <- gsub(":", ", ", df_4LtabbrNm.01$llongNm)
  df_4LtabbrNm.01$llongNm <- gsub(", ", " ", df_4LtabbrNm.01$llongNm)
  #paste long name and abbreviation together to use for table header
  df_4LtabbrNm.01$llongNm_and_abbr <- paste0(df_4LtabbrNm.01$llongNm, " (", df_4LtabbrNm.01$llabbr, ")")
  # paste all names and abbreviations into one string
  ablo9 <- paste(df_4LtabbrNm.01$llongNm_and_abbr, collapse = ", ")
  # make a table caption
  capt_tbl02 <-        paste0(
    "Table 6. Estimates for diversity and neutrality for the populations in the analysed ",
    GnNm,
    " regions of samples for Mnemiopsis leyidi. Abbreviations above columns are: Number of haplotypes per population (Nh), haplotypic diversity (h), and standard deviation for haplotypic diversity (h.sd), nucleotide diversity (nd), and standard deviation for nucleotide diversity (nd.sd), probability of Tajimas D (p.TD). Sampled locations are abbreviated: ",
    ablo9,
    ". Additional sequences were obtained from M. leidyi from NBCI GenBank as indicated by accession numbers."
  )
  #make a function that keeps 2 decimal places : see: https://stackoverflow.com/questions/48341878/increasing-decimal-positions-swirl-r-programming-environment-12-data-manip
  fSpr2 <- function(c)
    sprintf("%.2f", c)
  fSpr0 <- function(c)
    sprintf("%.0f", c)
  # apply the function to the data frame -  this makes it a matrix : see: https://www.tutorialkart.com/r-tutorial/r-apply-function-to-each-element-of-matrix/
  mtx_popinf02 <-  apply(df_popinf02, 2, fSpr2)
  mtx_popinf02[, 1:2] <- as.character(round(as.numeric(as.character(mtx_popinf02[, 1:2])), digits =
                                              2))
  
  # match to get coordinates for sampling positions
  smplcoord <- df_clo04$dec_latlon2[match(poplocNms, df_clo04$locality8)]
  smplcoord[is.na(smplcoord)] <- ""
  # write out the alignment as a fasta file
  outF_tmp_fasta <- paste0(wd00_wd05, "/tmp_algnM.fasta")
  alvrs.M <- lst_alvrs.M.ITS[[ng]]
  bioseq::write_fasta(
    alvrs.M,
    file = outF_tmp_fasta,
    line_length = Inf,
    block_length = Inf
  )
  # read in the fasta file again
  pth_inpf03 <- outF_tmp_fasta
  #read FASTA as multi.dna
  mltdn_pip03 <- apex::read.multiFASTA(pth_inpf03)
  # rowbind the nested lists to a dataframe
  df_pip03 <- data.frame(do.call
                         ('rbind', strsplit(as.character(
                           mltdn_pip03@labels
                         ), "_")))
  
  
  # assign  columns copies to havevolumns with more descriptive names
  df_pip03$locality <- df_pip03$X2
  df_pip03$smplNm <- df_pip03$X1
  # order data frame by location name and then by individual sample name
  df_pip03.2 <- df_pip03[order(df_pip03$locality, df_pip03$smplNm), ]
  # get  only the Mnelei sample numbers
  Mnl_smpls  <- df_pip03.2$smplNm[grepl("Mnelei",df_pip03.2$smplNm)]
  # get the corresponding NCBI accession numbers
  Mnl_smpls <- df_no_ncbi$NCBIAccNo[match(Mnl_smpls,df_no_ncbi$MneleiNo)]
  # add back in to the data frame
  df_pip03.2$smplNm[grepl("Mnelei",df_pip03.2$smplNm)] <- Mnl_smpls
  # get lettercode in sample and get number code in samples
  df_pip03.2$ltc.smpl <- gsub("^([A-Za-z]{+})([0-9]{+})$", "\\1", df_pip03.2$smplNm)
  df_pip03.2$noc.smpl <- gsub("^([A-Za-z]{+})([0-9]{+})$", "\\2", df_pip03.2$smplNm)
  # make the number numeric, to be albe to find the maximum and the minimum
  df_pip03.2$noc.smpl <- as.numeric(df_pip03.2$noc.smpl)
  
  
  
  # use dplyr to summarise upper and lower value of sample per group
  # to use later on in a table that holds all samples
  df_pip03.3 <- df_pip03.2 %>%
    # determine the columns to work on with dplyr
    dplyr::select(locality, ltc.smpl, noc.smpl) %>%
    # specify the columns to group by
    dplyr::group_by(locality, ltc.smpl) %>%
    # get the max value of the noc.smpl column per groups defined above, and
    # put in a new column called 'mxnb', and also do this for the minimum value
    dplyr::summarise(mxnb = max(noc.smpl), minb = min(noc.smpl))
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  df_pip03.3$mxnb <- ifelse(
    nchar(df_pip03.3$mxnb) < 3,
    stringr::str_pad(df_pip03.3$mxnb, 3, pad = "0"),
    df_pip03.3$mxnb
  )
  df_pip03.3$minb <- ifelse(
    nchar(df_pip03.3$minb) < 3,
    stringr::str_pad(df_pip03.3$minb, 3, pad = "0"),
    df_pip03.3$minb
  )
  # paste together to get range of sample numbers per sample location
  df_pip03.3$smplrng <- paste0(df_pip03.3$ltc.smpl,
                               df_pip03.3$minb,
                               "-",
                               df_pip03.3$ltc.smpl,
                               df_pip03.3$mxnb)
  # concatenate strings per group - see this example: https://datacornering.com/how-to-concatenate-text-by-group-in-r/
  # notice I added'dplyr::' to the dplyr functions, as I also have 'plyr' up and running
  require(dplyr)
  df_pip03.3 <- df_pip03.2 %>%
    dplyr::select(smplNm, locality) %>%
    dplyr::group_by(locality) %>%
    dplyr::mutate(all_spmls = paste(smplNm, collapse = ", "))  %>%
    dplyr::distinct(locality, all_spmls)
  # copy the column to a vector
  all_spmlsNos <- c(df_pip03.3$all_spmls, "")
  
  
  #get the local sampled label names
  ownsmplN2 <- names(alvrs.M)[grepl("Mnelei", names(alvrs.M))]
  lngNm <- ownsmplN2
  lngNmspl <- strsplit(as.character(lngNm), "_")
  smplyear4 <- as.factor(sapply(lngNmspl, "[[", 3))
  smplmnt4 <- as.character(as.factor(sapply(lngNmspl, "[[", 4)))
  # use substitute to get hte different elements
  MnsmNo <- gsub("^(Mnelei[0-9]{+})_.*", "\\1", ownsmplN2)
  # make data frame with dates and months and years for own samples
  df_ownsmpldt <- as.data.frame(cbind(MnsmNo, smplyear4, smplmnt4, lngNm))
  
  # get sample month for samples
  ncbismplmnth <- df_lN02$smplyear[match(df_pip03$smplNm[grepl("[A-Z]{2}[0-9]{5}", df_pip03$smplNm)], df_lN02$accession_nmb)]
  smplmnth.ncbi <- cbind(df_pip03$smplNm[grepl("[A-Z]{2}[0-9]{5}", df_pip03$smplNm)], ncbismplmnth)
  df_smplmnth.ncbi <- as.data.frame(smplmnth.ncbi)
  colnames(df_smplmnth.ncbi) <- c("smplNo", "mnth")
  ownsmplmnth <- df_ownsmpldt$smplmnt4[match(df_pip03$smp[grepl("Mnelei", df_pip03$smplNm)], df_ownsmpldt$MnsmNo)]
  smplmnth.Mnelei <- cbind(df_pip03$smp[grepl("Mnelei", df_pip03$smplNm)], ownsmplmnth)
  df_smplmnth.Mnelei <- as.data.frame(smplmnth.Mnelei)
  colnames(df_smplmnth.Mnelei) <- c("smplNo", "mnth")
  #put all data frames into list
  df_list <- list(df_smplmnth.ncbi, df_smplmnth.Mnelei)
  #merge all data frames together
  df_smpl.mnth <- Reduce(function(x, y)
    merge(x, y, all = TRUE), df_list)
  # match back to pip dataframe
  df_pip03$month <- df_smpl.mnth$mnth[match(df_pip03$smplNm, df_smpl.mnth$smplNo)]
  # make an empty column to fill
  df_pip03$month2 <- NA
  df_pip03$year <- df_pip03$X3
  # evaluate if year and month column are the same,
  # and then use the month or nothing
  df_pip03$month2[!df_pip03$month == df_pip03$year] <-
    df_pip03$month[!df_pip03$month == df_pip03$year]
  df_pip03$month2[df_pip03$month == df_pip03$year] <- ""
  # assign to vectors
  
  syear <- df_pip03$year[match(poplocNms, df_pip03$locality)]
  syear[is.na(syear)] <- ""
  smnth <- df_pip03$month2[match(poplocNms, df_pip03$locality)]
  smnth[is.na(smnth)] <- ""
  # get coordinates
  smplcoord <-  df_clo04$latlonpos3[match(poplocNms, df_clo04$loc4Lett)]
  
  # bind the population location names back on to the matrix array
  mtx_popinf02 <- cbind(poplocNms,
                        smplcoord,
                        mtx_popinf02,
                        syear,
                        smnth,
                        all_spmlsNos)
  # add back row names
  rownames(mtx_popinf02) <- poplocNms
  #View(mtx_popinf02)
  ncol(mtx_popinf02)
  # make a vector with the new column names
  ch_colNms <- c(
    "Population",
    "Coordinates",
    "Sample size",
    "Nh",
    "h",
    "h.sd",
    "nd",
    "nd.sd",
    "Tajima's D",
    "p.TD",
    "Fu's F",
    "year",
    "month",
    "sample No"
  )
  #
  length(ch_colNms)
  # change the column names
  colnames(mtx_popinf02) <- ch_colNms
  # count the number of columns in the matrix
  nclmtx <- ncol(mtx_popinf02)
  #make a vector for text adjustment
  vftxtadj <- c("l", rep("r", nclmtx - 1))
  
  # show the table
  t.HTML05 <- mtx_popinf02 %>%
    addHtmlTableStyle(align = vftxtadj) %>%
    #addHtmlTableStyle(css.cell = colourhtml) %>%
    #addHtmlTableStyle(css.cell = colourhtml_cell) %>%
    #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
    htmlTable(caption = capt_tbl02, rnames = FALSE)
  t.HTML05
  #pad with zeros to three characters for own Mnelei samples
  #see this website: https://stackoverflow.com/questions/5812493/adding-leading-zeros-using-r
  pnng <- ifelse(nchar(ng) < 2, stringr::str_pad(ng, 2, pad = "0"), ng)
  # hcange the row names so that it includes the gene analysed
  row.names(mtx_popinf02) <- paste0(GnNm,", ",row.names(mtx_popinf02))
  # collect the table in a list
  lst.nucdiv.tbls[[ng]] <- mtx_popinf02
  
  # ensure the package 'htmltools' is loaded to be able to
  # save the html table
  library(htmltools)
  # Make a filename to store the html table with Fst vaules
  filNm.for_html <- paste0(wd00_wd05,
                           "/Table06_v",
                           pnng,
                           "_html_table_haplotypes_for_",
                           GnNm,
                           ".html")
  
  htmltools::save_html(t.HTML05, file = filNm.for_html)
  
  
  
  # collect the Hptfreq table in the list of tables
  lst.nucdiv.tbls[[ng]] <-   mtx_popinf02
  
  
  # use the DNAbin object to get the frequency of haplotypes per population
  # this fucntion must have the rownames in the DNAbin object, otherwise,
  # it cannot work.
  mtx_HptFrq <- pegas::haploFreq(dnb_pip05,
                                 split = "_",
                                 what = 2,
                                 haplo = hpt.pip)
  
  
  hgstHptNmb <- max(as.numeric(as.roman(as.character(
    row.names(mtx_HptFrq)
  ))))
  row.names(mtx_HptFrq) <- as.numeric(as.roman(as.character(row.names(mtx_HptFrq))))
  
  df_HptFrq2 <- as.data.frame(mtx_HptFrq)
  totcnt.Hpt <- colSums(df_HptFrq2)
  df_totcnt.Hpt <- as.data.frame(t(totcnt.Hpt))
  # also add the total count of haplotypes
  df_HptFrq3 <- rbind(df_HptFrq2, df_totcnt.Hpt)
  
  df_HptFrq3$HptNmb <- rownames(df_HptFrq3)
  
  df_HptFrq3$HptNmb[as.numeric(df_HptFrq3$HptNmb) > hgstHptNmb] <- "Total"
  lastcol <- ncol(df_HptFrq3)
  df_HptFrq3 <- df_HptFrq3[, c(lastcol, seq(1, lastcol - 1))]
  colnames(df_HptFrq3)[1] <- c("Genotype")
  df_HptFrq3[df_HptFrq3 == 0] <- ""
  
  # make a table caption
  capt_tbl02 <-        paste0(
    "Table 7. Genotype frequency per population. For ",
    GnNm,
    ". Abbreviations for sampling locations are explained in table 1."
  )
  
  # show the table
  t.HTML06 <- df_HptFrq3 %>%
    addHtmlTableStyle(align = "r") %>%
    #addHtmlTableStyle(css.cell = colourhtml) %>%
    #addHtmlTableStyle(css.cell = colourhtml_cell) %>%
    #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
    htmlTable(caption = capt_tbl02, rnames = FALSE)
  t.HTML06
  
  # Make a filename to store the html table with Fst vaules
  filNm.for_html <- paste0(
    wd00_wd05,
    "/Table07_v",
    pnng,
    "_html_table_haplotypes_frequency_for_",
    GnNm,
    ".html"
  )
  
  htmltools::save_html(t.HTML06, file = filNm.for_html)
  # rename the entries in the genotype column
  df_HptFrq3$Genotype <- paste0(GnNm,", ",df_HptFrq3$Genotype)
  # collect the Hptfreq table in the list of tables
  lst.frq.tbls[[ng]] <-   df_HptFrq3
  

  # end iteration over alignments in list
}

#make the list a data frame
tibl_frq <- bind_rows(lst.frq.tbls)
# replace NAs in the tibble
tibl_frq <- tibl_frq %>% replace(is.na(.), "")
# make a table caption
capt_tbl02 <-        paste0(
  "Table 8. Genotype frequency per population. For ",
  "ITS1 and ITS2",
  ". Abbreviations for sampling locations are explained in table 1."
)
# show the table
t.HTML06 <- tibl_frq %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  #addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02, rnames = FALSE)
t.HTML06

# Make a filename to store the html table with Fst vaules
filNm.for_html <- paste0(
  wd00_wd05,
  "/Table08_v01",
  "_html_table_haplotypes_frequency_for_",
  "ITS1_and_ITS2",
  ".html"
)

htmltools::save_html(t.HTML06, file = filNm.for_html)

# Combine the two nucleotide diversity tables
df_nucdiv <- as.data.frame(do.call(rbind, lst.nucdiv.tbls))
# use the row names for the Populations 
df_nucdiv$Population <-  row.names(df_nucdiv)
# identify the columns not to keep
clidxN <- which(grepl(c("Coordinates|month|year|sample No"),colnames(df_nucdiv)))
df_locations <- df_nucdiv[,c(clidxN)]

# subset the data frame
df_nucdiv <- df_nucdiv[,-c(clidxN)]


# change the column that holds the gene region and the location sampled 
colnames(df_nucdiv)[(grepl("Population",colnames(df_nucdiv)))] <- "Gene, sampled location"
# make a table caption
capt_tbl02 <-        paste0(
  "Table 9. Estimates for diversity and neutrality for the populations in the analysed ",
  "ITS1 and ITS 2",
  " regions of samples for Mnemiopsis leyidi. Abbreviations above columns are: Number of haplotypes per population (Nh), haplotypic diversity (h), and standard deviation for haplotypic diversity (h.sd), nucleotide diversity (nd), and standard deviation for nucleotide diversity (nd.sd), probability of Tajimas D (p.TD). Sampled locations are abbreviated: ",
  ablo9,
  "."
)

# show the table
t.HTML06 <- df_nucdiv %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  #addHtmlTableStyle(css.cell = colourhtml_cell) %>%
  #  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  htmlTable(caption = capt_tbl02, rnames = FALSE)
t.HTML06

# Make a filename to store the html table with Fst vaules
filNm.for_html <- paste0(
  wd00_wd05,
  "/Table09_v01",
  "_html_nucleotide_div_",
  "ITS1_and_ITS2",
  ".html"
)

htmltools::save_html(t.HTML06, file = filNm.for_html)

# substitute in the rownames, and use this to make 
# a new column with regional location names
df_locations$locNm <- gsub("ITS., (.*$)","\\1",row.names(df_locations))
# remove duplicated rows
df_locations <- df_locations[!duplicated(df_locations$locNm),]

df_locations <- df_locations[!(grepl("Total",df_locations$locNm)),]
# identify the column with "locNm"
idxnloc_clm <- which(grepl("locNm",colnames(df_locations)))
ncl_df_loc <- ncol(df_locations)
sq_ncl_dfloc <- seq(1:ncl_df_loc)
cl_df_loc_tm <- sq_ncl_dfloc[!( sq_ncl_dfloc %in% idxnloc_clm)]
nclord_dfloc <- c(idxnloc_clm,cl_df_loc_tm)
# re order the columns
df_locations <- df_locations[,nclord_dfloc]
# replace the column name in the column with "locNm"
colnames(df_locations)[grepl("locNm",colnames(df_locations))] <- "sampled location"
# identify the column with sampling month and exlude the column
idxnloc_clm <- which(grepl("month",colnames(df_locations)))
df_locations <- df_locations[,-c(idxnloc_clm)]
# make a table caption
capt_tbl02 <-        paste0(
  "Table 10. Samples for Mnemiopsis leyidi. Sampled locations are abbreviated: ",
  ablo9,
  ". ‘Sample No’ refers to the NCBI Accession numbers. Sequences obtained in this study begins with ‘PP’."
)

# show the table
t.HTML06 <- df_locations %>%
  addHtmlTableStyle(align = "r") %>%
  htmlTable(caption = capt_tbl02, rnames = FALSE)
t.HTML06

# Make a filename to store the html table with Fst vaules
filNm.for_html <- paste0(
  wd00_wd05,
  "/Table10_v01",
  "_sample_locations_",
  "ITS1_and_ITS2",
  ".html"
)

htmltools::save_html(t.HTML06, file = filNm.for_html)

#
#
#
