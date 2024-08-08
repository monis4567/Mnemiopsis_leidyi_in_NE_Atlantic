
#remove everything in the working environment, without a warning!!
rm(list=ls())

library(devtools)
if(!require(bioseq)){
  # make sure you have Rtools installed
  if (!require('devtools')) install.packages('devtools')
  # install from GitHub
  remotes::install_github("fkeck/bioseq")
}
library(bioseq)

# #define input file as variable
inpf01 <- "Mnelei_seq_2024_aug_06_v01.fasta"
inpf02 <- "individual_seq_regions_for_ITS_and_SrRNA_in_Mnemiopsis_leydi_v01.fasta"


# read in the alignment with sequences
al01 <- bioseq::read_fasta(pth_inpf01)
# Read in a fasta file with gene regions
Fas_Ref_seq <- bioseq::read_fasta(pth_inpf02)

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
  return(list(st.pos.sta,st.pos.end))
  
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
# the first element in the list is the start position of the gene
# the second element is the end position of the gene
ITS2_s.p <- Get.st_en.of_geneR("ITS2",Fas_Ref_seq,12,al01)[[1]]
ITS2_e.p <- Get.st_en.of_geneR("ITS2",Fas_Ref_seq,12,al01)[[2]]

# use Bioseq 'seq_crop_pattern' function to cut out the block 
# with the coding ITS2 region
al01_ITS2 <- bioseq::seq_crop_position(al01,
                                       ITS2_s.p,
                                       ITS2_e.p)

