if(!require(GHap)){
  install.packages("GHap")
}
library("GHap")
# https://rdrr.io/cran/GHap/f/vignettes/GHap_vignette.Rmd
library(vcfR)


lstfwd05 <-list.files(wd00_wd05, full.names = T)
grpsubfix <- "_v01.fasta"
all_vars.fs.fl  <- lstfwd05[grepl(grpsubfix,lstfwd05)]

fl <- all_vars.fs.fl[1]
dnb_pip <- ape::read.dna(fl, format = "fasta")

#make the DNAbin object a genind object
geni_pip <- adegenet::DNAbin2genind(dnb_pip)
#make the  genind object a loci object, as the pegas package
# can write a vcf format file from a 'loci' object
loci_pip <- pegas::genind2loci(geni_pip)
# write a vcf format file from a 'loci' object
pegas::write.vcf(loci_pip,file=paste0(wd00_wd05,"/tmp_vcf.txt"))
# read in a  vcf format file 
vcf_pip <- vcfR::read.vcfR(paste0(wd00_wd05,"/tmp_vcf.txt"))


# split the labels in the DNAbin object 
Nms <- strsplit(as.character(labels(dnb_pip)), "_")
# get the 2 nd element per vector in the list - this holds the abbreviated 
# location name and the sampling year and the gene variant 
poplocNm <- sapply(Nms, "[[", 2)

vcfR::genetic_diff(vcf_pip, pops = as.factor(poplocNm))

chrR_pip <- vcfR::create.chromR(vcf_pip)
vcfR::vcfR2tidy(vcf_pip)

GH_pip <- GHap::ghap.vcf2phase(vcf.files=c(paste0(wd00_wd05,"/tmp_vcf.txt")),
                               out.file =c(paste0(wd00_wd05,"/tmp_ghap.txt")) )

