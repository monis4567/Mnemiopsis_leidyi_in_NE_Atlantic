
# following the instructions on this webpage
#https://stackoverflow.com/questions/25755930/how-to-plot-pie-charts-in-haplonet-haplotype-networks-pegas
library(pegas)
library(ape)

#_______________________________________________________________________________
#define input directories
wd01 <- "/suppmat01_inp_files"
wd00 <- "/home/hal9000/Documents/Documents/MS_Mnemiopsis/Mnemiopsis_leidyi_in_NE_Atlantic"
#paste directories together to get path for directory for input files
wd00_wd01 <- paste0(wd00,wd01)
#define inpus file
inpf01 <- "algn_Mnelei_18s_10.aligned.fasta.fas"
#paste path to directory together with filename
pth_inpf01 <- paste(wd00_wd01,"/",inpf01,sep="")
#read FASTA as dna.bin
x3 <- ape::read.dna(pth_inpf01, format = "fasta")
# make the haplotype object
h3 <- haplotype(x3)
net3 <- haploNet(h3)
#get location names
gl <- strsplit(as.character(rownames(x3)), "_")
#get first and second element of nested list
gnl <- sapply(gl, "[[", 2)
#replace row names
rownames(x3) <- gnl
# make the ind.hap object
ind.hap3<-with(
  stack(setNames(attr(h3, "index"), rownames(h3))), 
  table(hap=ind, pop=rownames(x3)[values])
)
# try plotting with different scale ratios
plot(net3, size=attr(net3, "freq"), scale.ratio = 2, cex = 0.8, pie=ind.hap3)
plot(net3, size=attr(net3, "freq"), scale.ratio = 20, cex = 0.8, pie=ind.hap3)
plot(net3, size=attr(net3, "freq"), scale.ratio = 0.20, cex = 0.8, pie=ind.hap3)
plot(net3, size=attr(net3, "freq"), scale.ratio = 500, cex = 0.8, pie=ind.hap3)
