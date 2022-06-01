
# following the instructions on this webpage
#https://stackoverflow.com/questions/25755930/how-to-plot-pie-charts-in-haplonet-haplotype-networks-pegas
library(pegas)
library(ape)
data(woodmouse)
x <- woodmouse[sample(15, size = 110, replace = TRUE), ]
h <- haplotype(x)
net <- haploNet(h)
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.8)

ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))), 
  table(hap=ind, pop=rownames(x)[values])
)
ind.hap[1:10, 1:9]  #print just a chunk

plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.8, pie=ind.hap)
legend(50,50, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)

wrong.pop<-rep(letters[1:5], each=22)
ind.hap2<-with(
  stack(setNames(attr(h, "index"), rownames(h))), 
  table(hap=ind, pop=wrong.pop[values])
)

plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.8, pie=ind.hap2)
legend(50,50, colnames(ind.hap2), col=rainbow(ncol(ind.hap2)), pch=20)

#Modify the input data
pipx <- x
# get the labels
labppix <- labels(pipx)
#Make the dnabin a genind object
gei_pipx <- adegenet::DNAbin2genind(pipx)
#Make the genind object a data frame
df_geipipx <- adegenet::genind2df(gei_pipx)
# replace any NAs with dash
df_geipipx[is.na(df_geipipx)] <- "-"
# repeat row number 1 and 9 for 234 times
df_geipipx2 <- df_geipipx[rep((c(1,9)), each = 234), ]
# change a single nucleotide
df_geipipx2[3,8] <- "c"
# repeat the sequence in the data frame
df_geipipx3 <- df_geipipx2[rep((c(3)), each = 87), ]
# change a single nucleotide
df_geipipx2[4,9] <- "g"
# repeat the sequence in the data frame
df_geipipx4 <- df_geipipx2[rep((c(4)), each = 234), ]
# also repeat the labels
lbx2 <- rep(labppix[c(1,9)],234)
lbx3 <- rep(labppix[c(3)],87)
lbx4 <- rep(labppix[c(4)],233)
# sort the labels
lbx2 <- lbx2[order(lbx2)]
# add all labels to a single vector
labppix2 <- c(labppix,lbx2,lbx3,lbx4)
# bind the extra sequence in individual dataframes to the main data frame
df_geipipx2 <- rbind(df_geipipx,df_geipipx2,df_geipipx3,df_geipipx4)
#replace rownames to get new names
rownames(df_geipipx2) <- paste(labppix2,rownames(df_geipipx2),sep="_")
#make the data frame a matrix and make this a dnabin object
dnb_pipx02 <- ape::as.DNAbin(as.matrix(df_geipipx2))
#
pipx2 <- dnb_pipx02
h2 <- haplotype(pipx2)
#attributes(h2)
net2 <- haploNet(h2)
# repeat letters to have extra location names
wrpop_01 <-rep(letters[1:5], each=22)
wrpop_02<-rep(letters[1], each=(2*234))
wrpop_03<-rep(letters[3], each=(87))
wrpop_04<-rep(letters[3], each=(233))
wrpop3 <- c(wrpop_01,wrpop_02,wrpop_03,wrpop_04)

ind.hap3<-with(
  stack(setNames(attr(h2, "index"), rownames(h2))), 
  table(hap=ind, pop=wrpop3[values])
)

plot(net2, size=attr(net2, "freq"), scale.ratio = 0.021, cex = 0.8, pie=ind.hap3)
#legend(50,50, colnames(ind.hap3), col=rainbow(ncol(ind.hap3)), pch=20)

