#https://cran.r-project.org/web/packages/Hapi/vignettes/Hapi.html

### Install dependencies
if(!require(Hapi)){
  install.packages("Hapi")
  install.packages('devtools')
  install.packages('HMM')
  devtools::install_github('Jialab-UCR/Hapi')
}
library("Hapi")
library("HMM")
#library(DT)

### load example data
data(gmt)
rownames(gmt) <- gmt$pos
head(gmt)
# 3.2 Haplotype phasing step by step
# 3.2.1 Data preprocessing
# 3.2.1.1 Convert genotype coding style
### covert A/T/C/G to 0/1
hetDa <- gmt[,1:4]
ref <- hetDa$ref
alt <- hetDa$alt

gmtDa <- gmt[,-(1:4)]
gmtDa <- base2num(gmt = gmtDa, ref = ref, alt = alt)
head(gmtDa)

### define HMM probabilities
hmm = initHMM(States=c("S","D"), Symbols=c("s","d"), 
              transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
              emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
              startProbs = c(0.5,0.5))
hmm
### filter out genotyping errors
gmtDa <- hapiFilterError(gmt = gmtDa, hmm = hmm)
## 35 hetSNPs with potential genotyping errors are filtered out !
# 3.2.1.3 Select high-quality hetSNPs to form a framework

### select a subset of high-quality markers
gmtFrame <- hapiFrameSelection(gmt = gmtDa, n = 3) ###

## Number of hetSNPs in the framework: 21074
# 3.2.1.4 Imputation of missing data in the framework
### imputation
imputedFrame <- hapiImupte(gmt = gmtFrame, nSPT = 2, allowNA = 0)
## Number of hetSNPs after imputation: 20994

head(imputedFrame)

# 3.2.2 Draft haplotype inference
# 3.2.2.1 Majority voting for draft haplotype inference
### majority voting
draftHap <- hapiPhase(gmt = imputedFrame) ###
head(draftHap)

### check positions with cv-links
draftHap[draftHap$cvlink>=1,]

# 3.2.2.2 Proofreading of draft haplotypes
### identification of clusters of cv-links
cvCluster <- hapiCVCluster(draftHap = draftHap, cvlink = 2)
cvCluster

### determine hetSNPs in small regions involving multiple cv-links
filter <- c()
for (i in 1:nrow(cvCluster)) {
  filter <- c(filter, which (rownames(draftHap) >= cvCluster$left[i] & 
                               rownames(draftHap) <= cvCluster$right[i]))
}

length(filter)

### filter out hetSNPs in complex regions and infer new draft haplotypes
if (length(filter) > 0) {
  imputedFrame <- imputedFrame[-filter, ]
  draftHap <- hapiPhase(imputedFrame)
} 

finalDraft <- hapiBlockMPR(draftHap = draftHap, gmtFrame = gmtFrame, cvlink = 1)

head(finalDraft)
# 3.2.3 High-resolution haplotype assembly
consensusHap <- hapiAssemble(draftHap = finalDraft, gmt = gmtDa)

head(consensusHap)

consensusHap <- hapiAssembleEnd(gmt = gmtDa, draftHap = finalDraft, 
                                consensusHap = consensusHap, k = 300)
### Haplotype 1
hap1 <- sum(consensusHap$hap1==0)
### Haplotype 2
hap2 <- sum(consensusHap$hap1==1)
### Number of unphased hetSNPs
hap7 <- sum(consensusHap$hap1==7)

### Accuracy
max(hap1, hap2)/sum(hap1, hap2)

### find hetSNP overlaps
snp <- which(rownames(hetDa) %in% rownames(consensusHap))

ref <- hetDa$ref[snp]
alt <- hetDa$alt[snp]

### convert back to A/T/C/G
consensusHap <- num2base(hap = consensusHap, ref = ref, alt = alt)
head(consensusHap)

### output all the information
hapOutput <- data.frame(gmt[snp,], consensusHap)
head(hapOutput)

# 3.2.4 Visualization of haplotypes in single gamete cells
### load haplotypes in each gamete cell
data(gamete11)
head(gamete11)

### load chromosome information of the genome
data(hg19)
head(hg19)

### view gamete cells 
hapiGameteView(chr = hg19, hap = gamete11)

# 4 Crossover analysis module

# 4.1 Identification of crossovers


### haplotypes
hap <- hapOutput[,10:11]
head(hap)


### gametes
gmt <- hapOutput[,5:9]
head(gmt)

### identify crossover
cvOutput <- hapiIdentifyCV(hap = hap, gmt = gmt)
cvOutput

#4.2 Crossover visualization

### load crossover table
data(crossover)
head(crossover)

# 4.2.1 Visualization of crossover resolution
hapiCVResolution(cv = crossover)
# 4.2.2 Visualization of crossover distances
hapiCVDistance(cv = crossover)

#4.2.3 Visualization of crossover map

hapiCVMap(chr = hg19, cv = crossover, step = 5)











