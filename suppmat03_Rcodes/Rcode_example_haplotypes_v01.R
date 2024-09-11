# https://search.r-project.org/CRAN/refmans/haplotypes/html/haplotypes-package.html


### Install dependencies
if(!require(haplotypes)){
  install.packages("haplotypes")
}
library("haplotypes")
## Read example FASTA file.	
f<-system.file("example.fas",package="haplotypes") 
# invalid character 'N' was replaced with '?' with a warning message
x<-read.fas(file=f)
# an object of class 'Dna'
x 

## or load DNA Sequence data set.
data("dna.obj") 
x<-dna.obj
## Not run: 
x

## End(Not run)

## Compute an absolute pairwise character difference matrix from DNA sequences.
# coding gaps using simple indel coding method
d<- distance(x,indels="sic") 
## Not run: 
d

## End(Not run)

## Infer haplotypes using the 'Dna' object.
# coding gaps using simple indel coding method
h<-haplotype(x,indels="s") 
## Not run: 
h

## End(Not run)

## Conduct statistical parsimony analysis with 95% connection limit.
#algortihmic method
## Not run: 
p<-parsimnet(x,prob=.95) 
p
# plot network
plot(p) 

## End(Not run)

## Plotting pie charts on the statistical parsimony network
## Not run: 
data("dna.obj")
x<-dna.obj
h<-haplotypes::haplotype(x)

## Statistical parsimony with 95
p<-parsimnet(x) 

#randomly generated populations
pop<-c("pop1","pop2","pop3","pop4","pop5","pop6","pop7","pop8")
set.seed(5)
pops<-sample(pop,nrow(x),replace=TRUE)

# Plotting with default parameters.
pieplot(p,h,1, pops)

#_______________________________________________________________________________
# Try out on own data
#_______________________________________________________________________________
# https://search.r-project.org/CRAN/refmans/haplotypes/html/haplotypes-package.html


### Install dependencies
if(!require(haplotypes)){
  install.packages("haplotypes")
}
library("haplotypes")
## Read example FASTA file.	
#f<-system.file("example.fas",package="haplotypes") 
lstfwd05 <-list.files(wd00_wd05, full.names = T)
grpsubfix <- "all_IUPACode_vars.fasta"
grpsubfix <- "_v01.fasta"
all_vars.fs.fl  <- lstfwd05[grepl(grpsubfix,lstfwd05)]

f <- all_vars.fs.fl[1]
# invalid character 'N' was replaced with '?' with a warning message
x<-read.fas(file=f)
# an object of class 'Dna'
x 

# ## or load DNA Sequence data set.
# data("dna.obj") 
# x<-dna.obj
# ## Not run: 
# x

## End(Not run)

## Compute an absolute pairwise character difference matrix from DNA sequences.
# coding gaps using simple indel coding method
d<- haplotypes::distance(x,indels="sic") 
## Not run: 
d

## End(Not run)

## Infer haplotypes using the 'Dna' object.
# coding gaps using simple indel coding method
h<-haplotypes::haplotype(x,indels="s") 
## Not run: 
h

## End(Not run)

## Conduct statistical parsimony analysis with 95% connection limit.
#algortihmic method
## Not run: 
p<-haplotypes::parsimnet(x,prob=.95) 
p
# plot network
plot(p) 

## End(Not run)

## Plotting pie charts on the statistical parsimony network
## Not run: 
# data("dna.obj")
# x<-dna.obj
h<-haplotypes::haplotype(x)
h
## Statistical parsimony with 95
p<-haplotypes::parsimnet(x) 
#randomly generated populations
set.seed(5)
Nms <- strsplit(as.character(labels(x)), "_")
# get the 2 nd element per vector in the list - this holds the abbreviated 
# location name and the sampling year and the gene variant 
poplocNm <- sapply(Nms, "[[", 2)
pops<-sample(poplocNm,nrow(x),replace=TRUE)

# Plotting with default parameters.
pieplot(p,h,1, pops)