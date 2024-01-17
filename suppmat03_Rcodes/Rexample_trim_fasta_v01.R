
library(tidyverse)
library(ape)
#install.packages("ape")
if(!require(ape)){
  install.packages("ape")
  library(ape)
}

# get the working dir
prwd<- getwd()
# set working directory to be able to read in functions
setwd(paste0(prwd,"/suppmat03_Rcodes"))
source("Rfunction_DNA_sequence_subset.R")
source("Rfunction_ReadFasta.R")
setwd(prwd)

pth_inpf01 <- "algn_Mnelei_18s_16.fas.aligned.fasta"
pth_inpf02 <- "nt_algn_kyphosus_and_outgrp_mtdna_cytb_2023oct27.fasta"
fraction_ok <- 0.90
max_gap_width<- 100
min_grp_width <- 100
# ------------- Load data using ReadFasta function -------------
df <- ReadFasta(pth_inpf01)

# ------------- Calculate positions to be included (OK)  -------------
res <- DNA_sequence_subset(df,
                           required_fraction=fraction_ok,
                           max_gap_width=max_gap_width,
                           min_grp_width=min_grp_width)

df_sum <- res[[1]]
sListOK <- res[[2]]

# join df_sum to the input data marking which postions should be included
df <- df %>%
  left_join(select(df_sum,Posn,OK),by="Posn")

df_samples <- df %>%
  distinct(Sample)

nsamples <- nrow(df_samples)
nrequired <- fraction_ok*nsamples

# ------------- Output dataframe -------------
# output dataframe with sequences where OK==TRUE
df_out <- df %>% 
  filter(OK==TRUE) %>%
  select(-OK)
df_out <- df_out %>%
  pivot_wider(values_from=Value,names_from=Posn)


# ------------- Plot count of OK data at each position -------------
xlabel <- paste0("Position (OK: ",sListOK,")")

# plot the number of OKs at each position
p1 <- ggplot(data=df_sum,aes(x=Posn,y=Count)) +
  geom_point()+ 
  geom_hline(yintercept=nrequired,color="red",linetype=2) +
  ylab("Count")+
  xlab(label="") + 
  geom_ribbon(aes(ymin=0, ymax=(1-OK)*nsamples),alpha=0.1) +
  theme_minimal()

p1

df_plot <- df %>%
  mutate(Value=ifelse(Value %in% valuelist,Value,othervalues))

df_plot$Sample <- factor(df_plot$Sample,levels=rev(df_samples$Sample),labels=rev(substr(df_samples$Sample,1,ncrop)))
df_plot <- df_plot %>%
  mutate(a=ifelse(OK==T,1,0.7))

df_plot$Value <- factor(df_plot$Value,levels=c(valuelist,othervalues))
plotlabels <- c(toupper(paste0(valuelist,"     ")),othervalues)

p2 <- ggplot() +
  geom_raster(data=df_plot,aes(x=Posn,y=Sample,fill=Value)) + #
  ylab("Sample")+
  xlab(label="") + 
  scale_fill_manual(values=plotcolors,name="",labels=plotlabels) +
  scale_alpha_continuous(guide="none") +
  theme_minimal() +
  theme(axis.text.y=element_text(size=5),
        legend.text = element_text(size=8),
        legend.position="bottom",
        legend.key=element_rect(color="black"),
        legend.key.size=unit(0.6,"lines")) +
  guides(fill=guide_legend(nrow=1))
p2

#____________________________________________________________
# trim the alignment
#https://www.biostars.org/p/158250/
# read in the alignment
fnm_toread <- pth_inpf01
#fnm_toread <- inpf01
# read in the alignment
al_dt01 <- ape::read.dna(fnm_toread,format="fasta", as.matrix=TRUE)
# get the upper and lower limits to cut the alignment by
uppbp <- as.numeric(gsub(".*-","",sListOK))
lowbp <- as.numeric(gsub("-.*","",sListOK))
if (uppbp==Inf){
  uppbp <- ncol(al_dt01)
  lowbp <- 1}
#class(al_dt02)
#trim the alignment
al_dt02 <- al_dt01[,lowbp:uppbp]
#
row.names(al_dt02) <- df_pip03$nwseqNm[match(row.names(al_dt02),df_pip03$oriseqNm)]
