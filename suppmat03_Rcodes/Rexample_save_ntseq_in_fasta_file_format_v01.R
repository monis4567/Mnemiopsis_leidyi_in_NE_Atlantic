#make names
nms <- paste(">seq_",LETTERS[1:4],sep="")
# make DNA seq
  vA <- rep("A", each=150)
  vC <- rep("C", each=150)
  vT <- rep("T", each=150)
  vG <- rep("G", each=150)
sqV <- c(vA,vC,vT,vG)
  sq <-  c( paste(sample(sqV,10),collapse = ""),
            paste(sample(sqV,5),collapse = ""),
            paste(sample(sqV,8),collapse = ""),
            paste(sample(sqV,9),collapse = "")
  )
# make it a dataframe
dfnt <-as.data.frame(cbind(nms,sq))

# dfnt er hvad jeg har til at begynde med
# jeg vil gerne ende med en data frame der kun har en eneste kolonne
# sådan her

# nmssq
# seq_A 
# CACAGGTGGC
# seq_B
# AGGTC
# seq_C
# TTTCGGCG
# seq_D
# GAATAAACT

library(tidyr)   # pivot_longer
library(dplyr)  # for at bruge pipes
dfnt_l <- dfnt %>% 
  pivot_longer(cols=c(nms,sq), values_to="rmssq", names_to = NULL)

# #write the sequences as a txt file , that is a fasta file
write.table(dfnt_l,  
            file = fnm_towrite, 
            quote = FALSE, sep ="\t" ,
            row.names = F, 
            col.names = F)

