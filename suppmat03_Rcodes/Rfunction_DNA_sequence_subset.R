DNA_sequence_subset<-function(df,PositionVariable=Posn,ValueVariable=Value,
                              MissingValue="-",required_fraction=0.5,
                              max_gap_width=100,min_grp_width=100,
                              multiple_blocks=F){
  
  df <- df %>%
    mutate(OK=ifelse({{ValueVariable}}==MissingValue,0,1))
  
  df_sum <- df %>%
    group_by({{PositionVariable}})%>%
    summarise(Count=sum(OK,na.rm=T))
  
  samplecount <- nrow(filter(df,{{PositionVariable}}==1))
  
  # number of samples required
  nrequired <- required_fraction*samplecount
  
  # mark postion as OK if the count of non-blank data is sufficient
  df_sum <- df_sum %>%
    mutate(OK=ifelse(Count>nrequired,T,F))
  
  # check if any positions are OK at all
  if(nrow(filter(df_sum,OK==T))==0){
    df_sum <- df_sum %>%
      select({{PositionVariable}},Count) %>%
      mutate(OK=F)
    sList <- paste0("NONE")
    return(list(df_sum,sList))
  }

  # add a variable indicating if the position just to the left is OK
  OKprev<- df_sum$OK[1:(nrow(df_sum)-1)]
  OKprev<-c(F,OKprev)
  df_sum$OKprev <- OKprev
  
  # find positions indicating the start of an "OK" block
  # previous position was not OK and this one is OK
  dfs1 <- df_sum %>%
    filter(OK==T & OKprev==F)
  dfs1 <- dfs1 %>%
    mutate(id=1:nrow(dfs1)) %>%
    select(id,Start={{PositionVariable}})
  
  # find positions indicating the end of an "OK" block
  # previous position was OK and this one is not OK
  dfs2 <- df_sum %>%
    filter(OK==F & OKprev==T)
  if(nrow(dfs2)>0){
    dfs2 <- dfs2 %>%
      mutate(id=1:nrow(dfs2)) %>%
      select(id,End={{PositionVariable}}) %>%
      mutate(End=End-1) # this is the first position which is not OK - we want the one just before
    if(nrow(dfs2)<nrow(dfs1)){
      dfs2 <- dfs2 %>%
        add_row(id=nrow(dfs1),End=nrow(df_sum))
    }
  }else{
    # if there are no OK=FALSE rows after the start of the last block of OK=TRUE
    # then we should stop at the 
    dfs2 <- data.frame(id=1,End=nrow(df_sum))
  }
  
  
  dfs <- dfs1 %>%
    left_join(dfs2,by="id")
  
  if(is.na(dfs[nrow(dfs),"End"])){
    dfs[nrow(dfs),"End"] <- nrow(df_sum)
  }
  
  
  if(nrow(dfs)==1){
    # only one block of data found
    pos1 <- dfs$Start[1]
    pos2 <- dfs$End[1]
    sList <- paste0(pos1,"-",pos2)
    df_sum <- df_sum %>%
      mutate(OK=ifelse({{PositionVariable}}>=pos1 & {{PositionVariable}}<=pos2,T,F))
    return(list(df_sum,sList))
  }
  
  # add a variable indicating end position of previous block
  EndPrev<- dfs$End[1:(nrow(dfs)-1)]
  EndPrev<- c(dfs$Start[1],EndPrev)
  dfs$EndPrev <- EndPrev
  
  dfs <- dfs %>% 
    mutate(Length=End+1-Start) %>%
    mutate(LenGap=Start-EndPrev)
  
  # we now have a table of distinct blocks of positions 
  # having more than the minimum number of non-blank data
  Grp<-1
  dfs$Grp<-NA
  for(i in 1:(nrow(dfs))){
    if(dfs$LenGap[i]>max_gap_width){
      Grp<-Grp+1
    }
    dfs$Grp[i]<-Grp
  }
  
  # we now assign each block to a group, starting at position 1
  # if a block is separated from the previous block 
  # by a number of blanks longer than max_gap_width,
  # then we start a new group
  dfgrp <- dfs %>%
    group_by(Grp) %>%
    summarise(Start=min(Start),End=max(End),Count=sum(Length)) %>%
    mutate(Length=End+1-Start) %>%
    ungroup()
  
  df_sum <- df_sum %>%
    select({{PositionVariable}},Count)
  
  # we now select only blocks longer than the specified minimum
  dflimits <- dfgrp %>%
    filter(Length>=min_grp_width)
  
  sList =""
  # if multiple blocks are allowed, each block meeting the minimum
  # length requirement is identified separately
  if(multiple_blocks==T){
    df_sum <- df_sum %>%
      mutate(OK=F)
    
    for(i in 1:nrow(dflimits)){
      pos1 <- dflimits$Start[i]
      pos2 <- dflimits$End[i]
      if(i==1){
        sList <- paste0(pos1,"-",pos2)
      }else{
        sList <- paste0(sList,", ",pos1,"-",pos2)
      }
      df_sum <- df_sum %>%
        mutate(OK=ifelse({{PositionVariable}}>=pos1 & {{PositionVariable}}<=pos2,T,OK))
    }
    
  }else{
    # if no multiple blocks are allowed
    # then one single blocks is defined covering all blocks
    # greater than the minimum width and including any blanks
    # between them
    dflimits <- dflimits %>%
      summarise(Start=min(Start),End=max(End)) %>%
      ungroup() 
    pos1 <- dflimits$Start[1]
    pos2 <- dflimits$End[1]
    sList <- paste0(pos1,"-",pos2)
    df_sum <- df_sum %>%
      mutate(OK=ifelse({{PositionVariable}}>=pos1 & {{PositionVariable}}<=pos2,T,F))
  }
  
  return(list(df_sum,sList))
}
