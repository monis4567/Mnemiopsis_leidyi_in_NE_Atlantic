# Function to read fasta format
# generates a data.frame with sample, position and value

ReadFasta<-function(filename,sample_identifier=">",sample_variable_name="Sample",position_variable_name="Posn",value_variable_name="Value"){

    # read the file
  data <- readLines(filename)
  
  # get the lines starting with the sample identifier character - the ae the sample names
  sample <- data[substr(data,1,1)==sample_identifier]
  sample <- substr(sample,2,nchar(sample))
  
  # get all the data lines (they still include the sample name)
  data  <- paste(data, collapse="")
  data <- strsplit(data,">")[[1]]
  data <- data[2:length(data)]
  
  # now remove the sample names from the data
  # for each sequence take the characters from n+1 to the end, where n is the length of the sample name 
  data <- substr(data,nchar(sample)+1,nchar(data))
  
  # create data frame with the sample name (Column 1) and sequence data (Column 2)
  df <- tibble(sample=sample,sequence=data)
  
  # get the length of the longest sequence
  n<-max(nchar(data))

  # separate sequence column into n separate columns, called C1 ... Cn
  df <- df %>% 
    separate(sequence, into = paste0("C",seq(1,n,1)), sep = seq(1,n-1,1))
  
  # pivot so that each letter is on a separate line
  # we now have 3 columns: sequence, position and value
  df <- df %>% 
    pivot_longer(c(1:n+1),names_to="position",values_to="value")
  
  # convert the position number from text "C[n]" to an integer n
  df <- df %>%
    mutate(position=as.integer(substr(position,2,nchar(position))))  
  
  # change the names of the columns to those specified in the function call
  names(df)[names(df)=="sample"]<-sample_variable_name
  names(df)[names(df)=="position"]<-position_variable_name
  names(df)[names(df)=="value"]<-value_variable_name
  
  # return the result
  return(df)
  
}