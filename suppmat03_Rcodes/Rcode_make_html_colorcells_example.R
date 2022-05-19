
library(htmlTable)
library(dplyr)
# make ex ddf
le <- c(letters[1:4])
r1 <- round(rnorm(1:4,2),1)
r2 <- seq(1,8,2)
r3 <- round(rnorm(1:4,0.2),1)
r4 <- round(rnorm(1:4,0.01),1)
d1 <- as.data.frame(cbind(r1,r2,r3,r4))
rownames(d1) <- le
colnames(d1) <- le
# re arrange
d2 <- d1 %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
# make heat map
ggplot(d2, aes(x = rowname,
               y = colname,
               fill = value)) +
  geom_tile()


#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)
# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# add necessary html around each colour value
colourhtmlvalues <- paste0("background-color:",colourhtmlvalues ,";")
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
# show the table
mdlo %>%
  addHtmlTableStyle(align = "r") %>%
  addHtmlTableStyle(css.cell = colourhtml) %>%
  htmlTable

#_______________________________________________________________________________

#https://cran.r-project.org/web/packages/htmlTable/vignettes/general.html
mdlo <- as.matrix(d1)
nr <- nrow(mdlo)
nc <- ncol(mdlo)
# normalize values to a scale from 0 to 1
minval <- min(mdlo)
maxval <- max(mdlo)
valuesnorm <- mapply(function(x) (x-minval)/(maxval-minval), t(mdlo))
# set up colour scale using default ggplot palette
colscale<-scale_fill_continuous()
# alternative palettes
# note that we only need to define colscale once - comment out the ones you don't want
library(paletteer)
colscale <- scale_fill_paletteer_c(palette="scico::tokyo", direction = 1)
colscale <- scale_fill_paletteer_c(palette="pals::kovesi.linear_kryw_5_100_c64", direction = 1)

# get colour values corresponding to normalized values in matrix
colourhtmlvalues <- colscale$palette(valuesnorm)
# evaluate for low values and assign these a different color to use for the font
coltxtval <- ifelse(valuesnorm>0.3,"black","white")
# add necessary html around each colour value
colourhtmlvalues  <- paste0("background-color:",colourhtmlvalues ,";")
# also make  values for colors for text
coltxtval         <- paste0("color:",coltxtval ,";")


colourhtml_cell <- paste0(colourhtml,colourtxthtml)
colourhtml_cell <- matrix(data=colourhtml_cell,nrow=nr+1,ncol=nc+1)
#
# convert back to matrix
colourhtmlvalues <- t(matrix(data=colourhtmlvalues,nrow=nr,ncol=nc))
colourhtmltxtvals <- t(matrix(data=coltxtval,nrow=nr,ncol=nc))
#View(colourhtmlvalues)
# we now have colour codes for the cells containing values
# make a matrix the size of the table including headers and row names
ncells <- (nc+1)*(nr+1)
colourhtml <- rep("background-color:transparent;",ncells)
#colourtxthtml <- rep("color:transparent;",ncells)
colourtxthtml <- rep("color:#000000;",ncells)
# or use for example: "background-color:#FFFFFF;"
colourhtml <- t(matrix(data=colourhtml,nrow=nr+1,ncol=nc+1))
colourtxthtml <- t(matrix(data=colourtxthtml,nrow=nr+1,ncol=nc+1))

# replace the values part of the table with the html for the colours
colourhtml[2:(nr+1),2:(nc+1)]<-colourhtmlvalues
colourtxthtml[2:(nr+1),2:(nc+1)]<-colourhtmltxtvals

# show the table
mdlo %>%
  addHtmlTableStyle(align = "r") %>%
  #addHtmlTableStyle(css.cell = colourhtml) %>%
  addHtmlTableStyle(css.cell = colourhtml_cell) %>%
#  addHtmlTableStyle(css.cell = colourtxthtml) %>%
  
  htmlTable
#______________________________________________________________________________
#______________________________________________________________________________
#https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
if(!require(tableHTML)){
  install.packages("tableHTML")
  library(tableHTML)
}
require(tableHTML)
  tableHTML(mtcars)
  tbl.h.m <-   tableHTML(mtcars, rownames = FALSE)
  tableHTML(mtcars, class = 'table1')
  
  tableHTML::write_tableHTML(tbl.h.m, file = paste(wd00_wd05,"/Table02a_PhiST_locality.html",sep=""))
  
  nrow(mtcars)
  tableHTML(mtcars) %>%
    add_css_rows_in_column(css = list('background-color', 
                                      rep(c('red', 'green'), each = 16)),
                           column = 'mpg') %>%
    add_css_rows_in_column(css = list('background-color', 
                                      rep(c('green', 'red'), each = 16)),
                           column = 'cyl')
  
  
#try the tableHTML with no border
tH <- d1 %>% 
    tableHTML(border = 0) 
#count the number of columns in the dataframe
l.s.MO <- length(d1)
tH <- tH %>% 
      add_css_conditional_column(columns = 2:l.s.MO, #make it work on column 2 to the last column
                                 #conditional = "contains",
                                 #value = word,
                                 css = list(c("yellow","red")))
tH
  
  t.HTML17 <- tH