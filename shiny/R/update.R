
updateFlashRPackages <- function(file='global.R')
{
   a <- readLines(file)
   x <- '(library)|(require)'
   b <- grep('(library)|(require)',a,value=TRUE)
   a <- gsub(x,"",b)
   b <- grep("#",gsub("[()]","",a),value=TRUE,invert=TRUE)
   if(length(b)>0) install.packages(b)
}
