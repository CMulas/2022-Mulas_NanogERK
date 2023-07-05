###### generate even peak intervals centered around 'centre' of peak
#### repeat for each file and save new trimmed file in a new folder

library(taRifx)
library(Rsamtools)

#input files = merged, sorted BAM files of all peaks, either at TSS or away
input_path <- "/combined/"

output_path <-"/clean/"
pattern_file<-"3-"
window.size<-1000



file_list <- list.files(path=input_path, pattern=pattern_file)
file_num <- length(file_list)
n <- 1

while (n <= file_num) {
  input_file<-paste(input_path, file_list[n], sep="")
  p2 <-read.delim(input_file, header = FALSE,  stringsAsFactors=FALSE)
  colnames(p2) <- c("chr", "start", "stop", "total", "peak", "maxSignalRegion")
  p2<-as.data.frame(p2)
  
  chpos<-regexpr(":", p2$maxSignalRegion)
  p2$peak.chr <- substr(p2$maxSignalRegion, 1, chpos-1)
  dashpos<-regexpr("-", p2$maxSignalRegion)
  p2$peak.start <- as.numeric(as.character(substr(p2$maxSignalRegion, chpos+1, dashpos-1)))
  p2$peak.end <- as.numeric(as.character(substr(p2$maxSignalRegion, dashpos+1, nchar(p2$maxSignalRegion))))
  
  p1 <- as.data.frame(p2[p2$peak > 10,])
  
  # # plot average features of peaks
  # hist(log(p1$stop - p1$start, 10))
  # median(p1$stop - p1$start)
  
  #find the peak size and the middle point of the peak
  p1$peak.dist <- p1$peak.end - p1$peak.start
  p1$peak.middle <- p1$peak.start + (round(p1$peak.dist / 2))
  
  #define new start and stop regions "windwo.ize" bp around peak
  p1$new.start <- p1$peak.middle - (window.size/2)
  p1$new.end <- p1$peak.middle + (window.size/2)
  p1<- p1[nchar(p1$peak.chr)<6,]
  
  #create new file with chromosome and new start and end coordinates
  p1.bam <- cbind(as.character(p1$chr), p1$new.start, p1$new.end)
  
  out_file <- paste(output_path, file_list[n], sep = "")
  write.table(p1.bam, out_file, sep = "\t", quote = F, row.names=FALSE, col.names = FALSE)

  n <- n+1
}
