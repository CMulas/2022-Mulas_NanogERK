library(dplyr)
library(data.table)

setwd("individual/")

marker<-"H3K27ac_peaks.txt" #select all individual BAM files, contain list of peaks for each file

files <- list.files(pattern = marker)

combined_files <- bind_rows(lapply(files, fread))

output_dir<-"combined/"

out_nam<-paste(output_dir, "combined_", marker, sep="")

write.table(combined_files, out_nam, sep = "\t", quote = F, row.names=FALSE, col.names = FALSE)