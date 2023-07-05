###################### Differential analysis of peaks 
### previous processing steps (steps 5-9)
## 1. concatenate individual peak files, trim sequences around centre, sort and merge, trim again
## 1b. compare segments of equal length to avoid huge size bias between comparisons.
## 2. calculate matrix
### in this R notebook
## 3. calculate TOTAL counts over regions, per sample
## 4. run EdgeR on count tables


library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(VennDiagram)
library(RColorBrewer)


########## 1. Import data and format tables for analysis. ####
#select sample to analyse
sample.name<-"H3K4me1_enhancers"

samples<-c("norm_CT46-2i_H3K4me1","norm_CT47-DMSO_H3K4me1","norm_CT48-PD03-siNanog_H3K4me1",
           "norm_CT49-2i_H3K4me1","norm_CT50-DMSO_H3K4me1","norm_CT51-PD03-siNanog_H3K4me1",
           "norm_CT52-2i_H3K4me1","norm_CT53-DMSO_H3K4me1","norm_CT54-PD03-siNanog_H3K4me1")


file.name<-paste(sample.name, "__Matrix.txt", sep="")

setwd("...") 

H3K37_all_matrix <- read.delim(file.name, header=FALSE, comment.char="@")

#reformat table and order for analysis
samples_clean<-samples
samples_clean[grepl("DMSO", samples_clean, fixed = TRUE)] <- "DMSO"
samples_clean[grepl("PD03", samples_clean, fixed = TRUE)] <- "PDsiNanog"
samples_clean[grepl("2i", samples_clean, fixed = TRUE)] <- "t2i"

num<-c(1,1,1,2,2,2,3,3,3)

sampleID<-paste(samples_clean, num, sep="")

mat_small<-H3K37_all_matrix[,7:ncol(H3K37_all_matrix)]


#sums all the reads over a certain region
t1<-apply(mat_small[1:nrow(mat_small),1:100], 1, sum)
t2<-apply(mat_small[1:nrow(mat_small),101:200], 1, sum)
t3<-apply(mat_small[1:nrow(mat_small),201:300], 1, sum)
t4<-apply(mat_small[1:nrow(mat_small),301:400], 1, sum)
t5<-apply(mat_small[1:nrow(mat_small),401:500], 1, sum)
t6<-apply(mat_small[1:nrow(mat_small),501:600], 1, sum)
t7<-apply(mat_small[1:nrow(mat_small),601:700], 1, sum)
t8<-apply(mat_small[1:nrow(mat_small),701:800], 1, sum)
t9<-apply(mat_small[1:nrow(mat_small),801:900], 1, sum)

regionID<-rep(1:nrow(mat_small),1)

means<-cbind(regionID, t1, t2, t3, t4, t5, t6, t7, t8, t9)
colnames(means)<-c("regionID", sampleID)

mm<-as.data.frame(means)
mm <- mm %>% 
  select(sort(names(.)))

meta_data<-data.frame(H3K37_all_matrix[,1:4], mm)

rm(t1, t2, t3, t4, t5, t6, t7, t8, t9, mat_small, regionID)

sample_data<-data.frame("Sample"=sampleID, "Genotype"=samples_clean, "Batch"=num)

#write.table(meta_data, "H3K27me3_metaGenes.txt", sep="\t", quote=F)
#write.table(sample_data, "H3K27me3_metaSamples.txt", sep="\t", quote=F)

####### 2. use EdgeR to get list of differentially bound regions ####

library(edgeR)
library(limma)
library(statmod)

dc <- DGEList(counts=means[,-1], group=sample_data$Genotype)

#generate contrast matrix for comparison
group <- sample_data$Genotype
design <- model.matrix(~0+group, data=sample_data$Sample)
colnames(design) <- levels(group)
rownames(design) <- make.names(sample_data$Sample)
design

contrasts <- makeContrasts(
  cont1 = t2i - DMSO,
  cont2 = t2i - PDsiNanog,
  cont3 = DMSO - PDsiNanog,
  levels = design
)

dc <- estimateDisp(dc, design, robust = TRUE)
#plotBCV(dc)

fit <- glmQLFit(dc, design, robust = TRUE)
#plotQLDisp(fit)

setA <- glmQLFTest(fit, contrast = contrasts[, "cont1"])
setB <- glmQLFTest(fit, contrast = contrasts[, "cont2"])
setC <-
  glmQLFTest(fit, contrast = contrasts[, "cont3"])


##### annotate each DE table with metadata

anno <- function(x) {
  xx <- as.data.frame(topTags(x, n = nrow(x), sort.by = "none"))
  regionID <- rownames(xx)
  extr <- cbind(regionID, xx)
  extr_anno <-
    merge(extr,
          meta_data,
          by = "regionID",
          all.x = TRUE,
          all.y = FALSE)
  rownames <-
    return(extr_anno)
}

setA_anno <- anno(setA)
setA_anno <- setA_anno[order(setA_anno$FDR),]

setB_anno <- anno(setB)
setB_anno <- setB_anno[order(setB_anno$FDR),]

setC_anno <- anno(setC)
setC_anno <- setC_anno[order(setC_anno$FDR),]

#generate smaller tables with each set of genes that go up and down in the pairwide comparisons
setA_up <- setA_anno[(setA_anno$logFC > 1), ]
setA_up <- setA_up[(setA_up$FDR < 0.01), ]
setA_up <- setA_up[order(setA_up$FDR), ]

setA_down <- setA_anno[(setA_anno$logFC < -1), ]
setA_down <- setA_down[(setA_down$FDR < 0.01), ]
setA_down <- setA_down[order(setA_down$FDR), ]

nam3 <- paste(sample.name, "_2i_PDsiNanog_DE.txt", sep = "")
write.table(setB_anno, nam3, sep = "\t", quote = F)

setB_up <- setB_anno[(setB_anno$logFC > 1), ]
setB_up <- setB_up[(setB_up$FDR < 0.01), ]
setB_up <- setB_up[order(setB_up$FDR), ]


setB_down <- setB_anno[(setB_anno$logFC < -1), ]
setB_down <- setB_down[(setB_down$FDR < 0.01), ]
setB_down <- setB_down[order(setB_down$FDR), ]

nam5 <- paste(sample.name, "DMSO_PDsiNanog_DE.txt", sep = "")
write.table(setC_anno, nam5, sep = "\t", quote = F)


setC_up <- setC_anno[(setC_anno$logFC > 1), ]
setC_up <- setC_up[(setC_up$FDR < 0.01), ]
setC_up <- setC_up[order(setC_up$FDR), ]


setC_down <- setC_anno[(setC_anno$logFC < -1), ]
setC_down <- setC_down[(setC_down$FDR < 0.01), ]
setC_down <- setC_down[order(setC_down$FDR), ]


### generate table for site that remain bound during differentiation 
# - i.e. sites high in both DMSO and 2i (>500 - ~10-11peak height), and not DE (FDR>0.1, abs(logFC)<1, PValue>0.05)
setA_nonDE <-setA_anno[(setA_anno$FDR > 0.1), ]
setA_nonDE <-setA_nonDE[(abs(setA_nonDE$logFC) <1), ]
setA_nonDE <-setA_nonDE[(abs(setA_nonDE$PValue) >0.05), ]
xx1<- apply(setA_nonDE[,11:13], 1, mean)
setA_nonDE <-setA_nonDE[(xx1>500) == TRUE,]
xx2<-apply(setA_nonDE[,17:19], 2, mean)
setA_nonDE <-setA_nonDE[(xx2>500)==TRUE,]


###### 5. plot venn diagrams between various fractions to identify things that change ####
# a. genes that lose binding during differentiation

venn.diagram(
  x = list(
    setA_up$regionID, #up in 2i vs DMSO
    setB_up$regionID), #up in 2i vs PD03siNanog
  category.names = c(
    "2i-high_DMSO" ,
    "2i-high_PDsiNanog"
  ),
  filename = '#1-lose binding during differentiation.pdf',
  output = TRUE
)
      
# b. genes that gain binding in differentiation


venn.diagram(
  x = list(
    setA_down$regionID, #up in DMSO vs 2i
    setB_down$regionID #up in PD03siNanog vs 2i
  ),
  category.names = c(
    "2i_DMSO-high",
    "2i_PDsiNanog-high"
  ),
  filename = '#2-gain binding during differentiation.pdf',
  output = TRUE
)

# c. genes that are maintained in differentiation


venn.diagram(
  x = list(
    setA_nonDE$regionID, #regions that maintein high binding
    setB_up$regionID, #higher in 2i vs PD03siNanog
    setC_up$regionID #higher in DMSO vs PD03siNanog
  ),
  category.names = c(
    "mainteined in 2ivsDMSO",
    "2i-high_PDsiNanog" ,
    "DMSO-high_PDsiNanog"
  ),
  filename = '#3-maintein binding during differentiation.pdf',
  output = TRUE)

#c. uniquely/aberrant gain or further increase binding in PDsiNanog conditions

FC2<-function(a,b){
  aa<-apply(a, 1, mean)
  bb<-apply(b, 1, mean)
  logFC2<- log(aa/bb,2)
  return(logFC2)
}

xx1<- FC2(setC_anno[,17:19], setC_anno[,14:16])
setC_uniquePD <-setC_anno[(xx1<(-0.5)) == TRUE,]
setC_uniquePD <-setC_uniquePD[(setC_uniquePD$logFC) < (-0.5), ]
setC_uniquePD <-setC_uniquePD[(setC_uniquePD$FDR) < 0.05, ]


##### 4. Export individual list of regions that match rules ####
write.table(setC_uniquePD, "group10-aberrant_high_PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)

#### export list of regions that change or maintain binding during differentiation
write.table(setA_up, "group1-lose binding during diff.txt", sep="\t", quote=F, row.names = FALSE)
write.table(setB_up, "group2-gain binding during diff.txt", sep="\t", quote=F, row.names = FALSE)
write.table(setA_nonDE, "group3-maintein binding during diff.txt", sep="\t", quote=F, row.names = FALSE)

tt<- intersect(setA_up$regionID,setB_up$regionID)
g4<-setA_up[setA_up$regionID %in% tt,]
write.table(g4, "group4-lose binding during diff AND PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)

t5<-setdiff(setA_up$regionID,setB_up$regionID)
g5<-setA_up[setA_up$regionID %in% t5,]
write.table(g5, "group5-lose binding during diff but NOT in PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)

t6<-intersect(setA_down$regionID,setB_down$regionID)
g6<-setA_down[setA_down$regionID %in% t6,]
write.table(g6, "group6-gain binding during diff AND PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)

t7<-setdiff(setA_down$regionID,setB_down$regionID)
g7<-setA_down[setA_down$regionID %in% t7,]
write.table(g7, "group7-gain binding during diff but NOT in PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)

t8<-union(intersect(setA_nonDE$regionID, setB_up$regionID), intersect(setA_nonDE$regionID, setC_up$regionID))
g8<-setA_nonDE[setA_nonDE$regionID %in% t8,]
write.table(g8, "group8-maintein binding during diff BUT lose binding in PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)

t9<-setdiff(setdiff(setA_nonDE$regionID, setB_up$regionID), setC_up$regionID)
g9<-setA_nonDE[setA_nonDE$regionID %in% t9,]
write.table(g9, "group9-maintein binding during diff AND PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)

g11<-rbind(g5, setC_uniquePD)
write.table(g11, "group11-overall higher binding than expected in PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)

g12<-rbind(g7, g8)
write.table(g12, "group12-overall lower binding than expected in PDsiNanog.txt", sep="\t", quote=F, row.names = FALSE)


