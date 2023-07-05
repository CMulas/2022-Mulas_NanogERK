#### 9. compute matrix to determine integrated signal over all possible regions


computeMatrix  reference-point \
--referencePoint center \
--upstream 750 \
--downstream 570 \
--binSize 15 \
--skipZeros \
-o H3K27me3promoter_matrix.gz \
--outFileSortedRegions regions1_H3K27me3promoter_matrix.bed \
-R /clean/5-H3K27me3promoter_peaks.txt \
-S /bigwig_files/norm_*_H3K27me3.bw  


