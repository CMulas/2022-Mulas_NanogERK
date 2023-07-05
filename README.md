# CUT&Tag
## Sample preparation.
Each batch consists of twelve samples processed in parallel. In each batch, the three samples [2i, N2B27/DMSO+siControl, N2B27/MEK(i)+siNanog] were processed in parallel for a given antibody, as well as 2 negative controls and one positive control (typically H3K27me3). Each sample consists of a single well of a 12-well plate. Cells were differentiated as described above and processed following the benchtop CUT&Tag protocol v2 (Kaya-Okur et al., 2019) dx.doi.org/10.17504/protocols.io.z6hf9b6) with minor modifications: all eppendorfs were pre-coated with PBS +0.1% BSA for 30min and all buffers were kept on ice prior to use. Nextera primers were used to generate libraries following the benchtop CUT&Tag recommended protocol. Antibodies and primers used are listed in Table S1. Libraries were pooled and sequenced 50bp pair-end Nextera.  

## Alignment and normalisation. 
Overview: Sample were processed in Galaxy (Afgan et al., 2022) workflows available in https://github.com/CMulas/CUT-Tag_analysis). FASTA files were aligned using Bowtie2 (Langmead and Salzberg, 2012) against mouse (mm10) and Escherichia coli K12 (spike in control). BAM files were then sorted, removing unaligned reads and PCR duplicates. For each antibody, all samples were normalised to 1x using bamCoverage (deepTools, (Ramírez et al., 2016).   

1) Align paired end-reads to mouse and E coli genome  
Code: 1_Galaxy-Workflow-CUTTag_alignment.ga  

2) Sort reads by coordinate, keep pair-end reads, remove PCR duplicates  
Code: 2_Galaxy-Workflow-CUTTag_filter_clean.ga  

3) Normalise coverage to 1x for each antibody across replicates  
Code: 3_Galaxy-Workflow-CUTTag_Normalise_per_antibody.ga  

## Differential binding analysis.  
Overview: For each sample, peaks were called using SEACR (Meers et al., 2019). To determine differential binding for a given antibody we first generated a file of all possible binding regions across samples. To generate this list, all individual lists of peaks were concatenated, sorted and merged. Next, we centred the peaks and trimmed all regions to 1000bp. Finally, we generate a count matrix of the total signal over the 1000 bp regions for each sample. Differential analysis was performed using EdgeR (Robinson et al., 2010). For H3K27me3, we divided regions into promoters and enhancers based on overlap with transcriptional start site (TSS). Peaks were annotated using PAVIS (Huang et al., 2013).   

4) Use SEACR to identify peaks in all samples  
Code: 4_Galaxy-Workflow-CUTTag_peaks.ga  

5) Generate a single combined list of peaks across all samples for each antibody  
Code: 5_R-concatenate-peak-BAM.R  

6) Sort and merge BAM files   
Code: 6_Python_sort_mergeBAM.txt  

7) Centre peaks reagions and trim to 1000bp.  
Code: 7_R-trim-centre-peaks.R  

8) Sort regions to obtain final list of regions (=regions where a peak was detected in at least one sample)  
Code: 8_Python_sort_mergeBAM copy.txt  

9) Compute matrix = measure signal across selected reagions in each sample  
Code: 9_Python_computeMatrix.py  

10) Perform differential binding analysis using EdgeR
Code: 10_R_DE-analysis.R   

# List of differentially bound regions  
group1 - lose binding during normal differentiation  
group2 - gain binding during normal differentiation  
group3 - maintain binding during normal differentiation  
group4 - lose binding during normal differentiation AND in MEK(i)+siNanog  
group5 - lose binding during normal differentiation but NOT in  MEK(i)+siNanog  
group6 - gain binding during normal differentiation AND in MEK(i)+siNanog  
group7 - gain binding during normal differentiation but NOT in  MEK(i)+siNanog  
group8 - maintain binding during normal differentiation BUT lose binding in MEK(i)+siNanog  
group9 - maintain binding during normal differentiation AND in MEK(i)+siNanog  
group10 - aberrant high binding MEK(i)+siNanog  
group11 - overall higher binding than expected in  MEK(i)+siNanog compared 2i or DMSO+siNeg  
group12 - overall lower binding than expected in  MEK(i)+siNanog compared 2i or DMSO+siNeg  
