## Tom Wolfe ##

library(pheatmap)
library(openxlsx)
tpm_table_path <- "" ## <--Provide the full path to the tpm count table 
                     ## (i.e. "C:/Users/[username]/Lab_7_2021_RNA_Seq/TPM_Across_8_Libraries_Ptep_3.xlsx")
tpm_8 <- read.xlsx(tpm_table_path)
tpm_8_summary <- summary(tpm_8)

all_counts <- as.vector(as.matrix(tpm_8[,-c(1,2)]))
min_nonzero_count <- min(all_counts[-which(all_counts == 0)])

## Identified nearby genes:
nearby_genes <- data.frame("Relative_Location"=c("Upstream",
                                        "Upstream",
                                        "Downstream",
                                        "Downstream",
                                        "Downstream",
                                        "Target_Gene"),
                           "Strand"=c("-","+","-","+","-","+"),
                           "Gene_Start_coord"=c(12905,3440,68712,92551,106486,59850),
                           "Gene_End_coord"=c(59672,14106,92557,105055,144602,68037),
                           "Transcript_Accession"=c("XM_016069198.3",
                                         "XM_016069200.2",
                                         "XM_016069202.2",
                                         "XM_016069203.2",
                                         "XM_016069194.2",
                                         "XM_016069201.3"),
                           "Gene_Symbol"=c("LOC107452653",
                                         "LOC107452654",
                                         "LOC107452656",
                                         "LOC107452657",
                                         "LOC107452650",
                                         "LOC107452655"),
                           "Description"=c("structural maintenance of chromosomes protein 4",
                                           "ribonuclease P protein subunit p29",
                                           "mediator of RNA polymerase II transcription subunit 17",
                                           "splicing factor YJU2",
                                           "very low-density lipoprotein receptor",
                                           "very-long-chain (3R)-3-hydroxyacyl-CoA dehydratase"))
counts <- as.data.frame(matrix(data = NA,
                 nrow = nrow(nearby_genes),
                 ncol = ncol(tpm_8)))
colnames(counts) <- colnames(tpm_8)
colnames(counts)[2] <- "Transcript_length"
for(i in 1:nrow(counts)){
  counts[i,] <- tpm_8[which(tpm_8$target_id == nearby_genes$Transcript_Accession[i]),]
}
counts <- counts[order(counts$target_id),]
nearby_genes <- nearby_genes[order(nearby_genes$Transcript_Accession),]
nearby_genes <- cbind(nearby_genes,counts)
nearby_genes <- nearby_genes[,-8]
medians <- t(as.data.frame(apply(tpm_8[,-c(1,2)], 2, median)))
row.names(medians) <- "median_count"



ratio_to_median <- as.data.frame(matrix(data = NA,
                          nrow = nrow(nearby_genes),
                          ncol = ncol(medians)))
colnames(ratio_to_median) <- colnames(medians)
row.names(ratio_to_median) <- nearby_genes$Gene_Symbol
for(i in 1:nrow(ratio_to_median)){
  ratio_to_median[i,] <- log(nearby_genes[i,-(1:8)]/medians)
}
for(i in 1:nrow(ratio_to_median)){
  ratio_to_median[i,which(is.infinite(as.matrix(ratio_to_median[i,])))] <- min(as.vector(as.matrix(ratio_to_median))[which(is.finite(as.vector(as.matrix(ratio_to_median))))])
}
pheatmap::pheatmap(t(ratio_to_median),
                   angle_col = 45,
                   main = sprintf("Transcript Count/Median Trancsript Count \n(logFC)"))
