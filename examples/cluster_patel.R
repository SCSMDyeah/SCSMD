### DATA
d <- read.table(".../data/Patel/GSE57872_GBM_data_matrix.txt")
# select 5 patients
d <- d[,grepl("MGH26_", colnames(d)) |
         grepl("MGH264_", colnames(d)) |
         grepl("MGH28_", colnames(d)) |
         grepl("MGH29_", colnames(d)) |
         grepl("MGH30_", colnames(d)) |
         grepl("MGH31_", colnames(d))]

### ANNOTATIONS
patients <- unlist(lapply(strsplit(colnames(d), "_"), "[[", 1))
patients[patients == "MGH264"] <- "MGH26"
ann <- data.frame(cell_type1 = patients)
rownames(ann) <- colnames(d)

### SINGLECELLEXPERIMENT
source(".../R/create_sce1.R")

sceset <- create_sce_from_logcounts(d, ann)
saveRDS(sceset, file = "patel.rds")


count_Patel <-readRDS(".../data/Patel/patel.rds")

counts <- exp(assay(count_Patel,"logcounts"))
assay(count_Patel, "counts") <- counts
c<-count_Patel@assays@data[1]
count_Patel@assays@data[1]<-count_Patel@assays@data[2]
count_Patel@assays@data[2]<-c
names(count_Patel@assays@data)[1]<-"counts"
names(count_Patel@assays@data)[2]<-"logcounts"

countmatrix_Patel <- assay(count_Patel)
data<-countmatrix_Patel

label <- as.numeric(factor(count_Patel$cell_type1))
n <- length(unique(label))

source('SCSMD.R')
cluster_results <-SCSMD(data, n)


labels_SC3 <- cluster_results[["labels_SC3"]]
labels_Seurat <- cluster_results[["labels_Seurat"]]
labels_SHARP <- cluster_results[["labels_SHARP"]]
labels_cidr <- cluster_results[["labels_cidr"]]
labels_SINCERA <- cluster_results[["labels_SINCERA"]]
labels_Rphenograph <- cluster_results[["labels_Rphenograph"]]
labels_RaceID <- cluster_results[["labels_RaceID"]]

###### calculate ARI values for the 7 methods ######
ARI_sc3 <- adjustedRandIndex(label,as.vector(labels_SC3))
ARI_Seurat <- adjustedRandIndex(label,as.vector(labels_Seurat))
ARI_SHARP <- adjustedRandIndex(label,as.vector(labels_SHARP))
ARI_cidr <- adjustedRandIndex(label,as.vector(labels_cidr))
ARI_SINCERA <- adjustedRandIndex(label,as.vector(labels_SINCERA))
ARI_Rphenograph <- adjustedRandIndex(label,as.vector(labels_Rphenograph))
ARI_RaceID <- adjustedRandIndex(label,as.vector(labels_RaceID))


source("select method.R")
source("the most k.R")
######the number of the cluster######
kr <- cluster_number_range(labels_SC3,
                           labels_Seurat,
                           labels_SHARP,
                           labels_cidr,
                           labels_SINCERA,
                           labels_Rphenograph,
                           labels_RaceID)
m <-min(kr)-1


######Screening clustering method######
select_method<- select_method(labels_SC3,labels_Seurat,labels_SHARP,labels_cidr,labels_SINCERA,labels_Rphenograph,labels_RaceID)
select_method[["rank"]]


cluster_result <- cluster_results[["label_7method"]]

#cluster_result <-cluster_result[-c(1,5,6),,drop=FALSE]

source("the most k.R")
library(SC3)
consensus_result<-consensus_matrix(t(cluster_result))
library('philentropy')
distance_cell<-distance(consensus_result, method = "chebyshev")

library("ctv")
library("kernlab")

k_select <- cluster_number(distance_cell,kr,m)
k_select

k <-kr[which.max(k_select)]
k





ARI_result <- list()
for (i in 1:200){
sc<-specc(distance_cell,centers= k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_patel <- adjustedRandIndex(label,as.vector(label_sc))
ARI_result[[i]] <- ARI_patel
}

ARI_result_vector <- as.vector(unlist(ARI_result))

ARI_patel <-getmode(ARI_result_vector)

ARI_patel

patel_ARI_overall <- c(ARI_sc3, ARI_Seurat, ARI_SHARP, ARI_cidr, ARI_SINCERA, ARI_Rphenograph, ARI_RaceID,ARI_patel)
patel_ARI_overall
patel_overall_rank <-rank(patel_ARI_overall)
patel_overall_rank

as.numeric(table(ARI_result_vector)) 
names(table(ARI_result_vector))      



source("distance_select.R")
dist_select <- distance_select(consensus_result, k)
dist_ARI <- dist_select[["distance_AVG"]]
dist_rank <- dist_select[["distance rank"]]
dist_ARI


sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_patel <- adjustedRandIndex(label,as.vector(label_sc))
ARI_patel


source("datapreprocess.R")
express <-data
express <-datapreprocess(data)


expression_label<-data.frame(label_sc,t(express))
cluster_list <- split(expression_label, expression_label$label_sc)
cluster1 <-data.frame(cluster_list[[1]][,2:5949])
colnames(cluster1)<-colnames(t(express))
cluster2 <-data.frame(cluster_list[[2]][,2:5949])
colnames(cluster2)<-colnames(t(express))

cluster1 <- data.frame(t(cluster_list[[1]][, 2:5949]))
cluster2 <- data.frame(t(cluster_list[[2]][, 2:5949]))

gene_expression_matrix <- cbind(cluster1, cluster2)


gene_p_values <- apply(cluster1, 1, function(gene_expression) {
  t_test_result <- t.test(gene_expression, cluster2[rownames(cluster1),])
  return(t_test_result$p.value)
})


significant_genes1 <- rownames(cluster1)[gene_p_values < 0.00001]




gene_p_values <- apply(cluster2, 1, function(gene_expression) {
  t_test_result <- t.test(gene_expression, cluster1[rownames(cluster2),])
  return(t_test_result$p.value)
})


significant_genes2 <- rownames(cluster2)[gene_p_values < 0.00001]



gene_p_values <- apply(gene_expression_matrix, 1, function(gene_expression) {
  t_test_result <- t.test(gene_expression[1:length(cluster1)], gene_expression[(length(cluster1)+1):length(gene_expression)])
  return(t_test_result$p.value)
})


results_df <- data.frame(Gene = rownames(gene_expression_matrix), PValue = gene_p_values)


results_df <- results_df[order(results_df$PValue), ]

significant_genes <- results_df$Gene[results_df$PValue < 0.00001]

cluster_info <- rep(c("Cluster1", "Cluster2"), c(ncol(cluster1), ncol(cluster2)))

top50_genes <- head(results_df$Gene, 50)


gene_expression_matrix <- as.matrix(gene_expression_matrix)


top50_gene_expression <- gene_expression_matrix[top50_genes, ]

cluster_colors <- c("Cluster1" = "#2E8B57", "Cluster2" = "#EE2C2C")


png(".../heatmap_output3.png", width = 2000, height = 1600, res = 160)

heatmap(top50_gene_expression, Rowv = NA, Colv = NA, scale = "row", dendrogram = "none", labCol = colnames(top50_gene_expression),
        ColSideColors = cluster_colors[cluster_info], col = colorRampPalette(c("white", "yellow", "red"))(n = 200),cexRow = 1.6, cexCol = 0.8)

dev.off()





common_genes <- intersect(intersect(significant_genes, significant_genes1), significant_genes2)


common_genes_1_2 <- intersect( significant_genes1, significant_genes2)


remaining_genes12_cluster1 <- setdiff(significant_genes1, common_genes_1_2)
remaining_genes12_cluster2 <- setdiff(significant_genes2, common_genes_1_2)
remaining_genes_cluster <- setdiff(significant_genes, common_genes_1_2)



significant_genes <- results_df$Gene[results_df$PValue < 0.00001]
remaining_genes_cluster <- setdiff(significant_genes, common_genes_1_2)

selected_rows <- results_df[results_df$Gene %in% remaining_genes_cluster, ]

cluster_info <- rep(c("Cluster1", "Cluster2"), c(ncol(cluster1), ncol(cluster2)))

top50_genes <- head(selected_rows$Gene, 50)
gene_expression_matrix <- as.matrix(gene_expression_matrix)

top50_gene_expression <- gene_expression_matrix[top50_genes, ]

cluster_colors <- c("Cluster1" = "#2E8B57", "Cluster2" = "#EE2C2C")


png(".../heatmap_output process.png", width = 2000, height = 1600, res = 160)

heatmap(top50_gene_expression, Rowv = NA, Colv = NA, scale = "row", dendrogram = "none", labCol = colnames(top50_gene_expression),
        ColSideColors = cluster_colors[cluster_info], col = colorRampPalette(c("white", "yellow", "red"))(n = 200),cexRow = 1.6, cexCol = 0.8)

dev.off()



























































