### DATA
d0 <- read.csv(".../klein/GSM1599494_ES_d0_main.csv", header = FALSE)
d2 <- read.csv(".../klein/GSM1599497_ES_d2_LIFminus.csv", header = FALSE)
d4 <- read.csv(".../klein/GSM1599498_ES_d4_LIFminus.csv", header = FALSE)
d7 <- read.csv(".../klein/GSM1599499_ES_d7_LIFminus.csv", header = FALSE)
d <- cbind(d0, d2[,2:ncol(d2)], d4[,2:ncol(d4)], d7[,2:ncol(d7)])
rownames(d) <- d[,1]
d <- d[,2:ncol(d)]
colnames(d) <- paste0("cell", 1:ncol(d))

### ANNOTATIONS
ann <- data.frame(
  cell_type1 = c(rep("d0", ncol(d0) - 1),
                 rep("d2", ncol(d2) - 1),
                 rep("d4", ncol(d4) - 1),
                 rep("d7", ncol(d7) - 1)))
rownames(ann) <- colnames(d)

### SINGLECELLEXPERIMENT
library(SingleCellExperiment)
source(".../create_sce.R")
sceset <- create_sce_from_counts(d, ann)
saveRDS(sceset, "klein.rds")

library(mclust)
library(SingleCellExperiment)

count_klein <-readRDS(".../data/klein/klein.rds")
data <- assay(count_klein)
label <- as.numeric(factor(count_klein$cell_type1))
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
cluster_result <-cluster_result[-c(1,2,4),,drop=FALSE]

library(SC3)
consensus_result<-consensus_matrix(t(cluster_result))
library('philentropy')
distance_cell<-distance(consensus_result, method = "chebyshev")

library("ctv")
library("kernlab")

source("the most k.R")
k_select <- cluster_number(distance_cell,kr,m)
k_select

k <-kr[which.max(k_select)]
k     #4

k<-6

ARI_result <- list()
label_result <- list()
for (i in 1:50){
sc<-specc(distance_cell,centers= k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_klein <- adjustedRandIndex(label,as.vector(label_sc))
ARI_result[[i]] <- ARI_klein
label_result[[i]] <- label_sc
}

ARI_result_vector <- as.vector(unlist(ARI_result))
#source("AIR_stable.R")
#ARI_baron_mouse1 <- AIR_stable(ARI_result_vector)
ARI_klein <-getmode(ARI_result_vector)
ARI_klein

klein_ARI_overall <- c(ARI_sc3, ARI_Seurat, ARI_SHARP, ARI_cidr, ARI_SINCERA, ARI_Rphenograph, ARI_RaceID,ARI_klein)
klein_ARI_overall
klein_overall_rank <-rank(klein_ARI_overall)
klein_overall_rank

as.numeric(table(ARI_result_vector)) 
names(table(ARI_result_vector))      






sc<-specc(distance_cell,centers= k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_klein <- adjustedRandIndex(label,as.vector(label_sc))
ARI_klein

source("distance_select.R")
dist_select <- distance_select(consensus_result, k)
dist_ARI <- dist_select[["distance_AVG"]]
dist_ARI
dist_rank<-dist_select[["distance rank"]]
dist_rank



library(DESeq2)
library(pheatmap)
data_matrix <- data
cluster_vector<-label_sc
design <- ~ cluster


selected_clusters <- c(1,  6)
subset_colData <- data.frame(cluster = factor(cluster_vector))
subset_colData <- subset(subset_colData, cluster %in% selected_clusters)


selected_clusters <- c(1,  6)
selected_cells <- which(cluster_vector %in% selected_clusters)

subset_data_matrix <- data_matrix[, selected_cells]

design<-~ cluster

subset_colData <- data.frame(cluster = factor(subset_colData$cluster))

if (ncol(subset_data_matrix) != nrow(subset_colData)) {
  stop("Error: The number of columns in subset_data_matrix is not equal to the number of rows in subset_colData.")
}


dds <- DESeqDataSetFromMatrix(countData = subset_data_matrix,
                              colData = subset_colData,
                              design = design)





data_matrix <- data
cluster_vector<-label_sc
design <- ~ cluster

dds <- DESeq(dds)

res <- results(dds)

significant_genes <- subset(res, padj < 0.001)

significant_genes <- significant_genes[order(-abs(significant_genes$log2FoldChange)), ]

top_50_genes <- head(significant_genes, 50)

expression_matrix <- counts(dds)[rownames(top_50_genes), ]

pheatmap(expression_matrix, cluster_rows = TRUE, cluster_cols = TRUE)




