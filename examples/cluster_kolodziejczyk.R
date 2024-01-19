d <- read.table(".../data/kolodziejczyk/counttable_es.csv")
d <- d[1:(nrow(d) - 5), ]

### ANNOTATIONS
ann <- data.frame(
  cell_type1 = unlist(lapply(strsplit(colnames(d), "_"), "[[", 3)),
  batch = paste(
    unlist(lapply(strsplit(colnames(d), "_"), "[[", 3)),
    unlist(lapply(strsplit(colnames(d), "_"), "[[", 4)),
    sep = "_"
  )
)
rownames(ann) <- colnames(d)
colnames(d) <- rownames(ann)

### SINGLECELLEXPERIMENT
source(".../create_sce.R")
sceset <- create_sce_from_counts(d, ann)

BiocManager::install("biomaRt")
BiocManager::install("SCESet")
BiocManager::install("scater")
library(scater)
library(biomaRt)
library("SCESet")

sceset <- getBMFeatureAnnos(
  sceset, filters="ensembl_gene_id",
  biomart="ensembl", dataset="mmusculus_gene_ensembl",host = "https://asia.ensembl.org/")
# remove features with duplicated names
sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
saveRDS(sceset, "kolodziejczyk.rds")

library(mclust)
library(SingleCellExperiment)

count_kolodziejczyk <-readRDS(".../kolodziejczyk/kolodziejczyk.rds")
countmatrix_kolodziejczyk <- assay(count_kolodziejczyk)
count_kolodziejczyk@rowRanges@elementMetadata@listData[["ensembl_gene_id"]]<-NULL
count_kolodziejczyk@rowRanges@elementMetadata@listData[["mgi_symbol"]]<-NULL
count_kolodziejczyk@rowRanges@elementMetadata@listData[["chromosome_name"]]<-NULL
count_kolodziejczyk@rowRanges@elementMetadata@listData[["gene_biotype"]]<-NULL
count_kolodziejczyk@rowRanges@elementMetadata@listData[["start_position"]]<-NULL
count_kolodziejczyk@rowRanges@elementMetadata@listData[["end_position"]]<-NULL
data<-countmatrix_kolodziejczyk

label <- as.numeric(factor(count_kolodziejczyk$cell_type1))
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
#cluster_result <-cluster_result[-c(1,2,7),,drop=FALSE]


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
k


ARI_result <- list()
label_result <- list()

for (i in 1:200){
sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_kolodziejczyk <- adjustedRandIndex(label,as.vector(label_sc))
ARI_result[[i]] <- ARI_kolodziejczyk
label_result[[i]] <- label_sc
}


ARI_result_vector <- as.vector(unlist(ARI_result))

ARI_kolodziejczyk <-getmode(ARI_result_vector)

ARI_kolodziejczyk

kolodziejczyk_ARI_overall <- c(ARI_sc3, ARI_Seurat, ARI_SHARP, ARI_cidr, ARI_SINCERA, ARI_Rphenograph, ARI_RaceID,ARI_kolodziejczyk)
kolodziejczyk_ARI_overall
kolodziejczyk_overall_rank <-rank(kolodziejczyk_ARI_overall)
kolodziejczyk_overall_rank
as.numeric(table(ARI_result_vector)) 
names(table(ARI_result_vector))      



sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_kolodziejczyk <- adjustedRandIndex(label,as.vector(label_sc))
ARI_kolodziejczyk


source("distance_select.R")
dist_select <- distance_select(consensus_result, k)
dist_ARI <- dist_select[["distance_AVG"]]
dist_rank <- dist_select[["distance rank"]]
dist_ARI


