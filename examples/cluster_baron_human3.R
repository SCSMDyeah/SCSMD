### DATA
# human3
h3 <- read.csv(".../baron/GSM2230759_human3_umifm_counts.csv", header = T)
rownames(h3) <- h3[,1]
labels_h3 <- as.character(h3$assigned_cluster)
h3 <- h3[,4:ncol(h3)]
h3 <- t(h3)


### ANNOTATIONS
# human
h_ann3 <- data.frame(
  human = c(
    rep(3, length(labels_h3))
  ),
  cell_type1 = c(labels_h3))
rownames(h_ann3) <- colnames(h3)


### SINGLECELLEXPERIMENT
source(".../create_sce.R")
h3_sceset <- create_sce_from_counts(h3, h_ann3)

saveRDS(h3_sceset, "baron-human3.rds")

library(SingleCellExperiment)
count_baron_human3 <-readRDS(".../baron/baron-human3.rds")

data<-assay(count_baron_human3)

label <- as.numeric(factor(count_baron_human3$cell_type1))
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

cluster_result <-cluster_result[-c(1,2,7),,drop=FALSE]


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
k<-7

ARI_result <- list()
label_result <- list()
for (i in 1:50){
sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_baron_human3 <- adjustedRandIndex(label,as.vector(label_sc))
ARI_result[[i]] <- ARI_baron_human3
label_result[[i]] <- label_sc
}

ARI_result_vector <- as.vector(unlist(ARI_result))

ARI_baron_human3 <-getmode(ARI_result_vector)

ARI_baron_human3
baron_human3_ARI_overall <- c(ARI_sc3, ARI_Seurat, ARI_SHARP, ARI_cidr, ARI_SINCERA, ARI_Rphenograph, ARI_RaceID,ARI_baron_human3)
baron_human3_ARI_overall
baron_human3_overall_rank <-rank(baron_human3_ARI_overall)
baron_human3_overall_rank
as.numeric(table(ARI_result_vector)) 
names(table(ARI_result_vector))      




sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_baron_human3 <- adjustedRandIndex(label,as.vector(label_sc))
ARI_baron_human3

source("distance_select.R")
dist_select <- distance_select(consensus_result, k)
dist_ARI <- dist_select[["distance_AVG"]]
dist_ARI
dist_rank<-dist_select[["distance rank"]]
dist_rank






