### DATA
# mouse2
m2 <- read.csv(".../baron/GSM2230762_mouse2_umifm_counts.csv", header = T)
rownames(m2) <- m2[,1]
labels_m2 <- as.character(m2$assigned_cluster)
m2 <- m2[,4:ncol(m2)]
m2 <- t(m2)


### ANNOTATIONS
# mouse
m_ann2 <- data.frame(
  mouse = c(
    rep(2, length(labels_m2))
  ),
  cell_type1 = c(labels_m2)
)
rownames(m_ann2) <- colnames(m2)

### SINGLECELLEXPERIMENT
source(".../create_sce.R")
m2_sceset <- create_sce_from_counts(m2, m_ann2)
saveRDS(m2_sceset, "baron-mouse2.rds")

source('SCSMD.R')
count_baron_mouse2 <-readRDS(".../baron/baron-mouse2.rds")

data<-assay(count_baron_mouse2)

label <- as.numeric(factor(count_baron_mouse2$cell_type1))
n <- length(unique(label))


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
cluster_result <-cluster_result[-c(1,6,7),,drop=FALSE]

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


ARI_result <- list()
label_result <- list()
for (i in 1:200){
sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_baron_mouse2 <- adjustedRandIndex(label,as.vector(label_sc))
ARI_result[[i]] <- ARI_baron_mouse2
label_result[[i]] <- label_sc
}

ARI_result_vector <- as.vector(unlist(ARI_result))
#source("AIR_stable.R")
ARI_baron_mouse2 <-getmode(ARI_result_vector)


ARI_baron_mouse2
baron_mouse2_ARI_overall <- c(ARI_sc3, ARI_Seurat, ARI_SHARP, ARI_cidr, ARI_SINCERA, ARI_Rphenograph, ARI_RaceID,ARI_baron_mouse2)
baron_mouse2_ARI_overall
baron_mouse2_overall_rank <-rank(baron_mouse2_ARI_overall)
baron_mouse2_overall_rank
as.numeric(table(ARI_result_vector)) 
names(table(ARI_result_vector))      



source("distance_select num.R")
dist_select <- distance_select(consensus_result, k)
dist_ARI <- dist_select[["distance_AVG"]]
dist_ARI
dist_rank<-dist_select[["distance rank"]]
dist_rank



sc<-specc(distance_cell,centers=4, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_baron_mouse2 <- adjustedRandIndex(label,as.vector(label_sc))
ARI_baron_mouse2
