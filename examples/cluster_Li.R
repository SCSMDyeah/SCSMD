### DATA
d <- read.csv(".../Li/GSE81861_Cell_Line_COUNT.csv")


### ANNOTATIONS
genes <- unlist(lapply(strsplit(as.character(d[,1]), "_"), "[[", 2))
d <- d[!duplicated(genes), ]
rownames(d) <- genes[!duplicated(genes)]
d <- d[,2:ncol(d)]
# metadata
ann <- data.frame(cell_type1 = unlist(lapply(strsplit(colnames(d), "__"), "[[", 2)))
rownames(ann) <- colnames(d)

### SINGLECELLEXPERIMENT
source(".../R/create_sce.R")
sceset <- create_sce_from_counts(d, ann)
saveRDS(sceset, "li.rds")

library(mclust)
library(SingleCellExperiment)
count_li <-readRDS(".../Li/li.rds")
write.table(data,".../Lili.csv",row.names=TRUE,col.names=TRUE,sep=",")


data<-assay(count_li)

label <- as.numeric(factor(count_li$cell_type1))
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


ARI_sc3 <- adjustedRandIndex(label,as.vector(labels_SC3))
ARI_Seurat <- adjustedRandIndex(label,as.vector(labels_Seurat))
ARI_SHARP <- adjustedRandIndex(label,as.vector(labels_SHARP))
ARI_cidr <- adjustedRandIndex(label,as.vector(labels_cidr))
ARI_SINCERA <- adjustedRandIndex(label,as.vector(labels_SINCERA))
ARI_Rphenograph <- adjustedRandIndex(label,as.vector(labels_Rphenograph))
ARI_RaceID <- adjustedRandIndex(label,as.vector(labels_RaceID))


source("select method.R")
source("the most k.R")

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
cluster_result <-cluster_result[-c(2,4,6),,drop=FALSE]

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
k                #9


ARI_result <- list()
label_result <- list()
for (i in 1:200){
sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_Li <- adjustedRandIndex(label,as.vector(label_sc))
ARI_result[[i]] <- ARI_Li
label_result[[i]] <- label_sc
}

#source("AIR_stable.R")

ARI_result_vector <- as.vector(unlist(ARI_result))
#ARI_Li <- AIR_stable(ARI_result_vector)
ARI_Li <-getmode(ARI_result_vector)

ARI_Li

Li_ARI_overall <- c(ARI_sc3, ARI_Seurat, ARI_SHARP, ARI_cidr, ARI_SINCERA, ARI_Rphenograph, ARI_RaceID,ARI_Li)
Li_ARI_overall
Li_overall_rank <-rank(Li_ARI_overall)
Li_overall_rank
as.numeric(table(ARI_result_vector)) 
names(table(ARI_result_vector))     

#8 5 1 6 2 3 4 7 centers=9 

sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_Li <- adjustedRandIndex(label,as.vector(label_sc))
ARI_Li

source("distance_select.R")
dist_select <- distance_select(consensus_result, k)
dist_ARI <- dist_select[["distance_AVG"]]
dist_ARI
dist_rank <- dist_select[["distance rank"]]
dist_rank




