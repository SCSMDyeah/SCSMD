library(mclust)
library(SingleCellExperiment)

furl<-"https://static-content.springer.com/esm/art%3A10.1038%2Fnature13173/MediaObjects/41586_2014_BFnature13173_MOESM31_ESM.txt"
download.file(furl,destfile="./41586_2014_BFnature13173_MOESM31_ESM.txt")
count_treutlein <-read.table("./41586_2014_BFnature13173_MOESM31_ESM.txt")

count_treutlein <- readRDS(".../treutlein/treutlein.rds")
count_treutlein <- t(data.frame(count_treutlein))
presetlabel <- count_treutlein[4:4,2:81]
counts_matrix <- count_treutlein[5:23275,2:81]
counts_matrix <- apply(counts_matrix,2,as.numeric)

colnames(counts_matrix)<-count_treutlein[1,2:81]
rownames(counts_matrix)<-count_treutlein[5:23275,1]
counts_matrix <-data.frame(counts_matrix)
counts_matrix <-as.matrix(counts_matrix)# must be a matrix object!
pretend.cell.labels <- as.factor(presetlabel)  #
pretend.gene.lengths <- as.factor(rownames(counts_matrix))


d_treutlein <- SingleCellExperiment(assays = list(
  normcounts = counts_matrix,
  logcounts = log(counts_matrix +1)
),
colData=DataFrame(cell_type1=pretend.cell.labels),
rowData=DataFrame(feature_symbol=pretend.gene.lengths)
)
data <- assay(d_treutlein)
label <- as.numeric(factor(d_treutlein$cell_type1))
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

cluster_result <-cluster_result[-c(1,5,7),,drop=FALSE]

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
k  #4




ARI_result <- list()
label_result <- list()
for (i in 1:200){
sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_treutlein <- adjustedRandIndex(label,as.vector(label_sc))
ARI_result[[i]] <- ARI_treutlein
label_result[[i]] <- label_sc
}

ARI_result <- as.vector(unlist(ARI_result))
ARI_treutlein <-as.numeric(names(table(ARI_result))[table(ARI_result) == max(table(ARI_result))])
ARI_treutlein
treutlein_ARI_overall <- c(ARI_sc3, ARI_Seurat, ARI_SHARP, ARI_cidr, ARI_SINCERA, ARI_Rphenograph, ARI_RaceID,ARI_treutlein)
treutlein_ARI_overall
treutlein_overall_rank <-rank(treutlein_ARI_overall)
treutlein_overall_rank
as.numeric(table(ARI_result)) 
names(table(ARI_result))     

source("distance_select.R")
dist_select <- distance_select(consensus_result, k)
dist_ARI <- dist_select[["distance_AVG"]]
dist_ARI



sc<-specc(distance_cell,centers=k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
label_sc<-sc@.Data
ARI_treutlein <- adjustedRandIndex(label,as.vector(label_sc))
ARI_treutlein