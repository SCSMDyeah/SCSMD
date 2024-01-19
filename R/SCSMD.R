Seurat_cluster <- function(data){
  rownames(data) <- paste0("GENE",c(1:dim(data)[1]))
  colnames(data) <- paste0("CELL",c(1:dim(data)[2]))
  seurat_obj <- CreateSeuratObject(counts = data, project = "SEURAT", assay = "RNA")
  all.genes <- rownames(data)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  seurat_obj <- RunPCA(seurat_obj,features = rownames(data))
  seurat_obj <- FindNeighbors(seurat_obj,dims = 1:10)
  seurat_obj <- FindClusters(object = seurat_obj)
  return(as.numeric(Idents(seurat_obj)))
}

CIDR_cluster <- function(data, n){
  sc_cidr <- scDataConstructor(data)
  sc_cidr <- determineDropoutCandidates(sc_cidr)
  sc_cidr <- wThreshold(sc_cidr)
  sc_cidr <- scDissim(sc_cidr)
  sc_cidr <- scPCA(sc_cidr,plotPC = FALSE)
  sc_cidr <- nPC(sc_cidr)
  sc_cidr <- scCluster(sc_cidr, nCluster = n)
  return(sc_cidr@clusters)
}



RaceID_cluster <- function(data,n){
  data <- as.data.frame(data)
  rownames(data) <- 1:dim(data)[1]
  colnames(data) <- 1:dim(data)[2]
  sc <- SCseq(data)
  sc <- filterdata(sc, mintotal=1, minexpr=5, minnumber=1)
  sc <- compdist(sc,metric="pearson")
  sc <- clustexp(sc,sat = FALSE,clustnr=20,bootnr=50, cln=n,rseed=17000)
  labels <- sc@cluster$kpart
  return(labels)
}

SCSMD <- function(data, n){
  library("SC3")
  library("SingleCellExperiment")
  library("Seurat")
  library("SHARP")
  library("cidr")
  library("SINCERA")
  library("Rphenograph")
  library("RaceID")
  library("spatstat.explore")

  #SC3
  sce <- SingleCellExperiment(assays = list(counts = data,logcounts = log2(data+1)))
  rowData(sce)$feature_symbol <- 1:nrow(data)
  dat <- sc3(sce, ks = n,gene_filter = FALSE,n_cores = 5,svm_max = ncol(data))
  labels_SC3 <- as.numeric(dat[[1]])
  print("SC3 is performed")
  
  #Seurat
  labels_Seurat <- Seurat_cluster(data)
  print("Seurat is performed")
  
  #SHARP
  res <- SHARP(data, prep = TRUE, n.cores = 5)
  labels_SHARP <- as.numeric(as.factor(res$pred_clusters))
  print("SHARP is performed")
  
  
  #cidr
  labels_cidr <- CIDR_cluster(data, n)
  print("cidr is performed")
  
  #SINCERA
  dat <- apply(data, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
  dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
  hc <- hclust(dd, method = "average")
  labels_SINCERA <- cutree(hc, k = n)
  print("SINCERA is performed")
  
  #Rphenograph
  Rphenograph_out <- Rphenograph(t(data),k = floor(sqrt(dim(data)[2]))+5)
  labels_Rphenograph <- as.numeric(membership(Rphenograph_out[[2]]))
  print("Rphenograph is performed")
  
  #RaceID
  labels_RaceID <- RaceID_cluster(data,n)
  print("RaceID is performed")
  
  
  label_7method <-rbind(labels_SC3,
                        labels_Seurat,
                        labels_SHARP,
                        labels_cidr,
                        labels_SINCERA,
                        labels_Rphenograph,
                        labels_RaceID)
  result <- list(label_7method = label_7method,
                 labels_SC3 = labels_SC3,
                 labels_Seurat = labels_Seurat,
                 labels_SHARP = labels_SHARP,
                 labels_cidr = labels_cidr,
                 labels_SINCERA = labels_SINCERA,
                 labels_Rphenograph = labels_Rphenograph,
                 labels_RaceID = labels_RaceID)
  return(result)
}
  
  
  
