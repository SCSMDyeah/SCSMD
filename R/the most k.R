######select the number of the cluster range

cluster_number_range <- function(labels_SC3,labels_Seurat,labels_SHARP,labels_cidr,labels_SINCERA,labels_Rphenograph,labels_RaceID){

######the number of the cluster######
n1 <- length(unique(as.vector(labels_SC3)))
n2 <- length(unique(as.vector(labels_Seurat)))
n3 <- length(unique(as.vector(labels_SHARP)))
n4 <- length(unique(as.vector(labels_cidr)))
n5 <- length(unique(as.vector(labels_SINCERA)))
n6 <- length(unique(as.vector(labels_Rphenograph)))
n7 <- length(unique(as.vector(labels_RaceID)))
N <-c(n1, n2, n3, n4, n5, n6, n7)
N
n_7AVG <- sum(N)/7
n_7AVG
n_3AVG <- sum(n2, n3, n6)/3
n_3AVG
min(N,n_7AVG,n_3AVG)
max(N,n_7AVG,n_3AVG)
k <- c(min(N,n_7AVG,n_3AVG):max(N,n_7AVG,n_3AVG))
return(k)
}




getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



cluster_number <- function(distance_cell,kr,m){
   k_most <-list()
   for(j in kr){

     #ARI_result <- list()
     #label_result <- list()
        #for (i in 1:200){
          sc<-specc(distance_cell,centers= j, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 200,mod.sample = 0.75, na.action =na.omit)
          label_sc<-sc@.Data
          #ARI <- adjustedRandIndex(label,as.vector(label_sc))
          #ARI_result[[i]] <- ARI
          #label_result[[i]] <- label_sc
         #}

      #ARI_result_vector <- as.vector(unlist(ARI_result))

      #ARI <-getmode(ARI_result_vector)

      #a<-as.matrix(which(ARI_result == ARI))

      #label_sc <- label_result[[a[1,1]]]


    #ARI_SC3_sc <- adjustedRandIndex(label_sc,as.vector(labels_SC3))
    ARI_Seurat_sc <- adjustedRandIndex(label_sc,as.vector(labels_Seurat))
    ARI_SHARP_sc <- adjustedRandIndex(label_sc,as.vector(labels_SHARP))
    ARI_cidr_sc <- adjustedRandIndex(label_sc,as.vector(labels_cidr))
    #ARI_SINCERA_sc <- adjustedRandIndex(label_sc,as.vector(labels_SINCERA))
    ARI_Rphenograph_sc <- adjustedRandIndex(label_sc,as.vector(labels_Rphenograph))
    ARI_RaceID_sc <- adjustedRandIndex(label_sc,as.vector(labels_RaceID))
    b <- c(ARI_Seurat_sc, ARI_SHARP_sc,ARI_cidr_sc, ARI_Rphenograph_sc, ARI_RaceID_sc)
    #b <- c(ARI_SC3_sc,ARI_Seurat_sc, ARI_SHARP_sc,ARI_cidr_sc, ARI_SINCERA_sc, ARI_Rphenograph_sc, ARI_RaceID_sc)
    k_most[[j-m]] <-sum(b)/5
    }

  return(k_most)
}









