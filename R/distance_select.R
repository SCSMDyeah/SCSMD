distance_select <-function(consensus_result, k){
  library('philentropy')
  
  tolerance <- 1e-3  
  #######chebyshev######
  distance_cell1<-distance(consensus_result, method = "chebyshev")
  
  
    sc<-specc(distance_cell1,centers= k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 300,mod.sample = 0.75, na.action =na.omit)
    label_sc<-sc@.Data
    ARI_chebyshev <- adjustedRandIndex(label,as.vector(label_sc))
    
  
  
  #######euclidean######
  distance_cell2<-distance(consensus_result, method = "euclidean") 
  
    sc<-specc(distance_cell2,centers= k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 300,mod.sample = 0.75, na.action =na.omit)
    label_sc<-sc@.Data
    ARI_euclidean  <- adjustedRandIndex(label,as.vector(label_sc))
   

  
  #######manhattan######
  distance_cell3<-distance(consensus_result, method = "manhattan")
  
  
    sc<-specc(distance_cell3,centers= k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 300,mod.sample = 0.75, na.action =na.omit,tol = tolerance)
    label_sc<-sc@.Data
    ARI_manhattan <- adjustedRandIndex(label,as.vector(label_sc))
    

  
  #############canberra######
  distance_cell4<-distance(consensus_result, method = "canberra") 
  
  
    sc<-specc(distance_cell4,centers= k, kernel ="linear", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 300,mod.sample = 0.75, na.action =na.omit,tol =1e-3 )
    label_sc<-sc@.Data
    ARI_canberra <- adjustedRandIndex(label,as.vector(label_sc))
   

  #######cosine######
  distance_cell5<-distance(consensus_result, method = "cosine")
  
    sc<-specc(distance_cell5,centers= k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 300,mod.sample = 0.75, na.action =na.omit)
    label_sc<-sc@.Data
    ARI_cosine <- adjustedRandIndex(label,as.vector(label_sc))
    
  
  #######sorensen######
  distance_cell6<-distance(consensus_result, method = "sorensen")
  
    sc<-specc(distance_cell6,centers= k, kernel ="rbfdot", kpar = "automatic",nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 300,mod.sample = 0.75, na.action =na.omit)
    label_sc<-sc@.Data
    ARI_sorensen <- adjustedRandIndex(label,as.vector(label_sc))
  
  
  
  
  ARI_overall <- c(ARI_chebyshev, ARI_euclidean, ARI_manhattan, ARI_canberra, ARI_cosine, ARI_sorensen)
  overall_rank <-rank(ARI_overall)  
  
  ARI_overall <-matrix(ARI_overall) 
  rownames(ARI_overall) <-c("chebyshev", "euclidean", "manhattan", "canberra", "cosine", "sorensen")
  overall_rank <-matrix(rank(overall_rank))   
  rownames(overall_rank) <-c("chebyshev", "euclidean", "manhattan", "canberra","cosine", "sorensen")
  
  select <- list(ARI_overall, overall_rank)
  names(select)[1] <- "distance_AVG"
  names(select)[2] <- "distance rank"
  return(select)
  
}
