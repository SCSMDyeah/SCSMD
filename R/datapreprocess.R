datapreprocess<-function(Express) {

  
  percentiles <- apply(Express, 1, function(x) quantile(x, c(0.05, 0.95)))
  
  
  selected_genes <- rownames(Express)[apply(percentiles, 1, function(x) x[1] < x[2])]
  
  expression_threshold <- 1
  high_expression_genes <- rowMeans(Express) >= expression_threshold
  Express <- Express[high_expression_genes, ]

  Express <- log2(Express + 1)
  
  
  return(Express)
}
