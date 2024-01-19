library(readxl)
library("ctv")
library("kernlab")

A <- cluster_result
ss <- 0
dd <- 0
ds <- 0
sd <- 0

row <- nrow(A)
col <- ncol(A)

similarity_matrix <- matrix(NA, nrow = row, ncol = row)

for (i in 1:(row - 1)) {
  for (m in (i + 1):row) {
    ss <- 0
    dd <- 0
    ds <- 0
    sd <- 0
    
    for (j in 1:(col - 1)) {
      for (k in (j + 1):col) {
        if (A[i, j] == A[i, k] & A[m, j] == A[m, k]) {
          ss <- ss + 1
        } else if (A[i, j] != A[i, k] & A[m, j] != A[m, k]) {
          dd <- dd + 1
        } else if (A[i, j] != A[i, k] & A[m, j] == A[m, k]) {
          ds <- ds + 1
        } else {
          sd <- sd + 1
        }
      }
    }
    
    denominator <- ss + dd + ds + sd
    
    if (denominator != 0) {
      result <- (ss + dd) / denominator
    } else {
      result <- NaN
    }
    
    similarity_matrix[i, m] <- result
    similarity_matrix[m, i] <- result
  }
}


similarity_matrix[is.na(similarity_matrix)] <- 1
colnames(similarity_matrix)<-c(1,2,3,4,5,6,7)
print(similarity_matrix)
#distance_cell<-distance(similarity_matrix, method = "chebyshev")
sc <- specc(similarity_matrix, centers=4, nystrom.red = FALSE, nystrom.sample =dim(x)[1]/6, iterations = 100)
sc




