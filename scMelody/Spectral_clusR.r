###This file contains source code for performing spectral clustering in R####
###Spectral clustering in R###

library(stats)
Affinity <- function(S, n.neighboors) {
  N <- length(S[,1])
  
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } 
  else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A
}


"%^%" <- function(M, power)
  with(eigen(M), vectors %*% (values^power * solve(vectors)))

###Calculating the normalized Laplacian matrix###
Z_fit <- function(A,k){
  D <- diag(apply(A, 1, sum))
  U <- D - A
  L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))
  evL <- eigen(U, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
  return (Z)
}

###Input the similarity matrix and the number of cluster classes to obtain the spectral cluster labels;###
SPEC <- function(S,k){
  N <- length(S[,1])
  
  if (k>= N) {  # fully connected
    A <- S
  } 
  else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:k*2]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  D <- diag(apply(A, 1, sum))
  U <- D - A
  L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))
  evL <- eigen(U, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
  clus_spec <- kmeans(Z,centers = k,nstart = k)
  return(clus_spec$cluster)
}
