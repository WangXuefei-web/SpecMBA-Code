#' @import    clue

# The main functions are SpecM(sample, alpha, K) and SpecMBA(sample, alpha, K)
# output NMI, HamAcc, loglik

### Functions for computing covariance matrices

tnsr.expand3 <- function(mat,indices){
  arr <- array(rep(mat,indices[3]),dim=indices)
  tnsr <- as.tensor(arr)
  return(tnsr)
}

meanXXT <- function(arr3d, transpose){ # arrd3: p * q; transpose: logical
  if (!transpose){ # X %*% t(X)
    mysum <- matrix(0, nrow = dim(arr3d)[1], ncol = dim(arr3d)[1])
    L <- dim(arr3d)[3]
    for (l in 1:L) {
      X <- arr3d[,,l]
      mysum <- mysum + X %*% t(X)  
    }
  }
  else{ # t(X) %*% X
    mysum <- matrix(0, nrow = dim(arr3d)[2], ncol = dim(arr3d)[2]) 
    for (l in 1:L) {
      X <- arr3d[,,l]
      mysum <- mysum + t(X) %*% X 
    }
  }

  result <- 1/dim(arr3d)[3] * mysum #  p * p or q * q
  return(result)
} 

### Estimating loading matrix
loading_mat_by_eigenratio <- function(covmat,K) { 

  # spectral decomposition
  eigcov <- eigen(covmat) 
  eigen_values <- eigcov$values # has sorted in decreasing order
  
  # eigenvalue ratio test for the number of factors
  # k_max <- floor((dim(covmat)[1])/3)
  # k_min <- min(K,dim(covmat)[1])
  # eigen_ratios <- eigen_values[-length(eigen_values)] / eigen_values[-1]  
  # eigen_ratios <- eigen_ratios[1:k_max]
  # max_ratio_index <- max(k_min,(which.max(eigen_ratios[2:k_max])+1))

  # estimating loading matrix
  # eigen_vectors <- eigcov$vectors[, 1:max_ratio_index]
  U <- eigcov$vectors[, 1:K]
  
  return(U)
}  

### Maximization likelihood method for selecting optimal alpha
compute_loglik <- function(A,label.hat){
  # calculate count statistics
  if (class(A)=="Tensor"){
    A <- A@data
  }
  N <- dim(A)[1]
  L <- dim(A)[3]
  K <- length(unique(label.hat))
  K.I <- diag(K)
  label.OneHot <- as.matrix(K.I[label.hat, ]) # N*K
  n.k <- as.matrix(apply(label.OneHot,2,sum)) # 1*K, size of communities
  
  if (K==1){
    diag.minus <- n.k
  }else{
    diag.minus <- diag(as.vector(n.k))
  }
  
  n.kk <- n.k %*% t(n.k) - diag.minus # number of possible edges
  Pi <- n.k/N
  
  m.kk <- array(NA, dim=c(K,K,L))
  P <- array(NA, dim=c(K,K,L))
  for (l in 1:L){
    m.mat <- t(label.OneHot)%*%A[,,l]%*%label.OneHot # K*K
    diag(m.mat) <- diag(m.mat)/2
    m.kk[,,l] <- m.mat
    P[,,l] <- m.kk[,,l]/n.kk
  } 
  
  # calculate MLE estimator
  n.kkl <- array(rep(n.kk,L),dim = c(K,K,L))
  l0.arr <- m.kk*(log(P+1e-30))+(n.kkl-m.kk)*(log(1-P+1e-30)) # K*K
  l.A <- sum(sum(sum(l0.arr)))
  l.Z <- sum(n.k*log(Pi))
  loglik <- l.A + l.Z
  return(loglik)
}

## Calculate NMI between two community assignments

nmi<-function(clus,comm)
{
  n=length(clus)
  entclus=sum((table(clus)/n)*log(table(clus)/n))
  entcomm=sum((table(comm)/n)*log(table(comm)/n))
  mi=sum((table(clus,comm)/n)*log(n*table(clus,comm)/(table(clus)%*%t(table(comm)))),na.rm=TRUE)
  nmi=-mi/((entclus+entcomm)/2)
  return(nmi)
}

### MCrate calculates misclassification rate between two community assignments
pMatrix.min <- function(A, B) { 
  n <- nrow(A) 
  D <- matrix(NA, n, n) 
  for (i in 1:n) { 
    for (j in 1:n) { 
      D[j, i] <- (sum((B[j, ] - A[i, ])^2)) 
    } } 
  vec <- c(solve_LSAP(D)) 
  list(A=A[vec,], pvec=vec) 
}
# estimate, true 
MCrate <- function (C, E){
  E = factor(E); C = factor(C)
  N = length(E); k = length(levels(C))
  E = factor(x = E, levels = 1:k)
  A = table(E,C)
  X <- pMatrix.min(A,diag(1,k)) 
  A = X$A
  N_err = sum(A) - sum(diag(A))
  P_err = N_err/N
  return(1-P_err)
}


### Function for applying SpecM
SpecM <- function(sample, alpha, K){ 
  # output NMI, HamAcc, loglik
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # compute covariance-type matrix
  A.bar <- as.matrix(modeMean(A,3,drop=TRUE)@data)
  A.bar.expand <- tnsr.expand3(mat=A.bar,indices=c(N,N,L))
  A.diff.tnsr <- A - A.bar.expand
  A.diff <- as.array(A.diff.tnsr@data) # N*N*L

  my.meanXXT <- meanXXT(A.diff,transpose=FALSE)
  Mu <- 1/N^2*((1+alpha)*A.bar%*%t(A.bar)+my.meanXXT)
  # Mv <- 1/N^2*((1+alpha)*t(A.bar)%*%A.bar+meanXXT(A.diff,transpose=TRUE)) # is same as Mu if A is directed
  
  
  # estimating U
  U <- loading_mat_by_eigenratio(Mu,K)
  
  # community detection 
  get.est.flag <- FALSE
  temp <- list()
  try.temp <- 1
  while(!get.est.flag){
    temp <- try(kmeans(U, K, nstart = 10),silent=TRUE)
    if('try-error' %in% class(temp)) # judge weather error occurs
    {
      try.temp <- try.temp + 1
      if (try.temp == 3)
      {
        stop("There seems to be an issue with the initialization within the kmeans step. Please run the algorithm again.") 
      }
      next
    }else{
      get.est.flag <- TRUE
    }
  }
  kmeans.re <- temp 
  label.hat <- kmeans.re$cluster
  
  # metric
  NMI <- nmi(label,label.hat)
  HamAcc <- MCrate(label, label.hat)
  loglik <- compute_loglik(A,label.hat)
  
  output <- list("label.hat"=label.hat, "NMI"=NMI, "HamAcc"=HamAcc, "loglik"=loglik)
  return(output)
}

SpecMBA <- function(sample, alpha, K){ 
  # output NMI, HamAcc, loglik
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # computations based on adjacency tensor
  A.bar <- as.matrix(modeMean(A,3,drop=TRUE)@data)
  A.bar.expand <- tnsr.expand3(mat=A.bar,indices=c(N,N,L))
  A.diff.tnsr <- A - A.bar.expand
  A.diff <- as.array(A.diff.tnsr@data) # N*N*L

  # debias
  d.bar <- colSums(A.bar)
  D.bar <- diag(d.bar)
  D.bias <- (1+alpha/L)*D.bar
  
  # compute covariance-type matrix
  my.meanXXT <- meanXXT(A.diff,transpose=FALSE)
  Mu <- 1/N^2*((1+alpha)*A.bar%*%t(A.bar)+my.meanXXT-D.bias)
  # Mv <- 1/N^2*((1+alpha)*t(A.bar)%*%A.bar+meanXXT(A.diff,transpose=TRUE)) # is same as Mu if A is directed
  
  # estimating U
  U <- loading_mat_by_eigenratio(Mu,K)
  
  # community detection 
  get.est.flag <- FALSE
  temp <- list()
  try.temp <- 1
  while(!get.est.flag){
    temp <- try(kmeans(U, K, nstart = 10),silent=TRUE)
    if('try-error' %in% class(temp)) # judge weather error occurs
    {
      try.temp <- try.temp + 1
      if (try.temp == 3)
      {
        stop("There seems to be an issue with the initialization within the kmeans step. Please run the algorithm again.") 
      }
      next
    }else{
      get.est.flag <- TRUE
    }
  }
  kmeans.re <- temp 
  label.hat <- kmeans.re$cluster
  
  # metric
  NMI <- nmi(label,label.hat)
  HamAcc <- MCrate(label, label.hat)
  loglik <- compute_loglik(A,label.hat)
  
  output <- list("label.hat"=label.hat,"NMI"=NMI, "HamAcc"=HamAcc, "loglik"=loglik)
  return(output)
}

### Application sample ###
# For real data without label
SpecMBA.A <- function(A, alpha, K){ 
  # output loglik
  # A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  # label <- sample$global_membership
  
  # computations based on adjacency tensor
  A.bar <- as.matrix(modeMean(A,3,drop=TRUE)@data)
  A.bar.expand <- tnsr.expand3(mat=A.bar,indices=c(N,N,L))
  A.diff.tnsr <- A - A.bar.expand
  A.diff <- as.array(A.diff.tnsr@data) # N*N*L
  
  # debias
  d.bar <- colSums(A.bar)
  D.bar <- diag(d.bar)
  D.bias <- (1+alpha/L)*D.bar
  
  # compute covariance-type matrix
  my.meanXXT <- meanXXT(A.diff,transpose=FALSE)
  Mu <- 1/N^2*((1+alpha)*A.bar%*%t(A.bar)+my.meanXXT-D.bias)
  # Mv <- 1/N^2*((1+alpha)*t(A.bar)%*%A.bar+meanXXT(A.diff,transpose=TRUE)) # is same as Mu if A is directed
  
  # estimating U
  U <- loading_mat_by_eigenratio(Mu,K)
  
  # community detection 
  get.est.flag <- FALSE
  temp <- list()
  try.temp <- 1
  while(!get.est.flag){
    temp <- try(kmeans(U, K, nstart = 10),silent=TRUE)
    if('try-error' %in% class(temp)) # judge weather error occurs
    {
      try.temp <- try.temp + 1
      if (try.temp == 3)
      {
        stop("There seems to be an issue with the initialization within the kmeans step. Please run the algorithm again.") 
      }
      next
    }else{
      get.est.flag <- TRUE
    }
  }
  kmeans.re <- temp 
  label.hat <- kmeans.re$cluster
  
  # metric
  # NMI <- nmi(label,label.hat)
  # HamAcc <- MCrate(label, label.hat)
  loglik <- compute_loglik(A,label.hat)
  
  output <- list("label.hat"=label.hat,"loglik"=loglik,"M"=Mu)
  return(output)
}

SpecMBA.app <- function(A,alpha.seq=seq(-1,5,0.2),K){
  # output alpha, label.hat
  result.seq <- array(NA, dim=c(length(alpha.seq),2)) # loglik alpha
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecMBA <- SpecMBA.A(A, alpha, K)
    result.seq[alpha.i,] <- c(re.SpecMBA$loglik,alpha)
    alpha.i <- alpha.i+1
  }
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('loglik','alpha')
  max.id <- which.max(result.seq$loglik)

  alpha.opt <- result.seq[max.id,'alpha']
  label.hat <- SpecMBA.A(A, alpha.opt, K)$label.hat
  
  output <- list("alpha"=alpha.opt,"label.hat"=label.hat)
  return(output)
}