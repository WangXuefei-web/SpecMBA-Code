### Functions for competitive methods ###
#' @import     rMultiNet # for TWIST
#' @import     rTensor # for HOOI
#' 
### Method: SpecM ###
SpecM.fun <- function(sample,alpha.seq,K){
  start <- Sys.time()
  result.seq <- array(NA, dim=c(length(alpha.seq),4)) # loglik alpha NMI HamAcc
  lb.list <- list()
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecM <- SpecM(sample, alpha, K)
    result.seq[alpha.i,] <- c(re.SpecM$loglik,alpha,re.SpecM$NMI,re.SpecM$HamAcc)
    lb.list[[alpha.i]] <- list('alpha'=alpha,'label.hat'=re.SpecM$label.hat)
    alpha.i <- alpha.i+1
  }
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('loglik','alpha','NMI','HamAcc')
  max.id <- which.max(result.seq$loglik)

  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("alpha"=result.seq[max.id,'alpha'],"NMI"=result.seq[max.id,'NMI'],
                 "HamAcc"=result.seq[max.id,'HamAcc'],"label.hat"=lb.list[[max.id]]$label.hat,"time"=as.numeric(time))
  return(output)
}

### Method: SpecMBA ###
SpecMBA.fun <- function(sample,alpha.seq,K){
  start <- Sys.time()
  result.seq <- array(NA, dim=c(length(alpha.seq),4)) # loglik alpha NMI HamAcc
  lb.list <- list()
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecMBA <- SpecMBA(sample, alpha, K)
    result.seq[alpha.i,] <- c(re.SpecMBA$loglik,alpha,re.SpecMBA$NMI,re.SpecMBA$HamAcc)
    lb.list[[alpha.i]] <- list('alpha'=alpha,'label.hat'=re.SpecMBA$label.hat)
    alpha.i <- alpha.i+1
  }
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('loglik','alpha','NMI','HamAcc')
  max.id <- which.max(result.seq$loglik)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("alpha"=result.seq[max.id,'alpha'],"NMI"=result.seq[max.id,'NMI'], 
                 "HamAcc"=result.seq[max.id,'HamAcc'],"label.hat"=lb.list[[max.id]]$label.hat,"time"=as.numeric(time))
  return(output)
}

SpecM.fun <- function(sample,alpha.seq,K){
  start <- Sys.time()
  result.seq <- array(NA, dim=c(length(alpha.seq),4)) # loglik alpha NMI HamAcc
  lb.list <- list()
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecM <- SpecM(sample, alpha, K)
    result.seq[alpha.i,] <- c(re.SpecM$loglik,alpha,re.SpecM$NMI,re.SpecM$HamAcc)
    lb.list[[alpha.i]] <- list('alpha'=alpha,'label.hat'=re.SpecM$label.hat)
    alpha.i <- alpha.i+1
  }
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('loglik','alpha','NMI','HamAcc')
  max.id <- which.max(result.seq$loglik)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("alpha"=result.seq[max.id,'alpha'],"NMI"=result.seq[max.id,'NMI'], 
                 "HamAcc"=result.seq[max.id,'HamAcc'],"label.hat"=lb.list[[max.id]]$label.hat,"time"=as.numeric(time))
  return(output)
}

SpecMBA.fast.fun <- function(sample,alpha.seq,K){
  start <- Sys.time()
  result.seq <- array(NA, dim=c(length(alpha.seq),4)) # loglik alpha NMI HamAcc
  lb.list <- list()
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecMBA <- SpecMBA(sample, alpha, K)
    result.seq[alpha.i,] <- c(re.SpecMBA$loglik,alpha,re.SpecMBA$NMI,re.SpecMBA$HamAcc)
    lb.list[[alpha.i]] <- list('alpha'=alpha,'label.hat'=re.SpecMBA$label.hat)
    alpha.i <- alpha.i+1
  }
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('loglik','alpha','NMI','HamAcc')
  max.id <- which.max(result.seq$loglik)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("alpha"=result.seq[max.id,'alpha'],"NMI"=result.seq[max.id,'NMI'], 
                 "HamAcc"=result.seq[max.id,'HamAcc'],"label.hat"=lb.list[[max.id]]$label.hat,"time"=as.numeric(time))
  return(output)
}

### Method: mean adjacency matrix ###
computeLaplacian <- function(Adj){
  d <- apply(Adj,1,sum)
  d[d==0] <- 1e-10
  D <- d^(0.5)
  D.inv <- D^(-1)
  G <- diag(D)
  G.inv <- diag(D.inv)
  L <- G.inv%*%Adj%*%G.inv
  return(L)
}

eigenratio_test <- function(mat){
  k_max <- floor((dim(mat)[1])/3)
  k_min <- min(3,dim(mat)[1])
  
  # spectral decomposition
  eigcov <- eigen(mat) 
  eigen_values <- eigcov$values # has sorted in decreasing order
  
  # eigenvalue ratio test for the number of factors
  eigen_ratios <- eigen_values[-length(eigen_values)] / eigen_values[-1]  
  eigen_ratios <- eigen_ratios[1:k_max]
  max_ratio_index <- max(k_min,which.max(eigen_ratios))
  
  # estimating loading matrix
  eigen_vectors <- eigcov$vectors[, 1:max_ratio_index]
  return(eigen_vectors)
}

MeanAdj.fun<-function(sample,K)
{
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # compute mean adjacency matrix and its Laplacian
  A.bar <- as.matrix(modeMean(A,3,drop=TRUE)@data)
  Lap <- computeLaplacian(A.bar)
  
  eigenvectors <- eigenratio_test(Lap)
  
  # community detection 
  get.est.flag <- FALSE
  temp <- list()
  try.temp <- 1
  while(!get.est.flag){
    temp <- try(kmeans(eigenvectors, K, nstart = 10),silent=TRUE)
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
  
  NMI <- nmi(label,label.hat)
  HamAcc <- MCrate(label, label.hat)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  
  return(output)
}

### Method: OLMF orthogonal linked matrix factorization ###
## The objective function 

lmffunctiono<-function(param,laplist,n,k){
  M=length(laplist)
  ustar<-matrix(param[1:(n*k)],n,k)
  lambda<-lapply(1:M,function(m){return(matrix(param[(n*k+(m-1)*k^2+1):(n*k+m*k^2)],k,k))})
  objloop<- sum(unlist(lapply(1:M,function(m){
    specobj<-norm(laplist[[m]]-ustar%*%lambda[[m]]%*%t(ustar),type="F")^2
    return(specobj)
  })))
  obj=objloop
  return(obj)
}


##  The gradients

lmfdero<-function(param,laplist,n,k){
  M=length(laplist)
  ustar<-matrix(param[1:(n*k)],n,k)
  lambda<-lapply(1:M,function(m){return(matrix(param[(n*k+(m-1)*k^2+1):(n*k+m*k^2)],k,k))})
  derlist1<-lapply(1:M,function(m){
    specobj= -(diag(n)-ustar%*%t(ustar))%*%laplist[[m]]%*%ustar%*%lambda[[m]]
    return(specobj)
  })
  derlist2<-lapply(1:M,function(m){
    specobj= -t(ustar)%*%(laplist[[m]]-ustar%*%lambda[[m]]%*%t(ustar))%*%ustar
    return(specobj)
  })
  der1<-Reduce("+",derlist1)
  der2<-unlist(derlist2)
  return(c(as.vector(der1),as.vector(der2)))
}


## The main function with BFGS optimization

OLMF.fun<-function(sample,K){
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # compute the Laplacian of each layer
  Laplist <- lapply(1:L,function(l){return(computeLaplacian(A@data[,,l]))})
  
  # Initialize with mean Laplacian
  Lap.bar <- Reduce('+',Laplist)
  spectra <- eigen(Lap.bar)
  ustar <- spectra$vectors[,1:K]
  lambda <- lapply(1:L,function(l){return(diag(spectra$values[1:K]))})
  
  # Optimization
  param <- c(as.vector(ustar),as.vector(unlist(lambda)))
  optimized <- optim(par=param,fn=lmffunctiono,gr=lmfdero,method="BFGS",control=list(reltol=0.0001,maxit=200),laplist=Laplist,n=N,k=K)
  param <- optimized$par
  
  # Community detection
  ustar <- matrix(param[1:(N*K)],N,K)
  lambda <- lapply(1:L,function(l){return(matrix(param[(N*K+(l-1)*K^2+1):(N*K+l*K^2)],K,K))})
  
  specstar <- kmeans(ustar,K)
  label.hat <- specstar$cluster
  
  # Compute metrics
  NMI <- nmi(label,label.hat)
  HamAcc <- MCrate(label, label.hat)
  
  end <- Sys.time()
  time <- end - start
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  
  return(output)
  
}

### Method: SoSBA ###
SoSBA.fun <- function(sample,K){
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # computations based on adjacency tensor
  A.l <- as.array(A@data) # N*N*L
  
  # debias
  A.bar <- as.matrix(modeMean(A,3,drop=TRUE)@data)
  d.bar <- colSums(A.bar)
  D.bar <- diag(d.bar)
  
  # compute covariance-type matrix
  AA <- meanXXT(A.l,transpose=FALSE)
  S <- L * (AA - D.bar)
  
  # estimating U
  U <- loading_mat_by_eigenratio(S,K)
  
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
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  return(output)
}

SoS.fun <- function(sample,K){
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # computations based on adjacency tensor
  A.l <- as.array(A@data) # N*N*L
  
  # debias
  A.bar <- as.matrix(modeMean(A,3,drop=TRUE)@data)
  
  # compute covariance-type matrix
  AA <- meanXXT(A.l,transpose=FALSE)
  S <- L * (AA)
  
  # estimating U
  U <- loading_mat_by_eigenratio(S,K)
  
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
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  return(output)
}


### Method: TWIST ###
InitializationMMSBM<-function(tnsr, ranks, estM=FALSE)
{
  num_modes <- tnsr@num_modes
  U_list <- vector("list", num_modes)
  temp_mat <-matrix(0,ncol=tnsr@modes[1], nrow=tnsr@modes[2])
  for(i in 1:tnsr@modes[3])
  {
    temp_mat=temp_mat+tnsr@data[,,i]
  }
  U_list[[1]] <- eigen(temp_mat,symmetric = T)$vector[,c(1:ranks[1])]
  U_list[[2]] <- U_list[[1]]
  outer_production=NULL
  for(i in 1:dim(U_list[[1]])[1])
  {
    row_matrix=NULL
    for (j in 1:dim(U_list[[1]])[2])
    {
      temp=U_list[[1]][i,j]*U_list[[2]]
      row_matrix=cbind(row_matrix,temp)
    }
    outer_production=rbind(outer_production,row_matrix)
  }
  temp_mat <- rs_unfold(tnsr, m = 3)@data %*% outer_production

  if (!estM) {# M is given
    U_list[[3]] <- svd(temp_mat, nu = ranks[3])$u
  }else{ # Estimate M
    svd_mode3 <- svd(temp_mat, nu = ranks[3])
    #eigenvalue ratio estimator
    d <- svd_mode3$d
    M <- which.max(d[1:length(d)-1]/d[2:length(d)])
    U_list[[3]] <- as.matrix(svd(temp_mat, nu = ranks[3])$u[,1:M])
  }
  return(U_list)
}

norm_vec <- function(x) sqrt(sum(x^2))
reg_vec <- function(x,delta) min(delta,norm_vec(x))/norm_vec(x)*x

PowerIteration<- function(tnsr, ranks=NULL, type="TWIST", U_0_list, delta1=1000, delta2=1000, max_iter = 20, tol = 1e-04)
{
  stopifnot(is(tnsr, "Tensor"))
  if (is.null(ranks))
    stop("ranks must be specified")
  if (sum(ranks > tnsr@modes) != 0)
    stop("ranks must be smaller than the corresponding mode")
  if (sum(ranks <= 0) != 0)
    stop("ranks must be positive")
  if(type == "TWIST")
  {
    num_modes <- tnsr@num_modes
    U_list <- U_0_list
    tnsr_norm <- fnorm(tnsr)
    curr_iter <- 1
    converged <- FALSE
    fnorm_resid <- rep(0, max_iter)
    CHECK_CONV <- function(Z, U_list) {
      est <- ttl(Z, U_list, ms = 1:num_modes)
      curr_resid <- fnorm(tnsr - est)
      fnorm_resid[curr_iter] <<- curr_resid
      if (curr_iter == 1)
        return(FALSE)
      if (abs(curr_resid - fnorm_resid[curr_iter - 1])/tnsr_norm <
          tol)
        return(TRUE)
      else {
        return(FALSE)
      }
    }
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
    while ((curr_iter < max_iter) && (!converged)) {
      cat("iteration", curr_iter, "\n")
      setTxtProgressBar(pb, curr_iter)
      modes <- tnsr@modes
      modes_seq <- 1:num_modes
      
      ##Regularization
      U_list_reg = U_list
      for(m in modes_seq)
      {
        if(m == 1 | m == 2)
        {
          U_list_reg[[m]] = as.matrix(apply(U_list_reg[[m]],1,reg_vec, delta=delta1))
          if(ranks[m]!=1){
            U_list_reg[[m]] = t(U_list_reg[[m]])
          }
        }
        if(m == 3)
        {
          U_list_reg[[m]] = as.matrix(apply(U_list_reg[[m]],1,reg_vec, delta=delta2))
          if(ranks[m]!=1){
            U_list_reg[[m]] = t(U_list_reg[[m]])
          }
        }
      }
      
      ##Iterate
      
      for (m in modes_seq) {
        mat_list <- lapply(U_list_reg[-m], t)
        X <- ttl(tnsr, mat_list, ms = modes_seq[-m])
        U_list[[m]] <- svd(rs_unfold(X, m = m)@data, nu = ranks[m])$u
      }
      
      Z <- ttm(X, mat = t(U_list[[num_modes]]), m = num_modes)
      
      if (CHECK_CONV(Z, U_list)) {
        converged <- TRUE
        setTxtProgressBar(pb, max_iter)
      }else {
        curr_iter <- curr_iter + 1
      }
    }
    close(pb)
    fnorm_resid <- fnorm_resid[fnorm_resid != 0]
    norm_percent <- (1 - (tail(fnorm_resid, 1)/tnsr_norm))*100
    est <- ttl(Z, U_list, ms = 1:num_modes)
    invisible(list(Z = Z, U = U_list, conv = converged, est = est,
                   norm_percent = norm_percent, fnorm_resid = tail(fnorm_resid,
                                                                   1), all_resids = fnorm_resid))
    network_embedding <- U_list[[3]]
    node_embedding <- U_list[[1]]
    output <- list("Z"=Z, "network_embedding"=network_embedding,"node_embedding"=node_embedding)
    return(output)
  }
  if(type == "TUCKER")
  {
    decomp=tucker(arrT,ranks,max_iter = 10000,tol=1e-05)
    nodes_embedding_Our=decomp[["U"]][[1]] #nodes' embedding
    network_embedding=decomp[["U"]][[3]] #layers' embedding
    Z = decomp[["Z"]]
    output <- list("Z"=Z, "network_embedding"=network_embedding,"node_embedding"=node_embedding)
    return(output)
  }
}

# main function for TWIST
TWIST.fun <- function(sample,K,M=NULL){
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # Tensor decomposition by TWIST
  if(is.null(M)){ # Estimate M
    U_0_list <- InitializationMMSBM(A, ranks=c(K,K,L),estM=TRUE)
    M <- dim(U_0_list[[3]])[2]
  }else{
    U_0_list <- InitializationMMSBM(A, ranks=c(K,K,M))
  }
  
  decomp <- PowerIteration(A, ranks=c(K,K,M), type="TWIST", U_0_list, delta1=1000, delta2=1000, max_iter = 20, tol = 1e-04)
  U <- decomp$node_embedding
  
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
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  return(output)
}

### Method: HOOI higher order orthogonal iteration ###
HOOI.fun <- function(sample,K){
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership

  tucker.re <- tucker(A,ranks=c(K,K,1))
  U <- tucker.re$U[[1]]

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

  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  return(output)
}