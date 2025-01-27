# Part of this code document is included in ``rMultiNet" provided by Li et al. (2023)
#' @name       GenerateMMSBM
#' @title      A function for generating a tensor
#' @description a method to generate a tensor for multilayer network
###Input:
#' @param      n the number of nodes
#' @param      m types of networks
#' @param      K number of communities
#' @param      L number of layers
#' @param      d average degree (a L-dimensional vector)
#' @param      r out-in ratio (q/p, a L-dimensional vector)
#' @import    rTensor
#' @import     glmnet
#' @import    lpSolve
#' @import     matlabr
#' @import    truncnorm
#' @import    Matrix
#' @export


generate_sbm<-function(n,K,r,d,membership)
{
  p=0.6*d/n
  q=p*r
  B=matrix(q,K,K)
  diag(B)=p

  Z=matrix(0,nrow=n,ncol=K)
  for(i in 1:n)
  {
    Z[i,membership[i]]=1
  }

  P=Z%*%B%*%t(Z)

  A=matrix(NA,ncol=n,nrow=n)
  diag(A)=0
  for (i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      A[i,j]=rbinom(1,1,P[i,j])
      A[j,i]=A[i,j]
    }
  }
  return(A)
}


generate_tensor<-function(n,m,K,L,r,d,Pi,Pi_network_type)
{
  ###Generate Membership matrix n by m
  local_membership=matrix(NA,nrow=n,ncol=m)
  colnames(local_membership)=c(1:m)
  for(i in 1:m)
  {
    local_membership[,i]=sample(x = K, size = n, replace = TRUE, prob = Pi)
  }

  if(is.null(Pi_network_type))
  {
    Pi_network_type=rep(1/m,m);
  }

  ###Assign network type for each layer
  network_type=sample(x = m, size = L, replace = TRUE, prob = Pi_network_type)

  mylist=list()

  for (i in 1:L)
  {
    mylist[[i]]=generate_sbm(n=n,K=K,r=r[i],d=d[i],membership=local_membership[,network_type[i]])
  }

  global_membership=rep(0,n)
  for(i in 1:m)
  {
    global_membership=global_membership+K^(i-1)*(local_membership[,i]-1)
  }
  indices <- c(n,n,L)
  arr<-array(as.numeric(unlist(mylist)), dim=indices)
  arrT <- as.tensor(arr)
  output<-list(arrT,global_membership,local_membership,network_type)
  return(output)
}


GenerateMMSBM <- function(n,m,L,K,Pi=NULL,d=NULL,r=NULL)
{

  ###Initialize the average degree list
  if (is.null(d))
  {
    d_list = abs(rnorm(L,5,5))
  }else if (length(d)==L){
    d_list = d
  }else{
    d_list = rep(d,L)
  }
   

  ###Initialize the out-in ratio of each layer list
  if (is.null(r))
  {
    r_list = rep(0.4,L)
  }else if (length(r)==L)
  {
    r_list = r
  }else
  {
    r_list = rep(r,L)
  }
 
  ### Initialize the community size
  Pi_network_type = rep(1/m,m)
  if (is.null(Pi)){
    Pi_network_community = rep(1/K,K) # balanced community
  }else{
    Pi_network_community = Pi/sum(Pi)
  }
  #library(rTensor)
  temp = generate_tensor(n=n,m=m,K=K,L=L,r=r_list,d=d_list,Pi=Pi_network_community,Pi_network_type=Pi_network_type)
  arrT = temp[[1]]
  Global_node_type_power = temp[[2]]
  number_global_community = length(unique(Global_node_type_power))
  list <- list('adjT' = arrT, 'global_membership' = Global_node_type_power,
               'number_global_community' = number_global_community)
  return(list)

}


### For special simulation: GenerateComplementaryMLSBM ###
# B is given in the following function
generate_sbm_givenB<-function(n,K,r,d,membership,typeB)
{
  # particularly designed for K = 3
  p=0.6*d/n
  q=r*p
  B=matrix(q,K,K)
  if (K>=3){
    if (typeB == 0){
      diag(B)=p
      B[1,2]=p
      B[2,1]=p
    }else if (typeB == 1){
      diag(B)=p
      B[2,3]=p
      B[3,2]=p
    }else if (typeB == 2){
      B[1,1]=p
      B[1,3]=p
      B[3,1]=p
      B[3,3]=p
    }else if (typeB == 3){
      B[1,1]=p
    }else if (typeB == 4){ #disassortative
      B[1:3,2]=p
      B[2,1:3]=p
    }else if (typeB == 5){ #disassortative
      B[1,3]=p
      B[3,1]=p
      B[2,3]=p#
      B[3,2]=p#
    }else if (typeB == 6){ #disassortative
      B[1,3]=p#
      B[3,1]=p#
      B[1,2]=p#
      B[2,1]=p#
      #B[upper.tri(B,diag = FALSE)]=p
      #B[lower.tri(B,diag = FALSE)]=p
    }else{
      diag(B)=p
    }
    diag(B)[-1:-3]=p
  }else{
    diag(B)=p
  }
    
  
  Z=matrix(0,nrow=n,ncol=K)
  for(i in 1:n)
  {
    Z[i,membership[i]]=1
  }
  
  P=Z%*%B%*%t(Z)
  
  A=matrix(NA,ncol=n,nrow=n)
  diag(A)=0
  for (i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      A[i,j]=rbinom(1,1,P[i,j])
      A[j,i]=A[i,j]
    }
  }
  return(A)
}


generate_tensor_givenB<-function(n,m,K,L,r,d,Pi,Pi_network_type,disass)
{
  ###Generate Membership matrix n by m
  local_membership=matrix(NA,nrow=n,ncol=m)
  colnames(local_membership)=c(1:m)
  for(i in 1:m)
  {
    local_membership[,i]=sample(x = K, size = n, replace = TRUE, prob = Pi)
  }
  
  if(is.null(Pi_network_type))
  {
    Pi_network_type=rep(1/M,M);
  }
  
  ###Assign network type for each layer
  network_type=sample(x = m, size = L, replace = TRUE, prob = Pi_network_type)
  # B_type=c(0,1,2,3) # designed for assortative complementary layers
  # B_type=c(0,1,2,4,5,6) # designed for mixed complementary layers
  disass.B=c(0,1,2,4,5,6)
  mylist=list()
  for (i in 1:L)
  {
    if (!disass)
    mylist[[i]]=generate_sbm_givenB(n=n,K=K,r=r[i],d=d[i],
                                    membership=local_membership[,network_type[i]],
                                    typeB=(i-1)%%4)
    else{ # three more disassortative layers
    mylist[[i]]=generate_sbm_givenB(n=n,K=K,r=r[i],d=d[i],
                                      membership=local_membership[,network_type[i]],
                                      typeB=disass.B[((i-1)%%6)+1]) 
    }
  }
  
  global_membership=rep(0,n)
  for(i in 1:m)
  {
    global_membership=global_membership+K^(i-1)*(local_membership[,i]-1)
  }
  indices <- c(n,n,L)
  arr<-array(as.numeric(unlist(mylist)), dim=indices)
  arrT <- as.tensor(arr)
  output<-list(arrT,global_membership,local_membership,network_type)
  return(output)
}


GenerateComplementaryMLSBM <- function(n,m,L,K,Pi=NULL,d=NULL,r=NULL,disass=FALSE)
{
  
  ### Initialize the average degree list
  if (is.null(d))
  {
    d_list = abs(rnorm(L,5,5))
  }else if (length(d)==L){
    d_list = d
  }
  else{
    d_list = rep(d,L)
  }
  
  ### Initialize the out-in ratio of each layer list
  if (is.null(r))
  {
    r_list = rep(0.4,L)
  }else if (length(r)==L)
  {
    r_list = r
  }
  else
  {
    r_list = rep(r,L)
  }
  
  ### Initialize the community size
  Pi_network_type = rep(1/m,m)
  if (is.null(Pi)){
    Pi_network_community = rep(1/K,K) # balabced community
  }
  else{
    Pi_network_community = Pi/sum(Pi)
  }
  
  #library(rTensor)
  temp = generate_tensor_givenB(n=n,m=m,K=K,L=L,r=r_list,d=d_list,Pi=Pi_network_community,Pi_network_type=Pi_network_type,disass)
  arrT = temp[[1]]
  Global_node_type_power = temp[[2]]
  number_global_community = length(unique(Global_node_type_power))
  list <- list('adjT' = arrT, 'global_membership' = Global_node_type_power,
               'number_global_community' = number_global_community)
  return(list)
  
}
