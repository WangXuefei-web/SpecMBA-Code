# Numerical study

library(rTensor)
library(ggplot2)
library(reshape2)
library(Matrix)
library(tidyverse)
library(plotly)
library(clue)
library(ggpubr)
library(RColorBrewer)
#library(rMultiNet)
#library(RSpectra)
source('SpecMBA_fun.R')
source('function_GenerateMMSBM.R')
source('compare_fun.R')

### Simulation 1: Properties of SpecMBA alpha=opt,fixed ### 
# (a) varying the number of layers L; 
set.seed(123)
N <- 600 #the number of nodes
M <- 1 #types of networks
L.seq <- seq(3,24,by=3) #number of layers
K <- 3 #number of communities
r <- 0.4 #out-in ratio (q/p, a L-dimensional vector)} default=0.4
rho <- 0.03
d <- rho*N #average degree (a L-dimensional vector)}
Pi <- NULL #balanced
alpha.seq <- seq(-1,5,by=0.2)
simul.max <- 200

simulresult <- data.frame("alpha"=NA,"HamAcc"=NA,"L"=NA)
simulresult <- simulresult[-1,]

for (simul in 1:simul.max){
  print(simul)
  for (L in L.seq){
  sample <- GenerateComplementaryMLSBM(N, M, L, K, Pi, d, r )
  result.seq <- array(NA, dim=c(length(alpha.seq),4)) # alpha L HamAcc loglik
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecMBA <- SpecMBA(sample, alpha, K)
    result.seq[alpha.i,] <- c(alpha,L,re.SpecMBA$HamAcc,re.SpecMBA$loglik) # alpha L HamAcc loglik
    alpha.i <- alpha.i+1
  }
  max.id <- which.max(result.seq[,4])
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('alpha','L','HamAcc','loglik')
  fix.rbind <- result.seq[,1:3] %>% filter(alpha %in% c(-1,0,5))
  re.rbind <- dplyr::bind_rows(fix.rbind,
                               data.frame('alpha'=-100,'L'=result.seq[max.id,2],'HamAcc'=result.seq[max.id,3])) # alpha L HamAcc
  simulresult <- dplyr::bind_rows(simulresult,re.rbind)
  }
}

mean.re <- simulresult %>% group_by(alpha,L) %>% summarise(HamAcc=mean(HamAcc))
mean.re$alpha <- as.character(mean.re$alpha)
mean.re$alpha[mean.re$alpha=='-100'] <- 'opt'
#save(mean.re,file="result/simul1_varyL.Rda")
load("result/simul1_varyL.Rda")

# plot result
#brewer.pal(9,"Set1")
colorbar <- c("#FF7F00","#377EB8","#4DAF4A","#660099")
shapebar <- c(16,17,18,15)
linetypebar <- c(2,4,6,1)
p1.L <- ggplot(data=mean.re,aes(x=L,y=HamAcc,color=alpha,shape=alpha,linetype=alpha))+
  geom_line(linewidth=0.9)+
  geom_point(size=1.5)+
  scale_shape_manual(values = shapebar)+
  scale_color_manual(values = colorbar)+
  scale_linetype_manual(values = linetypebar)+
  theme_gray()+
  theme(legend.position= c(0.85,0.22),legend.key.size = unit(0.7,'cm'))+
  #theme(legend.position = "bottom")+
  xlab('Number of layers')+
  ylab('Clustering accuracy')+
  scale_x_continuous(limits = c(3,24), breaks=L.seq)+
  scale_y_continuous(limits = c(0.3,1), breaks=seq(0.3,1,0.1))

png("result/simul1_varyL.png", height=700,width=600, res=168)
print(p1.L)
dev.off()

# (b) varying the number of nodes N; 
set.seed(123)
N.seq <- seq(300,1800,by=300) #the number of nodes
M <- 1 #types of networks
L <- 6 #number of layers
K <- 3 #number of communities
r <- 0.4 #out-in ratio (q/p, a L-dimensional vector)} default=0.4
rho <- 0.03
Pi <- NULL
alpha.seq <- seq(-1,5,by=0.2)
simul.max <- 200

simulresult <- data.frame("alpha"=NA,"HamAcc"=NA,"N"=NA)
simulresult <- simulresult[-1,]

for (simul in 1:simul.max){
  print(simul)
  for (N in N.seq){
    d <- rho*N
    sample <- GenerateComplementaryMLSBM(N, M, L, K, Pi,  d, r )
    result.seq <- array(NA, dim=c(length(alpha.seq),4)) # alpha N HamAcc loglik
    alpha.i <- 1
    for (alpha in alpha.seq){
      re.SpecMBA <- SpecMBA(sample, alpha, K)
      result.seq[alpha.i,] <- c(alpha,N,re.SpecMBA$HamAcc,re.SpecMBA$loglik) # alpha N HamAcc loglik
      alpha.i <- alpha.i+1
    }
    max.id <- which.max(result.seq[,4])
    result.seq <- as.data.frame(result.seq)
    colnames(result.seq) <- c('alpha','N','HamAcc','loglik')
    fix.rbind <- result.seq[,1:3] %>% filter(alpha %in% c(-1,0,5))
    re.rbind <- dplyr::bind_rows(fix.rbind,
                                 data.frame('alpha'=-100,'N'=result.seq[max.id,2],'HamAcc'=result.seq[max.id,3])) # alpha N HamAcc
    simulresult <- dplyr::bind_rows(simulresult,re.rbind)
  }
}

mean.re <- simulresult %>% group_by(alpha,N) %>% summarise(HamAcc=mean(HamAcc))
mean.re$alpha <- as.character(mean.re$alpha)
mean.re$alpha[mean.re$alpha=='-100'] <- 'opt'
#save(mean.re,file="result/simul1_varyN.Rda")
load("result/simul1_varyN.Rda")

# plot result

p1.N <- ggplot(data=mean.re,aes(x=N,y=HamAcc,color=alpha,shape=alpha,linetype=alpha))+
  geom_line(linewidth=0.9)+
  geom_point(size=1.5)+
  scale_shape_manual(values = shapebar)+
  scale_color_manual(values = colorbar)+
  scale_linetype_manual(values = linetypebar)+
  theme_gray()+
  theme(legend.position= c(0.85,0.22),legend.key.size = unit(0.7,'cm'))+
  #theme(legend.position = "bottom")+
  xlab('Number of nodes')+
  ylab('Clustering accuracy')+
  scale_x_continuous(limits = c(300,1800), breaks=seq(300,1800,300))+
  scale_y_continuous(limits = c(0.3,1), breaks=seq(0.3,1,0.1))

png("result/simul1_varyN.png", height=700,width=600, res=168)
print(p1.N)
dev.off()


# (c) varying the sparsity parameter rho
set.seed(123)
N <- 600 #the number of nodes
M <- 1 #types of networks
L <- 6 #number of layers
K <- 3 #number of communities
r <- 0.4 #out-in ratio (q/p, a L-dimensional vector)} default=0.4
rho.seq <- seq(0.01,0.09,0.02)
Pi <- NULL
alpha.seq <- seq(-1,5,by=0.2)
simul.max <- 200

simulresult <- data.frame("alpha"=NA,"HamAcc"=NA,"rho"=NA)
simulresult <- simulresult[-1,]

for (simul in 1:simul.max){
  print(simul)
  for (rho in rho.seq){
    d <- rho*N #average degree (a L-dimensional vector)}
    sample <- GenerateComplementaryMLSBM(N, M, L, K, Pi,  d, r )
    result.seq <- array(NA, dim=c(length(alpha.seq),4)) # alpha rho HamAcc loglik
    alpha.i <- 1
    for (alpha in alpha.seq){
      re.SpecMBA <- SpecMBA(sample, alpha, K)
      result.seq[alpha.i,] <- c(alpha,rho,re.SpecMBA$HamAcc,re.SpecMBA$loglik) # alpha rho HamAcc loglik
      alpha.i <- alpha.i+1
    }
    max.id <- which.max(result.seq[,4])
    result.seq <- as.data.frame(result.seq)
    colnames(result.seq) <- c('alpha','rho','HamAcc','loglik')
    fix.rbind <- result.seq[,1:3] %>% filter(alpha %in% c(-1,0,5))
    re.rbind <- dplyr::bind_rows(fix.rbind,
                                 data.frame('alpha'=-100,'rho'=result.seq[max.id,2],'HamAcc'=result.seq[max.id,3])) # alpha L HamAcc
    simulresult <- dplyr::bind_rows(simulresult,re.rbind)
  }
}

mean.re <- simulresult %>% group_by(alpha,rho) %>% summarise(HamAcc=mean(HamAcc))
mean.re$alpha <- as.character(mean.re$alpha)
mean.re$alpha[mean.re$alpha=='-100'] <- 'opt'
#save(mean.re,file="result/simul1_varyrho.Rda")
load("result/simul1_varyrho.Rda")

# plot result

p1.rho <- ggplot(data=mean.re,aes(x=rho,y=HamAcc,color=alpha,shape=alpha,linetype=alpha))+
  geom_line(linewidth=0.9)+
  geom_point(size=1.5)+
  scale_shape_manual(values = shapebar)+
  scale_color_manual(values = colorbar)+
  scale_linetype_manual(values = linetypebar)+
  theme_gray()+
  theme(legend.position= c(0.85,0.22),legend.key.size = unit(0.7,'cm'))+
  #theme(legend.position = "bottom")+
  xlab('Sparsity')+
  ylab('Clustering accuracy')+
  scale_x_continuous(limits = c(0.01,0.09), breaks=rho.seq)+
  scale_y_continuous(limits = c(0.3,1), breaks=seq(0.3,1,0.1))
png("result/simul1_varyrho.png", height=700,width=600, res=168)
print(p1.rho)
dev.off()

p1 <- ggarrange(p1.rho,p1.L,p1.N,ncol = 3, nrow = 1,common.legend = T)
png("result/simul1.png", height=500,width=1500, res=180)
print(p1)
dev.off()

# Compared Methods: HOOI, MeanAdj(only mean), OLMF(centered square), SoSBA, SpecM, SpecMBA, TWIST.T, TWIST.F
# Simulation 2: complementary information alpha=opt
colorbar <- c("#00C600","#FF7F00","#377EB8","#4DAF4A","#E41A00","#660099","#EE1289","#984E00")
shapebar <- c(9,16,17,18,6,15,7,13)
linetypebar <- c(3,2,4,6,5,1,2,3)
# (a) varying the number of layers L; 
set.seed(123)
N <- 600 #the number of nodes
M <- 1 #types of networks
L.seq <- seq(3,24,by=3) #number of layers
K <- 3 #number of communities
r <- 0.4 #out-in ratio (q/p, a L-dimensional vector)} default=0.4
rho <- 0.03
d <- rho*N #average degree (a L-dimensional vector)}
Pi <- NULL #balanced

alpha.seq <- seq(-1,5,by=0.2)
simul.max <- 200

simulresult <- data.frame("L"=NA,"SpecMBA"=NA,"MeanAdj"=NA,"OLMF"=NA,"SoSBA"=NA,
                          "HOOI"=NA,"TWIST.T"=NA,"TWIST.F"=NA,"SpecM"=NA)
simulresult <- simulresult[-1,]

for (simul in 1:simul.max){
  print(simul)
  for (L in L.seq){
    sample <- GenerateComplementaryMLSBM(N, M, L, K, Pi, d, r ) # indirected
    SpecMBA.re <- SpecMBA.fun(sample,alpha.seq,K)
    MeanAdj.re <- MeanAdj.fun(sample,K)
    OLMF.re <- OLMF.fun(sample,K)
    SoSBA.re <- SoSBA.fun(sample,K)
    HOOI.re <- HOOI.fun(sample,K)
    TWIST.T.re <- TWIST.fun(sample,K,M=min(4,L))
    TWIST.F.re <- TWIST.fun(sample,K,M=2)
    SpecM.re <- SpecM.fun(sample,alpha.seq,K)
    HamAcc.re <- data.frame('L'=L,'SpecMBA'=SpecMBA.re$HamAcc,'MeanAdj'=MeanAdj.re$HamAcc,
                            'OLMF'=OLMF.re$HamAcc,'SoSBA'=SoSBA.re$HamAcc,'HOOI'=HOOI.re$HamAcc,
                            'TWIST.T'=TWIST.T.re$HamAcc,'TWIST.F'=TWIST.F.re$HamAcc,'SpecM'=SpecM.re$HamAcc)
    simulresult <- dplyr::bind_rows(simulresult,HamAcc.re)
  }
}

mean.re <- simulresult %>% group_by(L) %>% summarise(SpecMBA=mean(SpecMBA),MeanAdj=mean(MeanAdj),
                                                     OLMF=mean(OLMF),SoSBA=mean(SoSBA),HOOI=mean(HOOI),
                                                     TWIST.T=mean(TWIST.T),TWIST.F=mean(TWIST.F),SpecM=mean(SpecM))
#save(mean.re,file="result/simul2_varyL.Rda")
load(file="result/simul2_varyL.Rda")
# plot result
mean.re.long <- mean.re %>% gather(key="Method",value="HamAcc",SpecMBA:SpecM)


p2.L <- ggplot(data=mean.re.long,aes(x=L,y=HamAcc,color=Method,shape=Method,linetype=Method))+
  geom_line(linewidth=0.9)+
  geom_point(size=1.5)+
  scale_shape_manual(values = shapebar)+
  scale_color_manual(values = colorbar)+
  scale_linetype_manual(values = linetypebar)+
  theme_gray()+
  #theme(legend.position = "bottom")+
  theme(legend.position= c(0.85,0.3),legend.key.size = unit(0.7,'cm'))+
  xlab('Number of layers')+
  ylab('Clustering accuracy')+
  scale_x_continuous(limits = c(3,24), breaks=L.seq)+
  scale_y_continuous(limits = c(0.3,1), breaks=seq(0.3,1,0.1))+
  ggtitle("Doubly-informative ")+
  theme(plot.title = element_text(hjust = 0.5,size=12))

png("result/simul2_varyL.png", height=700,width=600, res=168)
print(p2.L)
dev.off()

# (b) varying the sparsity parameter rho
set.seed(123)
N <- 600 #the number of nodes
M <- 1 #types of networks
L <- 6 #number of layers
K <- 3 #number of communities
r <- 0.4 #out-in ratio (q/p, a L-dimensional vector)} default=0.4
rho.seq <- seq(0.01,0.09,0.02)
Pi <- NULL #balanced
alpha.seq <- seq(-1,5,by=0.2)
simul.max <- 200

simulresult <- data.frame("rho"=NA,"SpecMBA"=NA,"MeanAdj"=NA,"OLMF"=NA,"SoSBA"=NA,
                          "HOOI"=NA,"TWIST.T"=NA,"TWIST.F"=NA,"SpecM"=NA)
simulresult <- simulresult[-1,]

for (simul in 1:simul.max){
  print(simul)
  for (rho in rho.seq){
    d <- rho*N #average degree (a L-dimensional vector)}
    sample <- GenerateComplementaryMLSBM(N, M, L, K, Pi,  d, r ) # indirected
    SpecMBA.re <- SpecMBA.fun(sample,alpha.seq,K)
    MeanAdj.re <- MeanAdj.fun(sample,K)
    OLMF.re <- OLMF.fun(sample,K)
    SoSBA.re <- SoSBA.fun(sample,K)
    HOOI.re <- HOOI.fun(sample,K)
    TWIST.T.re <- TWIST.fun(sample,K,M=3)
    TWIST.F.re <- TWIST.fun(sample,K,M=2)
    SpecM.re <- SpecM.fun(sample,alpha.seq,K)
    HamAcc.re <- data.frame('rho'=rho,'SpecMBA'=SpecMBA.re$HamAcc,'MeanAdj'=MeanAdj.re$HamAcc,
                            'OLMF'=OLMF.re$HamAcc,'SoSBA'=SoSBA.re$HamAcc,'HOOI'=HOOI.re$HamAcc,
                            'TWIST.T'=TWIST.T.re$HamAcc,'TWIST.F'=TWIST.F.re$HamAcc,'SpecM'=SpecM.re$HamAcc)
    simulresult <- dplyr::bind_rows(simulresult,HamAcc.re)
  }
}

mean.re <- simulresult %>% group_by(rho) %>% summarise(SpecMBA=mean(SpecMBA),MeanAdj=mean(MeanAdj),
                                                       OLMF=mean(OLMF),SoSBA=mean(SoSBA),HOOI=mean(HOOI),
                                                       TWIST.T=mean(TWIST.T),TWIST.F=mean(TWIST.F),SpecM=mean(SpecM))
#save(mean.re,file="result/simul2_varyrho.Rda")
load(file="result/simul2_varyrho.Rda")
# plot result
mean.re.long <- mean.re %>% gather(key="Method",value="HamAcc",SpecMBA:SpecM)

p2.rho <- ggplot(data=mean.re.long,aes(x=rho,y=HamAcc,color=Method,shape=Method,linetype=Method))+
  geom_line(linewidth=0.9)+
  geom_point(size=1.5)+
  scale_shape_manual(values = shapebar)+
  scale_color_manual(values = colorbar)+
  scale_linetype_manual(values = linetypebar)+
  theme_gray()+
  theme(legend.position= c(0.85,0.3),legend.key.size = unit(0.7,'cm'))+
  xlab('Sparsity')+
  ylab('Clustering accuracy')+
  scale_x_continuous(limits = c(0.01,0.09), breaks=rho.seq)+
  scale_y_continuous(limits = c(0.3,1), breaks=seq(0.3,1,0.1))+
  ggtitle("Doubly-informative ")+
  theme(plot.title = element_text(hjust = 0.5,size=12))

png("result/simul2_varyrho.png", height=700,width=600, res=168)
print(p2.rho)
dev.off()

# (c) varying the number of communities K; 
set.seed(123)
 #the number of nodes
M <- 1 #types of networks
L <- 6 #number of layers
K.seq <- seq(2,6,by=1) #number of communities
r <- 0.4 #out-in ratio (q/p, a L-dimensional vector)} default=0.4
rho <- 0.03
Pi <- NULL #balanced
alpha.seq <- seq(-1,5,by=0.2)
simul.max <- 200

simulresult <- data.frame("K"=NA,"SpecMBA"=NA,"MeanAdj"=NA,"OLMF"=NA,"SoSBA"=NA,
                          "HOOI"=NA,"TWIST.T"=NA,"TWIST.F"=NA,"SpecM"=NA)
simulresult <- simulresult[-1,]

for (simul in 1:simul.max){
  print(simul)
  for (K in K.seq){
    N <- 200*K #the number of nodes
    d <- rho*N #average degree (a L-dimensional vector)}
    sample <- GenerateMMSBM(N, M, L, K, Pi, d, r ) # indirected
    SpecMBA.re <- SpecMBA.fun(sample,alpha.seq,K)
    MeanAdj.re <- MeanAdj.fun(sample,K)
    OLMF.re <- OLMF.fun(sample,K)
    SoSBA.re <- SoSBA.fun(sample,K)
    HOOI.re <- HOOI.fun(sample,K)
    TWIST.T.re <- TWIST.fun(sample,K,M=1)
    SpecM.re <- SpecM.fun(sample,alpha.seq,K)
    HamAcc.re <- data.frame('K'=K,'SpecMBA'=SpecMBA.re$HamAcc,'MeanAdj'=MeanAdj.re$HamAcc,
                            'OLMF'=OLMF.re$HamAcc,'SoSBA'=SoSBA.re$HamAcc,'HOOI'=HOOI.re$HamAcc,
                            'TWIST.T'=TWIST.T.re$HamAcc,'SpecM'=SpecM.re$HamAcc)
    simulresult <- dplyr::bind_rows(simulresult,HamAcc.re)
  }
}

###
mean.re <- simulresult %>% group_by(K) %>% summarise(SpecMBA=mean(SpecMBA),MeanAdj=mean(MeanAdj),
                                                     OLMF=mean(OLMF),SoSBA=mean(SoSBA),HOOI=mean(HOOI),
                                                     TWIST.T=mean(TWIST.T),SpecM=mean(SpecM))
#save(mean.re,file="result/simul2_varyK.Rda")
load(file="result/simul2_varyK.Rda")
# plot result
mean.re.long <- mean.re %>% gather(key="Method",value="HamAcc",SpecMBA:SpecM)
mean.re.long <- add_row(mean.re.long, K = NULL, Method="TWIST.F", HamAcc=NULL)

p2.K <- ggplot(data=mean.re.long,aes(x=K,y=HamAcc,color=Method,shape=Method,linetype=Method))+
  geom_line(linewidth=0.9)+
  geom_point(size=1.5)+
  scale_shape_manual(values = shapebar)+
  scale_color_manual(values = colorbar)+
  scale_linetype_manual(values = linetypebar)+
  theme_gray()+
  #theme(legend.position = "bottom")+
  theme(legend.position= c(0.12,0.3),legend.key.size = unit(0.7,'cm'))+
  xlab('Number of communities')+
  ylab('Clustering accuracy')+
  scale_x_continuous(limits = c(2,6), breaks=K.seq)+
  scale_y_continuous(limits = c(0.3,1), breaks=seq(0.3,1,0.1))+
  ggtitle("Covariance-uninformative")+
  theme(plot.title = element_text(hjust = 0.5,size=12))

png("result/simul2_varyK.png", height=700,width=600, res=168)
print(p2.K)
dev.off()

# (d) mean-uninformative varying the sparsity parameter rho
set.seed(123)
N <- 600 #the number of nodes
M <- 1 #types of networks
L <- 6 #number of layers
K <- 3 #number of communities
r <- 0.4 #out-in ratio (q/p, a L-dimensional vector)} default=0.4
rho.seq <- seq(0.01,0.09,0.02)

alpha.seq <- seq(-1,5,by=0.2)
simul.max <- 200

simulresult <- data.frame("rho"=NA,"SpecMBA"=NA,"MeanAdj"=NA,"OLMF"=NA,"SoSBA"=NA,
                          "HOOI"=NA,"TWIST.T"=NA,"TWIST.F"=NA,"SpecM"=NA)
simulresult <- simulresult[-1,]

for (simul in 1:simul.max){
  print(simul)
  for (rho in rho.seq){
    d <- rho*N #average degree (a L-dimensional vector)}
    sample <- GenerateComplementaryMLSBM(N, M, L, K, Pi,d, r, NetType = 3 ) # indirected
    SpecMBA.re <- SpecMBA.fun(sample,alpha.seq,K)
    #MeanAdj.re <- MeanAdj.fun(sample,K)
    #OLMF.re <- OLMF.fun(sample,K)
    SoSBA.re <- SoSBA.fun(sample,K)
    #HOOI.re <- HOOI.fun(sample,K)
    #TWIST.T.re <- TWIST.fun(sample,K,M=6)
    #TWIST.F.re <- TWIST.fun(sample,K,M=3)
    SpecM.re <- SpecM.fun(sample,alpha.seq,K)
    HamAcc.re <- data.frame('rho'=rho,'SpecMBA'=SpecMBA.re$HamAcc,'MeanAdj'=MeanAdj.re$HamAcc,
                            'OLMF'=OLMF.re$HamAcc,'SoSBA'=SoSBA.re$HamAcc,'HOOI'=HOOI.re$HamAcc,
                            'TWIST.T'=TWIST.T.re$HamAcc,'TWIST.F'=TWIST.F.re$HamAcc,'SpecM'=SpecM.re$HamAcc)
    simulresult <- dplyr::bind_rows(simulresult,HamAcc.re)
  }
}

mean.re <- simulresult %>% group_by(rho) %>% summarise(SpecMBA=mean(SpecMBA),MeanAdj=mean(MeanAdj),
                                                       OLMF=mean(OLMF),SoSBA=mean(SoSBA),HOOI=mean(HOOI),
                                                       TWIST.T=mean(TWIST.T),TWIST.F=mean(TWIST.F),SpecM=mean(SpecM))
#save(mean.re,file="result/simul2_2.Rda")
load(file="result/simul2_2.Rda")
# plot result
mean.re.long <- mean.re %>% gather(key="Method",value="HamAcc",SpecMBA:SpecM)


p2.2 <- ggplot(data=mean.re.long,aes(x=rho,y=HamAcc,color=Method,shape=Method,linetype=Method))+
  geom_line(linewidth=0.9)+
  geom_point(size=1.5)+
  scale_shape_manual(values = shapebar)+
  scale_color_manual(values = colorbar)+
  scale_linetype_manual(values = linetypebar)+
  theme_gray()+
  theme(legend.position= c(0.15,0.7),legend.key.size = unit(0.7,'cm'))+
  xlab('Sparsity')+
  ylab('Clustering accuracy')+
  scale_x_continuous(limits = c(0.01,0.09), breaks=rho.seq)+
  scale_y_continuous(limits = c(0.3,1), breaks=seq(0.3,1,0.1))+
  ggtitle("Mean-uninformative")+
  theme(plot.title = element_text(hjust = 0.5,size=12))

png("result/simul2_2.png", height=700,width=500, res=168)
print(p2.2)
dev.off()

# Plot
p2 <- ggarrange(p2.K,p2.rho,p2.L,p2.2,ncol = 2, nrow = 2,common.legend = T)#
png("result/simul2.png", height=1200,width=1500, res=200)
print(p2)
dev.off()
