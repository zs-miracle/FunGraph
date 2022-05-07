#' @title FunMap 
#' @description FunMap is crucial to this model, for it excavates how specific QTLs 
#' determines the complex trait. The mean vector and covariance 
#' structure should be modeled according to the design of the experiment.
#' @import  mvtnorm parallel pbapply ggplot2 deSolve glmnet orthopolynom writexl
#' @importFrom  stats coef cor cov kmeans lm optim sd time
#' @importFrom  mvtnorm dmvnorm
#' @param funmap_data The phenotypic data.Remove or replace missing values in the data; Log-transformed the phenotypic data by command log() if the number in the data vary widely is recommended, such as numbers differ by more than five order of magnitude.
#' @param funmap_data_snp  Genotypic data.only accept numeric value for calculation, SNP genotype data should be converted into 0,1 matrix.
#' @param times A series of times for measuring phenotypic data
#' @param rownumber The default value is the total amount of genotype data.For simple tests, enter any total value between 1- total value of genotype data.
#' @return Type is list.The list contains H0 result、H1 result、All H1 value、LR of each SNP.
#' @export
Funmap_function <- function(funmap_data,funmap_data_snp,times,rownumber=dim(funmap_data_snp)[1])
{
  options (warn = -1)
  #获得均值向量
  get_miu<-function(par,times,options=list())
  {
    y <- par[1]/(1+par[2]*exp(-par[3]*times))
    return (y)
  }
  #SAD协方差矩阵
  SAD1_get_matrix = function(par, times = t, options=list()) {
    n <- ifelse (is.vector(times), length(times), NCOL(times) )
    phi<- par[1]
    v2 <- par[2]
    tmp <- (1-phi^2)
    sigma <- array(1, dim=c(n,n))
    for(i in 1:n)
    {
      sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i))/tmp
      sigma[i:n,i] <- sigma[i,i:n]
    }
    sigma <- sigma * abs(v2)
    return(sigma);
  }
  #获得SAD协方差矩阵初始参数
  SAD_init_par_fn<-function(dat)
  {
    m <- dim(dat)[2];
    dat.t0<-dat[ , 1] ;
    par.s2 <- sd(dat.t0, na.rm=T);
    
    dat.t1<-dat[ , 2 ];
    par.rho <- cor(dat.t1, dat.t0, use="complete.obs");
    par.phi <- par.rho/par.s2^2;
    return(c( par.phi, par.s2));
  }
  #获得均值向量初始参数
  LC_init_par_fn<-function(pheno_data, times){
    m <- dim(pheno_data)[2]
    um <- mean(pheno_data[,m], na.rm=T)
    um.1 <- mean(pheno_data[,m-1], na.rm=T)
    u2 <- mean(pheno_data[,2], na.rm=T)
    u1 <- mean(pheno_data[,1], na.rm=T)
    if (u2>u1)
      u0 <- u1*3/4
    else
      u0 <- u1*4/3
    par.a <- um * um / um.1
    par.b <- par.a / u0 -1
    if(par.b==Inf){
      par.b <- u2
    }
    par.r <- try( -log((par.a/um-1)/par.b)/m )
    if (is.na(par.r)|par.r==Inf|par.r==-Inf) par.r <- 0.6
    if(class(par.r)=="try-error")
      par.r <- 1
    return(c(par.a, par.b, par.r))
  }
  times=times
  SAD_init_par<-SAD_init_par_fn(funmap_data)
  LC_init_par<-LC_init_par_fn(funmap_data,times)
  get_funmap_H0 = function(data,times,par){
    miu = get_miu(par,times)
    sigma = SAD1_get_matrix(par = c(par[4],par[5]), times = times )
    L0 = c()
    L0 = sum(dmvnorm(data,miu,sigma,log = T))
    return(-L0)
  }
  H0_result<-optim(par=c(LC_init_par,0.7,2),get_funmap_H0,data=funmap_data,times=times,method="BFGS",control=list(maxit=200000))
  H0_result<-optim(par=c(H0_result$par),get_funmap_H0,data=funmap_data,times=times,method="BFGS",control=list(maxit=20000))
  
  H1_function <- function(funmap_data_snp_row){
    x0<-which(funmap_data_snp_row==0)
    x1<-which(funmap_data_snp_row==1)
    get_funmap_H1 <-  function(par,times){
      miu0 = get_miu(par[1:3],times)
      miu1 = get_miu(par[4:6],times)
      sigma = SAD1_get_matrix(par = c(par[7:8]),times = times )
      L0 = sum(dmvnorm(funmap_data[x0,],miu0,sigma,log = T))
      L1 = sum(dmvnorm(funmap_data[x1,],miu1,sigma,log = T))
      L2= L1+L0
      return(-L2)
    }
    H1_result<-optim(par=c(H0_result$par[1:3],H0_result$par[1:3],H0_result$par[4:5]),
                     times=times,get_funmap_H1,method="BFGS",control=list(maxit=20000))
    H1_value <- H1_result$value
    H1_par <- H1_result$par
    result <- list(H1_par=H1_par,H1_value=H1_value)
    return(result)
  }
  core.number <- detectCores()
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterExport(cl, c("dmvnorm"),envir=environment())
  All_H1_result <- pbapply(funmap_data_snp[1:rownumber,],1,H1_function,cl=cl)
  stopCluster(cl)
  All_H1_value <-sapply(1:dim(funmap_data_snp[1:rownumber,])[1],function(c)(All_H1_result[[c]][2]))
  LR<-(2*(H0_result$value-as.numeric(All_H1_value)))
  Fin_result <- list(H0_result=H0_result,All_H1_result=All_H1_result,
                     All_H1_value=All_H1_value,LR=LR)
  return(Fin_result)
}
