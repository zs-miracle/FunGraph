#' @title Get_effect
#' @description The net genetic effect of each gene locus (SNP) was calculated
#' @importFrom stats coef cor cov kmeans lm optim sd time
#'@param funmap_data_snp  Genotypic data.only accept numeric value for calculation, SNP genotype data should be converted into 0,1 matrix.
#' @param funmap_data_snp_par  The H1 parameters of each gene locus were calculated by FunMap function.code：
#' @param times A series of times for measuring phenotypic data.
#' @return The net genetic effect of each gene locus (SNP).Type is data frame.
#' @export

Get_effect<-function(funmap_data_snp,funmap_data_snp_par,times){
  diff_vg <- matrix(NA,dim(funmap_data_snp)[1],length(times))
  colnames(diff_vg) <- c(times)
  rownames(diff_vg) <- rownames(funmap_data_snp)
  get_miu<-function(par,times,options=list())
  {
    y <- par[1]/(1+par[2]*exp(-par[3]*times))
    return (y)
  }
  for (i in 1:dim(funmap_data_snp)[1]) {
    x0 <- as.numeric(which(funmap_data_snp[i,]==0))
    x1 <- as.numeric(which(funmap_data_snp[i,]==1))
    
    x0_num <- length(x0)
    x1_num <- length(x1)
    
    p0 <- (x0_num*2)/((x1_num+x0_num)*2)
    p1 <- (x1_num*2)/((x1_num+x0_num)*2)
    
    
    mean_x0 <- get_miu(as.numeric(funmap_data_snp_par[i,1:3]),times)
    mean_x1 <- get_miu(as.numeric(funmap_data_snp_par[i,4:6]),times)
    AE <- (mean_x1-mean_x0)/2
    
    Vg <- 2*p1*p0*(AE^2)
    Vg <- sqrt(Vg)
    diff_vg[i,] <- Vg
  }
  return(as.data.frame(diff_vg))
}
