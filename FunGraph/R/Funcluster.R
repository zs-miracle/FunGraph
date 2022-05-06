#' @title Funcluster 
#' @description According to modularity theory,functional clustering (FunClu) was 
#' introduced to cluster genetic effects into different functional modules. 
#' A hybrid of the EM and simplex algorithms were implanted to obtain the functional modules.
#' @importFrom stats coef cor cov kmeans lm optim sd time
#' @param data The net genetic effects of all loci were calculated by Get_effect function(remarker_effect).
#' @param times rep A series of times for measuring phenotypic data
#' @param rep Sets the number of funclusters specified, usually 1
#' @param k Specify the number of clusters.
#' @param legendre_order Specifies the order of Legendre polynomials
#' @return Type is list.The list contains BIC Maximum、likelihood、Legendre polynomials parameters、Mean curve parameters for each cluster.
#' @export

Funcluster <- function(data,times,rep,k,legendre_order){
  get_init_par <- function(data,k,legendre_order){
    get_legendre_par <- function(y,legendre_order,x){
      #lm_method
      get_legendre_matrix <- function(x,legendre_order){
        legendre_coef <- legendre.polynomials(n=legendre_order, normalized=F)
        legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
          polynomials=legendre_coef,x=scaleX(x, u=-1, v=1))))
        colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
        return(legendre_matrix[,2:(legendre_order+1)])
      }
      legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
      return(legendre_par)
    }
    #基于k-means获得初始参数--------------------------------------
    init_cluster <- kmeans(data,centers = k,iter.max = 1000)
    pro <- table(init_cluster$cluster)/nrow(data) #各子模型概率概率
    
    cuM <- init_cluster$centers
    cusd <- diag(cov(data))
    
    #init_curve_para <-  cuM #输入data为df_par
    #init_curve_para <- t(pbsapply(1:nrow(cuM),function(c)coef(lm(as.numeric(cuM[c,])~get_legendre_matrix(age,legendre_order)))))
    init_curve_para <- t(sapply(1:k,function(c)get_legendre_par(as.numeric(cuM[c,]),legendre_order,times)))
    init_sd_para <- c(mean(cusd),0.5) #SAD1
    init_pro <- pro
    
    return_object <- list(init_sd_para,init_curve_para,init_pro)
    names(return_object)<-c("init_sd_par","init_curve","init_pro")
    return(return_object)
  }
  get_cluster <- function(data,k,input,legendre_order){
    Delta <- 100; iter <- 0; itermax <- 100;
    SAD1_get_matrix <- function(par,data){
      p <-  ncol(data)
      v2 <- par[1]
      phi <- par[2]
      tmp <- (1-phi^2)
      sigma <- array(dim=c(p,p))
      for(i in 1:p){
        sigma[i,i:p] <- phi^( c(i:p) - i ) * (1-phi^(2*i ))/tmp
        sigma[i:p,i] <- sigma[i,i:p]}
      sigma <- sigma*abs(v2)
      return(sigma)
    }
    
    legendre_fit <- function(par,x_interpolation){
      legendre_formula <-  sapply(1:length(par),function(c)
        par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
      legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
        polynomials=legendre_formula,x=scaleX(x_interpolation, u=-1, v=1))))
      colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
      y_interpolation <- rowSums(legendre_matrix)
      return(y_interpolation)
    }
    
    mle <- function(par,data,prob){
      par1 <- par[1:2]
      par2 <- matrix(par[-c(1:2)],nrow = k,ncol = (legendre_order+1))
      temp_S <- sapply(1:k, function(c) dmvnorm(data,
                                                legendre_fit(par2[c,],times),
                                                SAD1_get_matrix(par1,data))*prob[c] )
      LL <- -sum(log(rowSums(temp_S)))
      return(LL)
    }
    cat(paste0("Start biFunClu Calculation ","\n","Cluster_number=",k," Legendre_order=", legendre_order))
    while ( Delta > 1 && iter <= itermax ) {
      # initiation
      if(iter == 0){
        init_sd_para <- input[[1]]
        init_curve_para <- input[[2]]
        pro <- input[[3]]
      }
      
      #E step, calculate the posterior probability
      old_par <- c(init_sd_para,init_curve_para)
      LL_mem <- mle(old_par,data,pro)
      mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                               legendre_fit(init_curve_para[c,],times),
                                               SAD1_get_matrix(init_sd_para,data))*pro[c] )
      omega <- mvn.c/rowSums(mvn.c)
      #M step, calculate parameters
      pro <- colSums(omega)/sum(omega)
      new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"
                           ,control = list(maxit=5000,trace=T)
      ))
      if ('try-error' %in% class(new_par))
        break
      L_Value <- new_par$value
      init_sd_para <- new_par$par[1:2]
      init_curve_para <- matrix(new_par$par[-c(1:2)],nrow = k)
      Delta <- abs(L_Value-LL_mem)
      if (Delta > 20000) break
      cat('\n',"iter=",iter,"LL=",L_Value,'\n')
      iter <- iter+1; LL_mem <- L_Value
    }
    
    BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
    #plot-----------
    cluster <- apply(omega,1,which.max)
    #clustered_df <- data.frame(cbind(row.names(data),df_fitted,cluster))
    clustered_df <- data.frame(cbind(row.names(data),data,cluster))
    colnames(clustered_df) <- c("row.names(data)",times,"cluster")
    long_df <- melt(clustered_df,id.vars=c("row.names(data)","cluster"))
    #long_df[,4] <- as.numeric(long_df[,4])
    #long_df[,3] <- as.numeric(long_df[,3])
    colnames(long_df) <- c("gene","cluster","time","fpkm")
    p <- ggplot()+geom_line(long_df,mapping=aes(as.numeric(as.character(time)),fpkm,group=gene,
                                                colour= as.character(cluster)))+
      facet_wrap(cluster,scales = "fixed")+
      theme(legend.position="none") + xlab("Total_fpkm")+ylab("individual_fpkm")
    
    clustered_df <- clustered_df[,-1]
    return_object <- list(init_sd_para,init_curve_para,pro,LL_mem,BIC,clustered_df,p)
    names(return_object)<-c("sd_par", "curve_par", "pro", "LL", "BIC", "clustered_data","plot")
    
    return(return_object)
  }
  get_BIC <- function(rep,k,legendre_order){
    input <- pblapply(1:rep, function(c) get_init_par(data,k,legendre_order))
    core.number <- detectCores()
    cl <- makeCluster(getOption("cl.cores", core.number))
    clusterEvalQ(cl, {library(mvtnorm)})
    clusterEvalQ(cl, {library(reshape2)})
    clusterEvalQ(cl, {library(pbapply)})
    clusterEvalQ(cl, {library(orthopolynom)})
    clusterEvalQ(cl, {library(ggplot2)})
    clusterExport(cl, c("get_init_par","get_cluster","input","k","rep","legendre_order",'times')
                  ,envir=environment()
    )
    result <- pblapply(1:rep,function(c)get_cluster(data,k,input[[c]],legendre_order))
    stopCluster(cl)
    
    out_df <- matrix(NA,nrow = 3,ncol = rep)
    out_df[1,] <- sapply(1:rep,function(c)length(table(result[[c]]$clustered_data[,(dim(data)[2]+1)])))
    out_df[2,] <- sapply(1:rep,function(c)result[[c]]$BIC)
    out_df[3,] <- sapply(1:rep,function(c)result[[c]]$LL)
    rownames(out_df) <- c("cluster_number","BIC","LL")
    colnames(out_df) <- paste0("rep",1:rep)
    return(list(t(out_df),result))
  }
  get_BIC(rep,k,legendre_order)
}
