#' @title Decomposition_curve
#' @importFrom stats coef cor cov kmeans lm optim sd time
#' @param exp_index vector of the time point
#' @param legendre_order scalar of LOP order
#' @param data Net effect value of the mean curve(matrix of observed data)
#' @param alpha Adjust the degree of punishment. The default value is 0. If necessary, the reasonable range is between 0 and 1
#' @return Type is list.The list contains LOP parameter、interactions between sites、the elements that make up a network.
#' @export


Decomposition_curve<- function(exp_index,legendre_order,data,alpha=0){
  get_legendre_par <- function(times,order,data,alpha) {
    set.seed(1)
    get_interaction <- function(data,col){
      n <- nrow(data)
      clean_data <- data
      gene_list <- list()
      m <- clean_data[,col]
      M <- clean_data[,-col]
      x_matrix <- M
      x_matrix <- as.matrix(x_matrix)
      #vec <- sapply(1:length(M[1,]),function(c)cor(m,M[,c]))
      #x_matrix <- M[,which( vec %in% -sort(-vec)[1:(n/log(n))] )]
      #x_matrix <- as.matrix(x_matrix)
      name <- colnames(clean_data)
      ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",
                             family="gaussian",nfolds = 10,alpha = 0)
      best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
      
      fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", family="gaussian",
                           nfolds = 10,alpha = 1,
                           penalty.factor = 1/best_ridge_coef,
                           keep = TRUE)
      best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
      
      gene_list_one <- list()
      gene_list_one[[1]] <- name[col]
      gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
      gene_list_one[[3]] <- best_alasso_coef1@x[-1]
      gene_list[[col]] <- gene_list_one
      
      return(gene_list_one)
    }
    cluster_mean <- data
    rownames(cluster_mean) <- rownames(data)
    module_relationship <- pblapply(1:nrow(data),function(c) get_interaction(t(cluster_mean),c))
    get_effect <- function(pars,effect,times,order){
      if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
      LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
      d_LOP_fit <-  sapply(1:length(pars),function(c)
        pars[c]*polynomial.derivatives(LOP)[[c+1]])
      h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
      #rk4 for legendre with step=h
      LOP_rk4 <- function(x0,y0){
        f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
        k1 <- f(x0,y0)
        k2 <- f(x0+h/2,y0+h/2*k1)
        k3 <- f(x0+h/2,y0+h/2*k2)
        k4 <- f(x0+h,y0+h*k3)
        y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
        return(y)
      }
      #dy_LOP, the increasment of each step
      dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[c],0))
      #dy_LOP*y= main or sub effect
      dy_fit <- effect*c(0,dy[1:(length(times)-1)])
      return(cumsum(dy_fit))
    }
    ode_optimize <- function(pars,ind,dep,times,data,order){
      ind_pars <- matrix(pars,ncol=order)[1,]
      dep_pars <- matrix(pars,ncol=order)[-1,]
      ind_effect <- get_effect(ind_pars,data[,ind],times,order)
      if ( is.null(nrow(dep_pars)) ) {
        dep_effect <- get_effect(dep_pars,data[,dep],times,order)
        y <- ind_effect+dep_effect+data[,ind][1]
      }else{
        dep_effect <- sapply(1:length(dep), function(c)
          get_effect(dep_pars[c,],data[,dep[c]],times,order))
        y <- ind_effect+rowSums(dep_effect)+data[,ind][2]
      }
      #ssr <- sum((data[,ind]-y)^2)
      #add penalty
      alpha=alpha
      ridge <- sum((data[,ind]-y)^2+alpha*(sum(ind_pars^2)+sum(dep_pars^2)))
      return(ridge)
    }
    get_value <- function(effect,data,times,order){
      #input
      ind <- data[[1]]
      dep <- data[[2]]
      ind_no <- as.numeric(which(colnames(effect)==ind))
      dep_no <- as.numeric(sapply(1:length(dep), function(c) which(colnames(effect)==dep[c])))
      init_pars <- rep(0.001,(length(ind_no)+length(dep_no))*order)
      result <- optim(init_pars,ode_optimize,ind=ind_no,dep=dep_no,
                      times=times,data=effect,order=order,
                      method = "BFGS", control=list(maxit=1000,trace=T))
      par_after <- matrix(result$par,length(ind)+length(dep),order)
      return(par_after)
    }
    core.number <- detectCores()
    cl <- makeCluster(getOption("cl.cores", core.number))
    clusterEvalQ(cl, {library(orthopolynom)})
    clusterExport(cl, c("get_value","data","ode_optimize","get_effect","times","get_interaction",
                        "cluster_mean","module_relationship","order"),envir=environment())
    lop_par <- pblapply(1:nrow(data),function(c)get_value(t(cluster_mean),
                                                          module_relationship[[c]],times,order),cl=cl)
    stopCluster(cl)
    return(list(lop_par,module_relationship))
  }
  all_lop_par <- get_legendre_par(exp_index,legendre_order,data,alpha)
  get_output <- function(relationship,par,effect,times,order){
    get_effect <- function(pars,effect,times,order){
      if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
      LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
      d_LOP_fit <-  sapply(1:length(pars),function(c)
        pars[c]*polynomial.derivatives(LOP)[[c+1]])
      h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
      #rk4 for legendre with step=h
      LOP_rk4 <- function(x0,y0){
        f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
        k1 <- f(x0,y0)
        k2 <- f(x0+h/2,y0+h/2*k1)
        k3 <- f(x0+h/2,y0+h/2*k2)
        k4 <- f(x0+h,y0+h*k3)
        y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
        return(y)
      }
      #dy_LOP, the increasment of each step
      dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[c],0))
      #dy_LOP*y= main or sub effect
      dy_fit <- effect*c(0,dy[1:(length(times)-1)])
      return(cumsum(dy_fit))
    }
    output <- list()
    output[[1]] <- relationship[[1]]
    output[[2]] <- relationship[[2]]
    output[[3]] <- par[1,]
    output[[4]] <- par[2:nrow(par),]
    ind_no <- as.numeric(which(colnames(effect)==output[[1]]))
    dep_no <- as.numeric(sapply(1:length(output[[2]]),
                                function(c) which(colnames(effect)==output[[2]][c])))
    inital_value <- effect[,ind_no][1]
    ind_effect <- get_effect(as.numeric(output[[3]]),effect[,ind_no],times,order)+inital_value
    if (length(dep_no)==1) {
      dep_effect <- get_effect(as.numeric(output[[4]]),effect[,dep_no],times,order)
    }else{
      dep_effect <- sapply(1:length(dep_no), function(c)
        get_effect(as.numeric(output[[4]][c,]),effect[,dep_no[c]],times,order))
      colnames(dep_effect) <- dep_no
    }
    #------------
    all_effect <- cbind(ind_effect,dep_effect)
    #effect_mean <- all_effect[5,]
    effect_mean <- apply(all_effect,2,mean)
    output[[5]] <- effect_mean
    output[[6]] <- all_effect
    return(output)
  }
  get_net <- function(data){
    times <- exp_index
    cluster_mean <- data
    module_relationship <- all_lop_par[[2]]
    net <- pblapply(1:nrow(data),function(c)
      get_output(module_relationship[[c]],all_lop_par[[1]][[c]],
                 t(cluster_mean),times=times,order=legendre_order))
    return(net)
  }
  all_net <- get_net(data)
  after_fn<-function(all_net_num){
    get_after <- function(i){
      temp <- matrix(NA,nrow = length(i[[2]]),ncol=3)
      temp[,1] <- i[[2]]
      temp[,2] <- i[[1]]
      temp[,3] <- i[[5]][2:(length(i[[2]])+1)]
      
      colnames(temp) <- c('from','to','dep_effect')
      temp <- data.frame(temp)
      temp[,3] <- as.numeric(as.character(temp[,3]))
      return(temp)
    }
    after <- do.call(rbind,lapply(all_net_num,get_after))
    after[,4] <-ifelse(after[,3]>=0,"+","-")
    colnames(after)[4] <- c("color")
    after[,3] <- abs(after[,3])
    return(after)
  }
  after <- after_fn(all_net)
  reslut <- list(all_lop_par=all_lop_par,all_net=all_net,after=after)
  return(reslut)
}
