## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  tidy = TRUE,
  cache = TRUE,
  collapse = TRUE,
  dev = "png",
  fig.width = 7, 
  fig.height =3.5
)

## ----load-packages, include=TRUE----------------------------------------------
library(FunGraph)

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  View(s)

## ----echo=TRUE----------------------------------------------------------------
knitr::kable(head(remarker_data_par[,1:8]))

## ----echo=TRUE----------------------------------------------------------------
knitr::kable(head(remarker_effect[,1:11]))

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  View(fourth_effect_cluster)

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  View(fourteen_cluster_mean_effect)

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  View(S_mo)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(S_mo[,1:11]))

## ---- eval=FALSE--------------------------------------------------------------
#  View(S_SNP)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(S_SNP[,1:10]))

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  s <- Funmap_function(funmap_data=S_mo,funmap_data_snp=S_SNP,times=c(0.5,1,1.5,2,4,6,8,10,12,16,20))

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  permutaion(funmap_data=S_mo,funmap_data_snp=S_SNP,
#                  times=c(0.5,1,1.5,2,4,6,8,10,12,16,20),
#                  permutaion_times=1000,rownumber=1000)
#  

## ----echo=TRUE----------------------------------------------------------------
plot_Manhattan_fn(funmap_data_snp=S_SNP,LR=s$LR,threshold_value=54.33574)

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  remarker_data_par <- t(sapply(1:length(s$LR),function(c)s$All_H1_result[c][[1]]$H1_par))
#  remarker_effect <- Get_effect(S_SNP,remarker_data_par,times=c(0.5,1,1.5,2,4,6,8,10,12,16,20))

## ----echo=TRUE----------------------------------------------------------------
library(ggplot2)
library(patchwork)
sig_SNP <- S_SNP[which(s$LR>54.33574),]
sig_data_par <- t(sapply(c(which(s$LR>54.33574)),function(c)s$All_H1_result[c][[1]]$H1_par))
sig_effect <- Get_effect(sig_SNP,sig_data_par,times=c(0.5,1,1.5,2,4,6,8,10,12,16,20))
effect_plot_fn <- function(effect_data,colnn,times){
  tmp <- as.data.frame(t(rbind(effect_data,times=times)))
  p1 <- ggplot(tmp)+geom_line(data=tmp,aes(x=times,y=tmp[,colnn]),size=2,color="#8CB8D4")+
    theme_classic()+labs(y=NULL,x=NULL)+
    theme(panel.background = element_rect(color = 'black', fill = 'transparent'),
          plot.margin = unit(c(0,0,0,0),"lines"))+
    scale_x_continuous(expand = c(0,0))+
    annotate("text",label=paste0("Q",rownames(effect_data)[colnn]),
             x=8.5,y=0.6,size=3.5)+
    scale_y_continuous(expand = c(0,0),limits = c(0,0.75))
  
  return(p1)
}
effect_plot_left_fn <- function(i){
  effect_plot_fn(
    effect_data=sig_effect,colnn=i,times=c(0,0.5,1,1.5,2,4,6,8,10,12,16,20))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.length.x = unit(-0.1,"cm"))
}
effect_plot_middle_fn <- function(i){
  effect_plot_fn(
    effect_data=sig_effect,colnn=i,times=c(0,0.5,1,1.5,2,4,6,8,10,12,16,20))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.length.x = unit(-0.1,"cm"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.length.y = unit(-0.1,"cm"))
}
effect_plot_down_fn <- function(i){
  effect_plot_fn(
    effect_data=sig_effect,colnn=i,times=c(0,0.5,1,1.5,2,4,6,8,10,12,16,20))+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.length.y = unit(-0.1,"cm"))
}

sig_effect_plot <- list()

for(i in c(1,6,11,16,21)){
  sig_effect_plot[[i]]<-effect_plot_left_fn(i)
}

for(i in c(2,3,4,5,7,8,9,10,12,13,14,15,17,18,19,20,22,23,24,25)){
  sig_effect_plot[[i]]<-effect_plot_middle_fn(i)
}
for(i in c(27,28,29,30)){
  sig_effect_plot[[i]]<-effect_plot_down_fn(i)
}
sig_effect_plot[[26]] <- effect_plot_fn(effect_data=sig_effect,colnn=26,
                                              times=c(0.5,1,1.5,2,4,6,8,10,12,16,20))

sig_effect_plot_all <- sig_effect_plot[[1]]
for(i in c(2:30)){
  sig_effect_plot_all <-sig_effect_plot_all+sig_effect_plot[[i]]
}
sig_effect_plot_all <- sig_effect_plot_all+plot_layout(ncol=5)
sig_effect_plot_all

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  bic_list <- list()
#  bic<- list()
#  for(k in 1:25){
#    bic <- Funcluster(remarker_effect,times=c(0.5,1,1.5,2,4,6,8,10,12,16,20),rep=1,k,legendre_order=4)
#    bic_list[[i]] <- bic
#  }
#  fourth_effect_cluster <- bic_list[[i]]

## ----echo=TRUE, warning=FALSE-------------------------------------------------
library(reshape2)
library(orthopolynom)
clustered_data <-fourth_effect_cluster
par1 <- clustered_data[[2]][[1]]$curve_par
times <- c(0.5,1,1.5,2,4,6,8,10,12,16,20)
clustered_df <- data.frame(cbind(row.names(clustered_data[[2]][[1]][["clustered_data"]]),clustered_data[[2]][[1]][["clustered_data"]]))
colnames(clustered_df)<- c("ID",times,"cluster")
long_df <- melt(clustered_df,id.vars=c("ID","cluster"))
long_df[,4] <- as.numeric(long_df[,4])
colnames(long_df) <- c("SNP","cluster","time","effect")
normalization <- function(x){((x-min(x))/(max(x)-min(x))*0.075)+0.02}

alpha <- 1/table(clustered_df$cluster)
alpha <- normalization(alpha)
get_mean_df <- function(data,times,legendre_order){
  legendre_fit <- function(par){
    x <- seq(min(times),max(times),length=30)
    fit <- sapply(1:length(par),function(c)
      par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
    legendre_fit <- as.matrix(as.data.frame(polynomial.values(
      polynomials=fit,x=scaleX(x, u=-1, v=1))))
    x_interpolation <- rowSums(legendre_fit)
    return(x_interpolation)
  }
  mean_df <- sapply(1:nrow(data),function(c)legendre_fit(as.numeric(data[c,])))
  colnames(mean_df) <- c(1:nrow(data))
  mean_df <- melt(mean_df)
  colnames(mean_df) <- c("time","cluster","effect")
  mean_df$time <- seq(min(times),max(times),length=30)
  return(mean_df)
}
mean_curve <- get_mean_df(par1,times=c(0.5,1,1.5,2,4,6,8,10,12,16,20),4)


#其中h是色相，范围越大，相邻颜色之间差异越大；c是饱和度，值越大色彩越浓艳饱满；l是亮度，大亮小暗。
library(RColorBrewer)
library(scales)

cols1 <- hue_pal(h = c(0, 720) + 15, c = 90, l = 60)(30)
darken <- function(color, factor=1.3){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
cluster_polt_fn<-function(i){
  cluster <- long_df[which(long_df$cluster==i),]
  cluster$time <- as.numeric(as.character(cluster$time))
  cluster_mean<- mean_curve[which(mean_curve$cluster==i),]
  p <- ggplot()+
    geom_line(cluster,mapping=aes(time,effect,group=SNP),color=cols1[i],alpha=alpha[i])+
    geom_line(cluster_mean,mapping=aes(time,effect),color=darken(cols1[i]),size=1.5)+
    theme_bw()+theme(panel.grid =element_blank())+
    xlab('')+ylab('')+
    scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0),limits = c(0,0.5))+
    annotate("text",label=paste0("M",unique(cluster$cluster),"(",dim(cluster)[1]/11,")"),
             x=10,y=0.45,size=5)+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))
  return(p)
}
cluster_polt <- list()
cluster_polt_left_fn <- function(i){
  p <- cluster_polt_fn(i)+theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.length.x = unit(-0.1,"cm"))
  return(p)
}

cluster_polt_middle_fn  <- function(i){
  p <- cluster_polt_fn(i)+theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.length.x = unit(-0.1,"cm"),
                                axis.title.y=element_blank(),
                                axis.text.y=element_blank(),
                                axis.ticks.length.y = unit(-0.1,"cm"))
  return(p)
}

cluster_polt_down_fn <- function(i){
  p <- cluster_polt_fn(i)+theme(axis.title.y=element_blank(),
                                axis.text.y=element_blank(),
                                axis.ticks.length.y = unit(-0.1,"cm"))
  return(p)
}
cluster_polt <- list()
for(i in c(1,6)){
  cluster_polt[[i]] <- cluster_polt_left_fn(i)
}

for(i in c(2,3,4,5,7,8,9,10)){
  cluster_polt[[i]] <- cluster_polt_middle_fn(i)
}

for(i in c(12,13,14)){
  cluster_polt[[i]] <- cluster_polt_down_fn(i)
}

cluster_polt[[11]] <- cluster_polt_fn(11)
cluster_polt_all <- cluster_polt[[1]]
  for (i in c(2:clustered_data[[1]][1])) {
    cluster_polt_all <- cluster_polt_all+cluster_polt[[i]]
  }

cluster_polt_all <- cluster_polt_all+plot_layout(ncol=5)
cluster_polt_all

## ----echo=TRUE, warning=FALSE-------------------------------------------------
  BIC <- c(-591045.6,-621432.9,-636578.2,-664730.9,-705232.3,
         -714223.3,-717334.2,-718324.2,-723442.1,-724443.2,
         -725674.3,-726223.2,-727881.3,-728135.5,-727281.0,
         -726923.3,-724106.0,-723333.7,-721561.9,-720992.1)
  BIC_data <- data.frame(bic=BIC,m = c(1:20))
  BIC_plot <- ggplot()+geom_line(data=BIC_data,aes(x=as.numeric(BIC_data$m),y=as.numeric(BIC_data$bic)),size=2)+
    theme_bw()+labs(x="No. of Modules",y="BIC")+
    geom_vline(xintercept =14,size=3,color="#F8766D",linetype="dotted")+
    scale_x_continuous(expand=c(0.05,0.05))+
    scale_y_continuous(expand=c(0.05,0.05))+
    theme(axis.text =element_text(size=8),
          axis.title = element_text(size = 12),
          plot.margin = unit(c(0,0,0,0),"lines"))+
    theme(legend.position = "none")+theme(panel.grid =element_blank())
  BIC_plot

## ----eval=FALSE, results="hide", fig.show="hide"------------------------------
#  times <- c(0.5,1,1.5,2,4,6,8,10,12,16,20)
#  clustered_data <-fourth_effect_cluster
#  par1 <- clustered_data[[2]][[1]]$curve_par
#  get_mean_effect <- function(data,times,legendre_order){
#    legendre_fit <- function(par){
#      x <- seq(min(times),max(times),length=30)
#      fit <- sapply(1:length(par),function(c)
#        par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
#      legendre_fit <- as.matrix(as.data.frame(polynomial.values(
#        polynomials=fit,x=scaleX(x, u=-1, v=1))))
#      x_interpolation <- rowSums(legendre_fit)
#      return(x_interpolation)
#    }
#    mean_df <- sapply(1:nrow(data),function(c)legendre_fit(as.numeric(data[c,])))
#    colnames(mean_df) <- c(1:nrow(data))
#    return(mean_df)
#  }
#  mean_effect <- t(get_mean_effect(par1,times,4))
#  fourteen_cluster_mean_effect <- Decomposition_curve(exp_index=seq(min(times),max(times),length=30),
#                      legendre_order = 4,data=mean_effect,alpha=1e-6)
#  fourteen_cluster_mean_effect

## ----echo=TRUE, warning=FALSE-------------------------------------------------
fourteen_cluster_mean_effect$after

## ----pressure, echo=FALSE, out.width = '100%'---------------------------------
knitr::include_graphics("D:/cytoscape.png")

