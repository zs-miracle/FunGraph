#' @title plot_Manhattan 
#' @importFrom stats coef cor cov kmeans lm optim sd time
#' @param funmap_data_snp  Genotypic data.only accept numeric value for calculation, SNP genotype data should be converted into 0,1 matrix.
#' @param LR LR values for each gene locus.
#' @param threshold_value The threshold calculated by the Permutaion function.
#' @param x_axis_name The name of the X-axis can be customized.
#' @return Manhattan plot
#' @export
plot_Manhattan_fn <- function(funmap_data_snp,LR,threshold_value,
                              x_axis_name="S. aureus Genome (Mb)")
{
  tmp_Manhattan<- data.frame(as.numeric(as.numeric(rownames(funmap_data_snp))/1000),LR)
  sig_Manhattan<- data.frame(tmp_Manhattan[which(LR>=threshold_value),])
  unsigsig_Manhattan<- data.frame(tmp_Manhattan[which(LR<1.3&LR>1.298),])
  colnames(tmp_Manhattan) <- c(x_axis_name,"LR")
  ggplot(tmp_Manhattan)+
    geom_point(aes(x=tmp_Manhattan[,1],y=LR),color="#00A3CC")+
    geom_point(data=unsigsig_Manhattan,aes(x=unsigsig_Manhattan[,1],y=LR),color="#CC0066")+
    geom_point(data=sig_Manhattan,aes(x=sig_Manhattan[,1],y=LR),color="red")+
    labs(x=x_axis_name)+scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0),breaks = c(0,20,40,60),limits = c(0,80))+
    geom_hline(yintercept =threshold_value,linetype=2,color="red")+
    theme_classic()+theme(panel.background = element_rect(color = 'black', fill = 'transparent'),
                          axis.title =element_text(size=15),
                          axis.text =element_text(size=12))
}