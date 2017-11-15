#' Plot Inferred v.s. Average Ploidy
#'
#' Plot the inferred ploidy from FPI to average ploidy.
#' ggplot2 is needed for this function
#'
#' @param data the CNV dataset (see \code{\link[graphics]{load_data}})
#' @param var The optimal value for all variables (see \code{\link[graphics]{init_var}})
#' @param g_sep The separation value for inferred ploidy
#' 
#' 
#' @return ggplot2 graph object
#'
#'
#' @examples
#' 
#' library(ggplot2)
#' p<-plot_ploidy(wkdata,FPI_result$var,g_sep)
#'
#' @export
#' 
#'


plot_ploidy<-function(data,var,g_sep){
  
  C_m=data$minor
  C_p=data$major
  K=data$nsample
  N=data$nloci
  
  inf_g=var$inf_g
  
  ave_ploidy<-c()
  for(i in 1:K){
    ave_ploidy<-c(ave_ploidy,(sum(C_m[,i])+sum(C_p[,i]))/N)
  }
  
  wkdata<-data.frame(ave_ploidy,inf_ploidy=inf_g)
  fit<-lm(inf_g~ave_ploidy)
  
  p_a<-ggplot(wkdata)+geom_abline(slope=0.5,intercept=0,color="blue",size=2.2,alpha=0.6)+geom_point(aes(x=ave_ploidy,y=inf_ploidy),alpha=0.8,size=2)
  p_a<-p_a+geom_abline(slope=0,intercept=g_sep,color="gold",lty=2,size=2)
  p_a<-p_a+geom_abline(slope=fit$coefficients[2],intercept=fit$coefficients[1],color="darkgreen",size=2.7,alpha=0.6)+
    theme(axis.text = element_text(color = "blue", size = 22))+
    theme(axis.title.x = element_text(color = "black", size = 25))+
    theme(axis.title.y = element_text(color = "black", size = 25))+
    labs(x = "Average Ploidy",y = "Predicted Ploidy")
  return(p_a)
}