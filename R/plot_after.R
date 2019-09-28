#' Plot Recurrent CNV After WGD for each loci
#'
#' Plot the values for the common events happenning after WGD for each loci
#' ggplot2 is needed for this function
#'
#' @param data the CNV dataset (see \code{\link[graphics]{load_data}})
#' @param var The optimal value for all variables (see \code{\link[graphics]{Run_lasso}})
#' 
#' 
#' @return ggplot2 graph object
#'
#'
#' @examples
#' 
#' library(ggplot2)
#' p_after<-plot_after(wkdata,Lasso_res$var)
#'
#' @export
#' 
#'

plot_after<-function(data,var){
  
  coordinate_m=data$loci
  N=data$nloci
  inf_X=var
  
  inf_X<-as.matrix(inf_X)
  
  chr_c<-c()
  for (i in 2:N){
    if(coordinate_m[i,1]>coordinate_m[i-1,1]){
      chr_c<-c(chr_c,i)
    }
  }
  chr_center=c(0,chr_c[1:length(chr_c)])
  chr_center=(chr_center+c(chr_c,N))/2
  chrmosome=data.frame(chr_center,Chr=factor(1:length(chr_center)))
  
  wkdata<-data.frame(SNP_locus=1:N,Before=2*inf_X[,1],After=2*inf_X[,2],Chr=factor(coordinate_m[,1]))
  p_b<-ggplot(data=wkdata)+geom_point(aes(x=SNP_locus,y=After,color=Chr),size=1.3,show.legend = F)+
    geom_text(data=chrmosome,aes(x=chr_center,y=rep(min(2*inf_X[,2]*1.001),length(chr_center)),label=c(1:length(chr_center)),color=Chr),show.legend = F,size=5)+
    geom_vline(xintercept =chr_c-0.5)
  p_b<-p_b+
    labs(x="SNP locus", y = "After WGD")+theme(axis.text = element_text(color = "blue", size = 20))+
    theme(axis.title.x = element_text(color = "black", size = 20))+
    theme(axis.title.y = element_text(color = "black", size = 20))+
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold")
    )
  p_b
  return(p_b)
}