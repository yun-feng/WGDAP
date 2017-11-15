#' Plot Recurrent CNV Before v.s. After WGD for each loci
#'
#'Plot the values for the common events happenning before v.s. after WGD for each loci, by scatterplot
#' ggplot2 is needed for this function
#'
#' @param data the CNV dataset (see \code{\link[graphics]{load_data}})
#' @param var The optimal value for all variables (see \code{\link[graphics]{Run_lasso}})
#' @param loci The particular loci to be emphasized
#' 
#' @return ggplot2 graph object
#'
#'
#' @examples
#' 
#' library(ggplot2)
#' p_scatter<-plot_scatter(wkdata,Lasso_res$var,loci=c(3606,864,4858,5261,5570))
#'
#' @export
#' 
#'
plot_scatter<-function(data,var,loci=NULL){
  
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
  p_c<-ggplot(data=wkdata)+geom_point(aes(x=Before,y=After,color=Chr),size=1.3)
  if(length(loci)){
    p_c<-p_c+
      geom_rect(data=wkdata[loci,],aes(xmin=Before-0.015,xmax=Before+0.015,ymin=After-0.02,ymax=After+0.02),color="black",fill=NA,show.legend=F,size=1.0)+
      geom_text(data=wkdata[loci,],aes(x=Before,y=After,label=paste("Locus",SNP_locus)),nudge_y=-0.025,size=4)
  }
  
  p_c<-p_c+ labs(x = "Before WGD",y = "After WGD")+
    theme(axis.text = element_text(color = "blue", size = 20))+
    theme(axis.title.x = element_text(color = "black", size = 20))+
    theme(axis.title.y = element_text(color = "black", size = 20))+
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold")
    )
  return(p_c)
}