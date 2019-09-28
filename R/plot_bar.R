#' Plot the Original Copy number for Selected Loci
#'
#' To study some particular loci, this function count the number of samples 
#' containing certain number of loci in samples with and without WGD
#'
#' @param data The CNV dataset (see \code{\link[graphics]{load_data}})
#' @param g Inferred ploidy value (see \code{\link[graphics]{}})
#' @param loci The particular loci to be emphasized
#' @param max The max number of copies for all loci 
#' 
#' @return ggplot2 graph object
#'
#'
#' @examples
#' 
#' library(ggplot2)
#' p_bar<-plot_bar(wkdata,g_int,loci=c(3606,864,4858,5261,5570),max=7)
#'
#' @export
#' 
#'

plot_bar<-function(data,g,loci,max){
  C_m=data$minor
  C_p=data$major
  
  g<-rep(g,2)
  C_phase=cbind(C_p,C_m)
  num<-c()
  ploidy<-c()
  CNV<-c()
  locus<-c()
  for(i in 0:max){
    for(j in 1:length(loci)){
      num<-c(num,length(which(C_phase[loci[j],which(g<1.5)]==i)))
      CNV<-c(CNV,i)
      ploidy<-c(ploidy,1)
      locus<-c(locus,loci[j])
      
      num<-c(num,length(which(C_phase[loci[j],which(g>1.5)]==i)))
      CNV<-c(CNV,i)
      ploidy<-c(ploidy,2)
      locus<-c(locus,loci[j])
    }
    
  }
  loci_comp<-data.frame(Number_of_samples=num,ploidy,CNV,SNP_locus=factor(locus))
  library(ggplot2)
  p_d<-ggplot(data=loci_comp)+geom_bar(aes(x=CNV,y=Number_of_samples,fill=SNP_locus),stat="identity",pos="dodge")+
    facet_wrap(~ploidy)+theme(axis.text = element_text(color = "blue", size = 20))+
    theme(axis.title.x = element_text(color = "black", size = 20))+
    theme(axis.title.y = element_text(color = "black", size = 20))+
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20, face = "bold"),
      strip.text = element_text(size=20)
    )+
    labs(fill="SNP locus",x="CNV",y="Number of Samples")
  return(p_d)
}