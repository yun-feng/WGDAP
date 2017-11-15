#' Calculate Integer Ploidy Value for Post-processing
#'
#' Set (haplotype) ploidy value to either 1 or 2.
#'
#' To get the interger ploidy value, 
#' we set all the infered ploidy value from FPI to be either 1 or 2, 
#' based on their value relative to the seperation value.
#' 
#' This gives us (haplotype) ploidy values for all samples, based on which we could run a lasso-like post-processing step.
#' The post-processing could make the prediction for recurrent CNV before and after WGD more biologically meaningful.
#' 
#' @param var The current value for all variables (see \code{\link[graphics]{init_var}})
#'
#' @param sep The separation value (see \code{\link[graphics]{get_sep}})
#'
#' @return Cector containing integer-valued ploidy
#'
#' @examples
#' 
#' g_int=int_ploidy(FPI_result$var,g_sep)
#'
#' @export
#' 
#'

int_ploidy<-function(var,sep){
  g<-rep(1,length(var$inf_g))
  g[which(var$inf_g>sep)]<-2
  return(g)
}