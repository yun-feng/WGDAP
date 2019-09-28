#' Prepare Variables for Post-processing
#'
#' Set up initial values for variables used for lasso problem.
#' 
#' Since in most cases, fixed-point method is run before post-processing. 
#' It saves a bit time if we start from that result. 
#' We can just replicate that result by the factor which was used before during down-sampling, to get an initial value for post-processing.
#'  
#'
#' @param var The value for all variables from FPI (see \code{\link[graphics]{init_var}})
#' 
#' @param factor The factor to which you want to up-sample the variable from FPI (should be the same as the down-sample factor used before)
#' 
#' @return matrice containing the recurrent CNV before and after WGD
#'
#'
#' @examples
#' 
#' var<-var_lasso(FPI_result$var)
#'
#' @export
#' 
#'

var_lasso<-function(var,factor=10){
  inf_X_m<-as.matrix(var$inf_X_m)
  inf_X_p<-as.matrix(var$inf_X_p)
  inf_X=rep(0.5*inf_X_m+0.5*inf_X_p,each=factor)
  N=length(inf_X)/2
  dim(inf_X)<-c(N,2)
  return(inf_X)
}