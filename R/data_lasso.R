#' Prepare Dataset for Post-processing
#'
#' Modify the original dataset in order to save some steps in the approximal gradient descent iteration for lasso.
#'
#' @param data The CNV dataset (see \code{\link[graphics]{load_data}})
#' 
#' @param g The integer-valued ploidy
#' 
#' @param par_a The constant indicating how ploidy infecting the final copy number profile
#' 
#' @return List of values or matrice for lasso problem
#'
#'
#' @examples
#' 
#' wkdata<-data_lasso(wkdata,g_int)
#'
#' @export
#' 
#'

data_lasso<-function(data,g,par_a=1){
  C_m<-as.matrix(data$minor)
  C_p<-as.matrix(data$major)
  N=nrow(C_m)
  K=ncol(C_m)
  
  C_m=C_m-par_a*rep(1,N)%*%t(g)
  C_p=C_p-par_a*rep(1,N)%*%t(g)
  
  return(list(minor=C_m,major=C_p,nsample=K,nloci=N, loci=data$loci))
}