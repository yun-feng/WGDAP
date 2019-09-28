#' Calculate the Separation Value for Ploidy
#'
#' In order to make the results more biologically meaningful, we want to set the (haplotype) ploidy to either 1 or 2, 
#' indicating whether such sample has genome duplication or not.
#' This function finds the valley for the distribution of inferref ploidy from FPI, and then we could separate the ploidy based on this value (see \code{\link[graphics]{int_ploidy}}).
#'
#' @param var The current value for all variables (see \code{\link[graphics]{init_var}})
#'
#' @return Value for the separation of ploidy
#'
#' @examples
#' 
#' g_sep=get_sep(var)
#'
#' @export
#' 
#'

get_sep<-function(var){
  inf_g=var$inf_g
  h<-density(inf_g)
  for(i in 2:length(h$y)){
    if(h$y[i]<h$y[i-1] && h$y[i]<h$y[i+1] ) break
  }
  g_sep<-h$x[i]
  return(g_sep)
}