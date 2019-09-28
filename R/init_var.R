#' Set the Initial Value for Variables
#' 
#' In order to run the fixed-point method, we need to start from some initial value.
#' This function set the initial values for all the variables to their mean values.
#'
#' @param data The CNV dataset (see \code{\link[graphics]{load_data}})
#' 
#' @param par The parameter used for fixed point method. (see \code{\link[graphics]{set_par}})
#'
#' @return List of three matrice containing initial values for ploidy and recurrent CNV before and after WGD.
#'
#'
#' @examples
#' 
#' common_var<-init_var(wkdata,CNV_par)
#'
#' @export
#'

init_var<-function(data,par){
  v<-list( inf_g=matrix(rep(1,data$nsample),ncol=1),
           inf_X_m=matrix(rep(par$mean_m,each=data$nloci),nrow=data$nloci,ncol=length(par$mean_m)),
           inf_X_p=matrix(rep(par$mean_p,each=data$nloci),nrow=data$nloci,ncol=length(par$mean_p)) 
  )
  return(v)
}