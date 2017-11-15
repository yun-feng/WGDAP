#' Evaluate Loss
#'
#' Evaluate current loss for fixed-point method based on current variable values.
#'
#' @param data the CNV dataset (see \code{\link[graphics]{load_data}})
#' 
#' @param par The parameter used for fixed point method. (see \code{\link[graphics]{set_par}})
#' 
#' @param var The current value for all variables (see \code{\link[graphics]{init_var}})
#'
#' @return Value for the current loss
#'
#' @examples
#' 
#' init_loss<-feval(wkdata,CNV_par,common_var)
#'
#' @export
#' 
#'

feval<-function(data,par,var){
  left_m=data$minor-par$par_a_m*rep(1,data$nloci)%*%t(var$inf_g)-var$inf_X_m%*%par$par_A_m%*%t(var$inf_g)-var$inf_X_m%*%par$par_B_m%*%t(rep(1,data$nsample))
  left_p=data$minor-par$par_a_p*rep(1,data$nloci)%*%t(var$inf_g)-var$inf_X_p%*%par$par_A_p%*%t(var$inf_g)-var$inf_X_p%*%par$par_B_p%*%t(rep(1,data$nsample))
  loss=sum(diag(crossprod(left_m)))+sum(diag(crossprod(left_p)));
  return(loss)
}