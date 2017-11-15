#' Quality Control for Post-processing
#'
#' Check the dual problem of the lasso problem helps us to know how far we are from optimal.
#' 
#' This function would print out 5 values.
#' 
#' The first 2 are the dual variable and the 3rd and 4th should be your lambda parameters, which were used for \code{\link[graphics]{Run_lasso}}.
#' To achieve optimality, we should have the first 2 smaller or equal to the 3rd and 4th respectively.
#' However, due to numeric error, sometimes it is slightly larger (less than 0.01%), which is also acceptable.
#' 
#' The last one is a duality gap, but it is a valid value only when the requirement above is satisfied.
#' And it is strictly positive when it is valid.
#' However, it should be a small value in most cases.
#' 
#'
#' @param data The CNV dataset prepared for lasso (see \code{\link[graphics]{data_lasso}})
#' 
#' @param par The parameter prepared for lasso (see \code{\link[graphics]{par_lasso}})
#' 
#' @param res The result after solving lasso problem  (see \code{\link[graphics]{Run_lasso}})
#' 
#' @param g Integer-valued ploidy
#' 
#' @param lambda L1-Penalty for all the variables (see \code{\link[graphics]{Run_lasso}})
#' 
#' @return None
#'
#' @examples
#' 
#' Check_dual(wkdata,par,Lasso_res,g_int)
#'
#' @export
#' 
#'

Check_dual<-function(data,par,res,g,lambda=c(8,4)){
  C_m=data$minor
  C_p=data$major
  K=data$nsample
  
  inf_X=res$var
  
  par_a=par$par_a
  par_A=par$par_A
  par_B=par$par_B
  A=par$A
  X_mean=par$X_mean
  
  new=res$loss_array[length(res$loss_array)]
  
  for (i in 1:2){print(max((0.5*C_m+0.5*C_p-inf_X%*%par_A%*%t(g)-inf_X%*%par_B%*%t(rep(1,K)))%*%A[i,]));
    print(lambda[i]);}
  
  dual<-1/4*norm(C_m-X_mean%*%par_A%*%t(g)-X_mean%*%par_B%*%t(rep(1,K)),type="F")^2+1/4*norm(C_p-X_mean%*%par_A%*%t(g)-X_mean%*%par_B%*%t(rep(1,K)),type="F")^2-1/2*norm((inf_X-X_mean)%*%par_A%*%t(g)+(inf_X-X_mean)%*%par_B%*%t(rep(1,K)),type="F")^2
  gap=new-dual
  print(gap)
  
}