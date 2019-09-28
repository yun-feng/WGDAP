#' Fixed Point Iteration
#'
#' Run fixed point iteration.
#' 
#' This is the main part for WGD analysis. 
#' This function uses the fixed point method to iteratively calculate the optimal value for all variables and store the loss for each iteration in an array.
#'
#' It runs until the improvement in loss is sufficiently small or the number of iteration reaches the max defined by user.
#' It then returns a list, containing two objects: var, which is a list for the optimal value of the variables; loss_array, which is an array of loss values.
#'
#' @param data The CNV dataset (see \code{\link[graphics]{load_data}})
#' 
#' @param par The parameter used for fixed point method. (see \code{\link[graphics]{set_par}})
#' 
#' @param var The current value for all variables (see \code{\link[graphics]{init_var}})
#' 
#' @param new The initial loss value
#' 
#' @param max_iter The maximum number of iterations that the method is allowed to run
#'
#' @return list of the optimal values for variables and the loss for each iteration
#'
#'
#' @examples
#' 
#' FPI_result=FPI(wkdata,CNV_par,common_var,init_loss)
#'
#' @export
#' 
#'


FPI<-function(data,par,var,new,max_iter=400){
  remain=c(new)
  
  i=0
  last=new+1;
  
  C_m=data$minor
  C_p=data$major
  N=data$nloci
  K=data$nsample
  
  inf_g=var$inf_g
  inf_X_m=var$inf_X_m
  inf_X_p=var$inf_X_p
  
  par_a_m=par$par_a_m
  par_A_m=par$par_A_m
  par_B_m=par$par_B_m
  mean_m= par$mean_m
  penal_m=par$penal_m
  
  par_a_p=par$par_a_p
  par_A_p=par$par_A_p
  par_B_p=par$par_B_p
  mean_p=par$mean_p
  penal_p=par$penal_p
  
  penal_g=par$penal_g
  weight_L=par$weight_L
  
  while(abs(new-last)>10^(-5)){
    
    last=new;
    
    temp_C_m=C_m-par_a_m*rep(1,N)%*%t(inf_g)
    temp_C_m=temp_C_m%*%inf_g%*%t(par_A_m)+temp_C_m%*%rep(1,K)%*%t(par_B_m)
    temp_C_m=temp_C_m+matrix(rep(mean_m*penal_m,each=N),nrow=N,ncol=length(mean_m))
    temp_m=diag(penal_m)+par_A_m%*%(t(inf_g)%*%inf_g)%*%t(par_A_m)+par_A_m%*%(t(inf_g)%*%rep(1,K))%*%t(par_B_m)
    temp_m=temp_m+par_B_m%*%(t(rep(1,K))%*%rep(1,K))%*%t(par_B_m)+par_B_m%*%(t(rep(1,K))%*%inf_g)%*%t(par_A_m)
    inf_X_m=temp_C_m%*%solve(temp_m)
    
    temp_C_p=C_p-par_a_p*rep(1,N)%*%t(inf_g)
    temp_C_p=temp_C_p%*%inf_g%*%t(par_A_p)+temp_C_p%*%rep(1,K)%*%t(par_B_p)
    temp_C_p=temp_C_p+matrix(rep(mean_p*penal_p,each=N),nrow=N,ncol=length(mean_p))
    temp_p=diag(penal_p)+par_A_p%*%(t(inf_g)%*%inf_g)%*%t(par_A_p)+par_A_p%*%(t(inf_g)%*%rep(1,K))%*%t(par_B_p)
    temp_p=temp_p+par_B_p%*%(t(rep(1,K))%*%rep(1,K))%*%t(par_B_p)+par_B_p%*%(t(rep(1,K))%*%inf_g)%*%t(par_A_p)
    inf_X_p=temp_C_p%*%solve(temp_p)
    
    temp_C=t(C_m)%*%(weight_L%*%(inf_X_m%*%par_A_m+par_a_m*rep(1,N)))+t(C_p)%*%(weight_L%*%(inf_X_p%*%par_A_p+par_a_p*rep(1,N)))
    temp_C=temp_C-rep(1,K)%*%((t(par_B_m)%*%t(inf_X_m))%*%weight_L%*%(inf_X_m%*%par_A_m+par_a_m*rep(1,N)))-rep(1,K)%*%((t(par_B_p)%*%t(inf_X_p))%*%weight_L%*%(inf_X_p%*%par_A_p+par_a_p*rep(1,N)))
    temp_C=temp_C+penal_g*rep(1,K)
    temp_g=t((inf_X_m%*%par_A_m+par_a_m*rep(1,N)))%*%weight_L%*%(inf_X_m%*%par_A_m+par_a_m*rep(1,N))+t((inf_X_p%*%par_A_p+par_a_p*rep(1,N)))%*%weight_L%*%(inf_X_p%*%par_A_p+par_a_p*rep(1,N))
    temp_g=temp_g+penal_g
    inf_g=temp_C/temp_g[1,1]
    
    left_m=C_m-par_a_m*rep(1,N)%*%t(inf_g)-inf_X_m%*%par_A_m%*%t(inf_g)-inf_X_m%*%par_B_m%*%t(rep(1,K))
    left_p=C_p-par_a_p*rep(1,N)%*%t(inf_g)-inf_X_p%*%par_A_p%*%t(inf_g)-inf_X_p%*%par_B_p%*%t(rep(1,K))
    new=sum(diag(crossprod(left_m)))+sum(diag(crossprod(left_p)));
    remain=c(remain,new)
    
    i=i+1
    if(i>max_iter){break}
  }
  return(list(var=list(inf_g=inf_g,inf_X_m=inf_X_m,inf_X_p=inf_X_p),loss_array=remain))
}
