#' Post-processing
#'
#' Solve the lasso-like problem for post-processing
#' 
#' This step makes the solution more biologically interpretable. 
#' And the lasso-like L1 regularizer makes the result sparser which is crucial for researchers to find region of interest.
#' 
#' This function utilizes a proximal gradient descent method to solve the problem,
#' and it stops when the improvement is sufficiently small and the dual variable is close to feasible (or it reaches the max number of iterations specified by user).
#' 
#' In post-processing, we don't have separate values for minor and major copies, 
#' and thus, the result would be an intermediate between changes for minor and major copy.
#' 
#' We find this easier to understand, but if the user indeed wants to have a set of separate values for minor and major copy,
#' he could copy major copy data file to minor copy and re-load the dataset.
#' Namely, he could change major and minor copy data file to both contain either minor or major data, and run the lasso.
#' 
#' @param data The CNV dataset prepared for lasso (see \code{\link[graphics]{data_lasso}})
#' 
#' @param par The parameter prepared for lasso (see \code{\link[graphics]{par_lasso}})
#' 
#' @param var The initialized variables for lasso  (see \code{\link[graphics]{var_lasso}})
#' 
#' @param g Integer-valued ploidy
#' 
#' @param lambda L1-Penalty for all recurrent CNV
#' 
#' @param max The maximum number of cycles to run
#' 
#' @return A list containing the optimal values for variables and the loss array
#'
#'
#' @examples
#' 
#' Lasso_res<-Run_lasso(wkdata,par,var,g_int)
#'
#' @export
#' 
#'

Run_lasso<-function(data,par,var,g,lambda=c(8,4),max=1000){
  
  C_m=data$minor
  C_p=data$major
  K=data$nsample
  N=data$nloci
  
  inf_X=var
  
  par_a=par$par_a
  par_A=par$par_A
  par_B=par$par_B
  A=par$A
  X_mean=par$X_mean
  mean_X=par$mean_X
  steplen=par$steplen
  
  error_array=c()
  
  new=1/4*norm(C_m-inf_X%*%par_A%*%t(g)-inf_X%*%par_B%*%t(rep(1,K)),type="F")^2+1/4*norm(C_p-inf_X%*%par_A%*%t(g)-inf_X%*%par_B%*%t(rep(1,K)),type="F")^2+sum(abs(inf_X-X_mean)%*%diag(lambda))
  error_array=c(new)
  old=new+1
  
  cycle=0
  while(abs(new-old)>10^(-5) || flag){
    for(i in 1:length(steplen)){
      flag=F
      temp=(0.5*C_m+0.5*C_p-inf_X%*%par_A%*%t(g)-inf_X%*%par_B%*%t(rep(1,K)))%*%A[i,]
      if(max(abs(temp))>lambda[i]*1.000001){flag=T}
      inf_X[,i]=inf_X[,i]+steplen[i]*temp
      off_set=inf_X[,i]-mean_X[i]
      temp<-abs(off_set)-lambda[i]*steplen[i]
      inf_X[,i]<-ifelse(temp<0,mean_X[i],mean_X[i]+temp*sign(off_set))
    }
    old=new
    new=1/4*norm(C_m-inf_X%*%par_A%*%t(g)-inf_X%*%par_B%*%t(rep(1,K)),type="F")^2+1/4*norm(C_p-inf_X%*%par_A%*%t(g)-inf_X%*%par_B%*%t(rep(1,K)),type="F")^2+sum(abs(inf_X-X_mean)%*%diag(lambda))
    error_array=c(error_array,new)
    cycle=cycle+1;
    if(cycle>max){break}
  }
  return(list(var=inf_X,loss_array=error_array))
}
