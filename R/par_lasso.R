#' Prepare Parameters for Post-processing
#'
#' Prepare all the parameters used in lasso problem.
#' 
#' In this function, we allow users to change the type of CNV happenning before and after WGD,
#' but these changes should be the same for minor and major copies, for this makes interpretation of these changes much easier.
#' 
#' The default one is set to changes before WGD affecting all samples and changes after WGD only affecting samples without WGD.
#' The mean for these changes are set to 0 and -0.5 respectively, for in most samples with WGD are likely to lose one copy after WGD. 
#' These should be the best settings for most datasets (if minor and major copies are not considered seperately), adding or changing the types will usually make interpretation much harder.
#' 
#' This function also calculate a multiplicative matrix and step length used in optimization.
#' 
#'
#' @param N Number of loci
#' @param K Number of samples
#' @param g Integer-valued ploidy
#' @param par_a The constant indicating how ploidy infecting the final copy number profile
#' @param mean_X A vector containing mean values of all the recurrent CNV
#' @param a A vector containing the coeffient of the effect of the recurrent CNV with respect to the ploidy value
#' @param b A vector containing the coeffient of the effect of the recurrent CNV with respect to the constant value
#' 
#' 
#' @return a list containing all the parameters used in solving lasso problem
#'
#'
#' @examples
#' 
#' par<-par_lasso(wkdata$nloci,wkdata$nsample,g_int)
#' 
#' ##Add CNV that only affct samples without WGD to the model.
#' 
#' par<-par_lasso(wkdata$nloci,wkdata$nsample,g_int,mean_X=c(0,-0.5,0),a=c(1,1,-1),b=c(0,-1,2))
#'
#' @export
#' 
#'

par_lasso<-function(N,K,g,
                    par_a=1,
                    mean_X=c(0,-0.5),
                    a=c(1,1),
                    b=c(0,-1)){
  X_mean<-rep(mean_X,each=N)
  dim(X_mean)<-c(N,2)
  
  #Multiplier
  par_A<-matrix(a,ncol=1)
  par_B<-matrix(b,ncol=1)
  
  A<-par_A%*%t(g)+par_B%*%t(rep(1,K))
  steplen<-diag(1/(A%*%t(A)))
  
  return(list(par_a=par_a,
              mean_X=mean_X,
              X_mean=X_mean,
              par_A=par_A,
              par_B=par_B,
              A=A,
              steplen=steplen)
  )
}