#' Set Up Parameters for Analysis
#'
#' Set up the parameters for the fixed point method. 
#' 
#' The 2nd and 3rd input are the regularizers for the L2 penalty of CNV before and after WGD, and (haplotype) ploidy from their mean.
#' The 4th input r_discontinuity is a penalty for the difference between two adjacent loci. If this is set to 0, then identity matrix would be used as matrix L in the model. 
#' 
#' The following parameters indicate the types of recurrent CNV that should be included into the model.
#' The default setting is the changes before WGD affecting all samples and changes after WGD affecting only samples with WGD.
#' However, for changes after WGD, the mean value is set to 0 for major copy, but -1 to minor copy, because it is common to see 1 copy number loss after WGD on minor copy.
#' A_m,B_m and Mean_m should have the same lenth and so do the ones for major copy. And they should all be smaller or equal to the length of r_CNV.
#' 
#' This function will return all these parameters in a list for fixed-ponit method.
#' 
#' 
#' @param data The CNV dataset (see \code{\link[graphics]{load_data}})
#' 
#' @param r_CNV Regularizer penalty for all recurrent CNV
#' 
#' @param r_ploidy Regularizer penalty for (haplotype) ploiidy
#' 
#' @param r_discontinuity Regularizer penalty for the discotinuity between two adjacent loci
#'
#' @param par_a_m,par_a_p The constant indicating how ploidy infecting the final copy number profile for minor and major copy respectively
#' 
#' @param A_m,A_p A vector containing the coeffient of the effect of the recurrent CNV with respect to the ploidy value for minor and major copy respectively
#' 
#' @param B_m,B_p A vector containing the coeffient of the effect of the recurrent CNV with respect to the constant value for minor and major copy respectively
#'
#' @param Mean_m,Mean_p A vector containing mean values of all the recurrent CNV for minor and major copy respectively
#'
#' @return List of parameters used for analysis
#'
#'
#' @examples
#' 
#' CNV_par<-set_par(wkdata)
#'
#' @export
#' 
#'


set_par<-function(data,r_CNV=c(10,10),r_ploidy=5,r_discontinuity=0.5,
                  par_a_m=1,par_a_p=1,
                  A_m=c(1,1),B_m=c(0,-1),Mean_m=c(0,-1),
                  A_p=c(1,1),B_p=c(0,-1),Mean_p=c(0,0)){
  weight_L=diag(data$nloci)
  D<-matrix(0*(data$nloci-1)*data$nloci,nrow=data$nloci-1,ncol=data$nloci)
  for(i in 1:(data$nloci-1)){
    if(data$loci[i,1]==data$loci[i+1,1]){
      D[i,i]=(3*10^9/data$nloci)/(data$loci[i+1,2]-data$loci[i,2])
      D[i,i+1]=(3*10^9/data$nloci)/(data$loci[i,2]-data$loci[i+1,2])
    }
  }
  D<-t(D)%*%D
  penal_discontinuity=r_discontinuity
  weight_L=weight_L+penal_discontinuity*D
  
  param<-list(  par_a_m=par_a_m,
                par_A_m=matrix(A_m,ncol=1),
                par_B_m=matrix(B_m,ncol=1),
                mean_m=matrix(Mean_m,ncol=1),
                penal_m=r_CNV[1:length(A_m)],
                
                par_a_p=par_a_p,
                par_A_p=matrix(A_p,ncol=1),
                par_B_p=matrix(B_p,ncol=1),
                mean_p=matrix(Mean_p,ncol=1),
                penal_p=r_CNV[1:length(A_p)],
                
                penal_g=r_ploidy,
                weight_L=weight_L
              )
  return(param)
}