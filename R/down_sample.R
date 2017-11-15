#' Down Sample The Original Dataset
#' 
#' Down-sample the dataset in order to run faster.
#' 
#' For the copy number dataset is usually quite large to work with, and the fixed-point iteration is quite computationally expensive,
#' to save some time, it is helpful to down-sample the number of loci (no change to the number of samples).
#' 
#' This function down-samples the original one to make the algorithms run faster, and because most loci have similar values to the adjacent one,
#' this doesn't cause much inaccuracy.
#'
#' @param CNV_data The CNV dataset (see \code{\link[graphics]{load_data}})
#' 
#' @param factor The factore to which the dataset to do down sampled
#'
#' 
#' @return Dataset with the minor and major copy, as well as the number of loci being down sampled (see \code{\link[graphics]{load_data}})
#'
#' @seealso \code{\link{load_data}}
#'
#' @examples
#' 
#' wkdata<-down_sample(wkdata,10)
#'
#' @export
#'


down_sample<-function(CNV_data,factor=10){
  C_m<-as.matrix(CNV_data$minor[seq(1,nrow(CNV_data$minor),factor),])
  C_p<-as.matrix(CNV_data$major[seq(1,nrow(CNV_data$major),factor),])
  coordinate_m<-as.matrix(CNV_data$loci[seq(1,nrow(CNV_data$loci),factor),])
  K=ncol(C_m)
  N=nrow(C_m)
  data<-list(nsample=K,nloci=N,minor=C_m,major=C_p,loci=coordinate_m)
  return(data);
}