#' Load Dataset Containing Copy Number Profile
#'
#' Read text files containing copy number profiles called for multiple samples.
#' 
#' There are three different csv files containing minor copy, major copy and coordinates.
#'
#' This returns list containing all the input data we need for down-stream analysis.
#' The list consists of two values, which indicate the number of samples and loci respectively,
#' and three data frame which contain minor copy, major copy, and coordinates.
#'
#' @param mainDir The directory which contains all the files.
#' 
#' @param cancer_type The type of the cancer: the files are named as "cancer_type"_minorcn.csv or "cancer_type"_majorcn.csv
#'
#' @param coor The name for the coordinate file. Default as coords.txt.
#' 
#' @return A list containing: the number of samples, number of loci, all minor copies, all major copies, all loci
#'
#' @examples
#' 
#' wkdata<-load_data("C:\\Double","BLCA")
#'
#' @export
#'


load_data<-function(mainDir,cancer_type,coor="coords.txt"){
  setwd(mainDir)
  minor_file=paste(cancer_type,"_minorcn.csv",sep="")
  major_file=paste(cancer_type,"_majorcn.csv",sep="")
  
  coordinate_file=coor
  C_m=read.table(minor_file,sep=",",header=F)
  C_p=read.table(major_file,sep=",",header=F)
  coordinate_m=read.table(coordinate_file,sep=",",header=F)
  K=ncol(C_m)
  N=nrow(C_m)
  data<-list(nsample=K,nloci=N,minor=C_m,major=C_p,loci=coordinate_m)
}
