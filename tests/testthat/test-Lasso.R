context("Post-processing")


test_that("lasso solvers works correctly", {
  
  wkdata<-data_lasso(wkdata,g)
  var<-list( inf_g=g,
             inf_X_m=matrix(rep(c(0,-1),each=N),nrow=N,ncol=2),
             inf_X_p=matrix(rep(c(0,0),each=N),nrow=N,ncol=2) )
  var<-var_lasso(var,1)
  par<-par_lasso(N,K,g)
  Lasso_res<-Run_lasso(wkdata,par,var,g,max=5)
  expect_lt(Lasso_res$loss_array[6],Lasso_res$loss_array[1])
})
