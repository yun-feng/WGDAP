context("Fixed Point Iterations")


test_that("down_sampling works correctly", {
  
  wkdata<-down_sample(wkdata)
  expect_equal( wkdata$nloci, 10 )
  expect_equal( wkdata$nsample, 20 )
  
})

test_that("Fixed Point Iterations", {
  
  CNV_par<-set_par(wkdata)
  common_var<-init_var(wkdata,CNV_par)
  init_loss<-feval(wkdata,CNV_par,common_var)
  FPI_result=FPI(wkdata,CNV_par,common_var,init_loss,max_iter=4)
  expect_equal( FPI_result$loss_array[1], init_loss )
  expect_lt( FPI_result$loss_array[5], init_loss )
})
