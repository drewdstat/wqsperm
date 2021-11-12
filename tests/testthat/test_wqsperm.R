context("Testing the wqsperm function")

test_that("Function fails with non-Gaussian input", {
  
  library(gWQS)
  
  PCBs <- names(wqs_data)[1:34]
  
  wqs_perm <- gWQS::gwqs(ybin ~ wqs, mix_name = PCBs, data = wqs_data, q = 10, validation = 0, 
                         b = 5, b1_pos = T, bl_constr = F, family = "binomial", seed = 16)
  
  expect_error(wqsperm(wqs_perm, niter = 10, boots = 5, b1_pos = T))
})
