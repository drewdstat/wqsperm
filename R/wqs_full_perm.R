#' Full wrapper WQS permutation test 
#' 
#' \code{wqs_full_perm} is a full wrapper function that is a full implementation of 
#' the Weighted Quantile Sum (WQS) method followed by the permutation test to determine 
#' the significance of the WQS coefficient. 
#' 
#' @param formula An object of class formula. The wqs term must be included in the
#' formula (e.g. y ~ wqs + ...).
#' @param data The \code{data.frame} to be used in the WQS run. 
#' @param mix_name A vector with the mixture column names. 
#' @param q An integer to indicate the number of quantiles to split the mixture variables. 
#' @param b_main The number of bootstraps for the main WQS run. 
#' @param b_perm The number of bootstraps for the permutation test WQS runs. 
#' @param b1_pos A logical value that indicates whether beta values should be positive 
#' or negative.
#' @param rs A logical value indicating whether random subset implementation should be 
#' performed. 
#' @param niter Number of permutation test iterations. 
#' @param seed An integer to fix the seed. This will only impact the the initial WQS run and
#' not the permutation test iterations. 
#' @param plan_strategy Evaluation strategy for the plan function. You can choose among "sequential",
#' "transparent", "multisession", "multicore", "multiprocess", "cluster" and "remote." See gWQS documentation
#' for full details. 
#' @param ... Other parameters to put into the gwqs function call
#'
#' @return \code{wqs_full_perm} returns an object of class `wqs_perm`, which contains 
#' three sublists: 
#' 
#' \item{perm_test}{Contains three objects: (1) `pval`: p-value for the proportion of 
#' permuted WQS coefficient values greater than the reference value, (2) `testbeta1`: reference WQS coefficient beta1 value, 
#' (3) `betas`: Vector of beta values from each permutation test run.}
#' \item{gwqs_main}{main gWQS object (same as model input)}
#' \item{gwqs_perm}{permutation test reference gWQS object (NULL if same number of bootstraps
#' as main gWQS object)}
#' @import gWQS
#' @export wqs_full_perm
#'
#' @examples
#'  perm_test_res <- wqs_full_perm(formula = yLBX ~ wqs, data = wqs_data, mix_name = PCBs,
#'                                q = 10, b_main = 10, b_perm = 5, b1_pos = T, niter = 10,
#'                                 seed = 16, plan_strategy = "multicore", returnbetas = TRUE)
#' 
#'  # Note: The default values of b_main = 1000, b_perm = 200, and niter = 200 are the recommended
#'  # parameter values. This example has a lower b_main, b_perm, and niter in order to serve 
#'  # as a shorter test run. 
#'  
wqs_full_perm <- function(formula, data, mix_name, q = 10, b_main = 1000, b_perm = 200,
                          b1_pos = TRUE, rs = FALSE, niter = 200, seed = NULL, 
                          plan_strategy = "multicore", ...){
  
  # run main WQS 
  gwqs_res_main <- gWQS::gwqs(formula = formula, data = data, mix_name = mix_name, q = q, 
                              b = b_main, b1_pos = b1_pos, rs = rs, seed = seed, validation = 0,
                              family = "gaussian", plan_strategy = plan_strategy, ...) 
  
  # run permutation test (using wqs_perm function) 
  results <- wqs_perm(gwqs_res_main, niter = niter, boots = b_perm, b1_pos = b1_pos, 
                      rs = rs, plan_strategy = plan_strategy, seed = seed)
  
  class(results) <- "wqs_perm"
  
  results
}