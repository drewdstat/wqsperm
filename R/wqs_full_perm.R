#' Full wrapper WQS permutation test 
#' 
#' \code{wqsfullperm} is a full wrapper function that is a full implementation of 
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
#' @param seed An integer to fix the seed.
#' @param plan_strategy (Taken from gWQS documentation) A character value that allows to 
#' choose the evaluation strategies for the plan function. You can choose among "sequential",
#' "transparent", "multisession", "multicore", "multiprocess", "cluster" and "remote."
#' @param ... 
#'
#' @return \code{wqsperm} returns three sublists: 
#' 
#' \item{gwqs_res_main}{Full output of all results from the main \code{wqs} run.}
#' \item{gwqs_res_perm}{Full output of all results from the permutation test \code{wqs} run.}
#' \item{perm_test_res}{Results from the permutation test, including \code{pval, 
#' testbeta1} and \code{betas}. For more information, see the outputs for the \code{wqsperm}
#' function.}
#' @import gWQS
#' @export wqs_full_perm
#'
#' @examples
wqs_full_perm <- function(formula, data, mix_name, q = 4, b_main = 1000, b_perm = 200,
                          b1_pos = TRUE, rs = FALSE, niter = 200, seed = NULL, 
                          plan_strategy = "multicore", ...){
  
  # run main WQS 
  gwqs_res_main <- gWQS::gwqs(formula = formula, data = data, mix_name = mix_name, q = q, 
                              b = b_main, b1_pos = b1_pos, rs = rs, seed = seed, validation = 0,
                              family = "gaussian", plan_strategy = plan_strategy, ...) 
  
  # run permutation test (using wqsperm function) 
  results <- wqs_perm(gwqs_res_main, niter = niter, boots = b_perm, b1_pos = b1_pos, 
                      rs = rs, plan_strategy = plan_strategy, seed = seed)
  
  class(results) <- "wqs_perm"
  
  results
}