#' Full wrapper WQS permutation test 
#' 
#' \code{wqs_full_perm} is a full wrapper function that is a full implementation 
#' of the Weighted Quantile Sum (WQS) regression method followed by the 
#' permutation test to determine the significance of the WQS coefficient. 
#' 
#' @param formula An object of class formula. The wqs term must be included in 
#' the formula (e.g. y ~ wqs + ...).
#' @param data The \code{data.frame} to be used in the WQS regression run. 
#' @param mix_name A vector with the mixture column names. 
#' @param q An integer to indicate the number of quantiles to split the mixture 
#' variables. 
#' @param b_main The number of bootstraps for the main WQS regression run. 
#' @param b_perm The number of bootstraps for the permutation test WQS 
#' regression runs. 
#' @param b1_pos A logical value that indicates whether beta values should be 
#' positive or negative.
#' @param rs A logical value indicating whether random subset implementation 
#' should be performed. 
#' @param niter Number of permutation test iterations. 
#' @param seed An integer to fix the seed. This will only impact the the initial 
#' WQS regression run and not the permutation test iterations. 
#' @param plan_strategy Evaluation strategy for the plan function. You can choose 
#' among "sequential", "transparent", "multisession", "multicore", "multiprocess", 
#' "cluster" and "remote." See gWQS documentation for full details. 
#' @param ... Other parameters to put into the gwqs function call.
#' @param b1_constr Logical value that determines whether to apply positive or 
#' negative constraints in the optimization function for the weight optimization.
#' @param family A character value: `gaussian` for linear regression, `binomial` 
#' for logistic regression.
#' @param stop_if_nonsig if TRUE, the function will not proceed with the 
#' permutation test if the main WQS regression run produces nonsignificant 
#' p-value.
#' @param stop_thresh numeric p-value threshold required in order to proceed 
#' with the permutation test, if `stop_if_nonsig = TRUE`.
#'
#' @return \code{wqs_full_perm} returns an object of class `wqs_perm`, which 
#' contains three sublists: 
#' 
#' \item{perm_test}{List containing: (1) `pval`: permutation test p-value, (2) (linear 
#' regression only) `testbeta1`: reference WQS regression coefficient beta1 value, (3) 
#' (linear regression only) `betas`: Vector of beta values from each 
#' permutation test run, (4) (logistic regression only) `testpval`: test reference 
#' p-value, (5) (logistic regression only) `permpvals`: p-values from the null 
#' models.}
#' \item{gwqs_main}{Main gWQS object (same as model input).}
#' \item{gwqs_perm}{Permutation test reference gWQS object (NULL if model 
#' `family = "binomial"` or if same number of bootstraps are used in permutation 
#' test WQS regression runs as in the main run).}
#' @import gWQS
#' @export wqs_full_perm
#'
#' @examples
#' library(gWQS)
#'
#' # mixture names
#' PCBs <- names(wqs_data)[1:34]
#' 
#'  perm_test_res <- wqs_full_perm(formula = yLBX ~ wqs, data = wqs_data, 
#'                                 mix_name = PCBs, q = 10, b_main = 10, 
#'                                 b_perm = 5, b1_pos = TRUE, b1_constr = FALSE, 
#'                                 niter = 10, seed = 16, plan_strategy = "multicore", 
#'                                 stop_if_nonsig = FALSE)
#' 
#'  # Note: The default values of b_main = 1000, b_perm = 200, and niter = 200 
#'  # are the recommended parameter values. This example has a lower b_main, 
#'  # b_perm, and niter in order to serve as a shorter test run. 
#'  
wqs_full_perm <- function(formula, data, mix_name, q = 10, b_main = 1000, 
                          b_perm = 200, b1_pos = TRUE, b1_constr = FALSE, 
                          rs = FALSE, niter = 200, seed = NULL, 
                          family = "gaussian", plan_strategy = "multicore",
                          stop_if_nonsig = FALSE, stop_thresh = 0.05, ...){
  
  # run main WQS regression
  gwqs_res_main <- gWQS::gwqs(formula = formula, data = data, mix_name = mix_name, 
                              q = q, b = b_main, b1_pos = b1_pos, 
                              b1_constr = b1_constr, rs = rs, seed = seed, 
                              validation = 0, family = family, 
                              plan_strategy = plan_strategy, ...) 
  
  naive_p <- summary(gwqs_res_main)$coefficients["wqs", 4]

  if (stop_if_nonsig == TRUE & naive_p > stop_thresh){
    message(sprintf("The main WQS regression run did not give a significant 
                    result (p = %s)", 
                    naive_p))
    
    results <- list(gwqs_main = gwqs_res_main, 
                    family = gwqs_res_main$family$family,
                    gwqs_perm = NULL, 
                    perm_test = NULL)
  }
  else{
    # run permutation test (using wqs_perm function) 
    results <- wqs_perm(gwqs_res_main, niter = niter, boots = b_perm, 
                        b1_pos = b1_pos, b1_constr = b1_constr, rs = rs, 
                        plan_strategy = plan_strategy, seed = seed)
    
  }
  
  class(results) <- "wqs_perm"
  
  results
}