#' WQS permutation test
#' 
#' \code{wqsperm} takes a gwqs object as an input and runs the permutation test (Day 
#' et al, 2021) to obtain an estimate for the p-value significance for the WQS coefficient.  
#' 
#' Note that to use this function, there are some restrictions that users should be aware of. 
#' For complete details, please reference the vignette.
#'
#' @param model A \code{gwqs} object as generated from the \code{gWQS} package.  
#' @param niter Number of permutation test iterations. 
#' @param boots Number of bootstrap samples for each permutation test \code{wqs} run.  
#' @param b1_pos A logical value that indicates whether beta values should be positive 
#' or negative.
#' @param rs A logical value indicating whether random subset implementation should be 
#' performed. 
#' @param plan_strategy (Taken from gWQS documentation) A character value that allows to 
#' choose the evaluation strategies for the plan function. You can choose among "sequential",
#' "transparent", "multisession", "multicore", "multiprocess", "cluster" and "remote."
#' @param seed Random seed for the permutation test WQS reference run.  
#'
#' @return \code{wqsperm} returns three objects: 
#' 
#' \item{pval}{The p-value obtained from the permutation test.}
#' \item{testbeta1}{Reference WQS coefficient value.}
#' \item{betas}{Vector of beta values from each permutation test run.}
#' @import gWQS
#' @export wqsperm
#'
#' @examples
wqsperm <- function(model, niter = 200, boots = 200, b1_pos = TRUE, rs = FALSE, 
                    plan_strategy = "multicore", seed = NULL) {
  
  if (class(model) == "gwqs") {
    if (model$family$family != "gaussian" | model$family$link != "identity"){
      stop("The permutation test is currently only set up to accomodate the Gaussian family with an 
           identity link.")
    }
  } else stop("'model' must be of class 'gwqs' (see gWQS package).")
  
  mm <- model$fit
  formchar <- as.character(formula(mm))
  
  if (!is.null(model$stratified) | grepl("wqs:", formchar[3], fixed = T))
  {
    # TODO: We should be able to accomodate stratified weights though we haven't tested that yet,
    # and I'm not sure it makes sense to have stratified weights without a WQS interaction term.
    stop("This permutation test is not yet set up to accomodate stratified weights or 
         WQS interaction terms.")
  }  
  
  cl = match.call()
  Data <- model$data[model$vindex, -which(names(model$data) %in% c("wqs", "wghts"))]
  yname <- as.character(formula(mm))[2]
  mix_name <- names(model$bres)[names(model$bres) %in% model$final_weights$mix_name]

  if (!is.null(model$qi)) {
    nq <- max(sapply(model$qi, length)) - 1
  } else {
    # this is for cases when there is no quantile transformation or it's already been
    # done in the data frame
    nq <- NULL
  }

  # reference WQS run 
  if (is.null(boots)){
    boots <- length(model$bindex)
  }
  
  if (boots == length(model$bindex)){
    perm_ref_wqs <- model
    ref_beta1 <- mm$coef[2]
  }
  
  else{
    perm_ref_wqs <- gwqs(formula = formula(mm), data = Data, mix_name = mix_name, 
                         q = nq, b = boots, rs = rs, validation = 0, plan_strategy = plan_strategy,
                         b1_pos = b1_pos, seed = seed)
    
    ref_beta1 <- perm_ref_wqs$fit$coef[2]
  }
  
  
  if (length(mm$coef) > 2) {
    # This is the permutation test algorithm when there are multiple independent variables in
    # the model
    lm_form <- formula(paste0(formchar[2], formchar[1], gsub("wqs + ", "", formchar[3], fixed = T)))
    fit.partial <- lm(lm_form, data = Data)
    partial.yhat <- predict(fit.partial)
    partial.resid <- resid(fit.partial)
    reorgmat <- matrix(NA, dim(Data)[1], niter)
    reorgmat <- apply(reorgmat, 2, function(x) partial.yhat + sample(partial.resid, replace = F))
  } else {
    # This is the permutation test algorithm when there is only one independent variable in
    # the model
    reorgmat <- matrix(NA, dim(Data)[1], niter)
    reorgmat <- apply(reorgmat, 2, function(x) sample(Data[, yname]))
  }
  
  getbetas <- function(x) {
    
    newDat <- Data
    newDat[, yname] <- x
    names(newDat) <- c(names(Data))
    formchar <- as.character(formula(mm))
    
    if (length(mm$coef) > 2) {
      form1 <- formula(paste0(formchar[2], formchar[1], formchar[3]))
    } else {
      form1 <- formula(paste0(formchar[2], formchar[1], "wqs"))
    }
    
    gwqs1 <- tryCatch({
      suppressWarnings(gwqs(formula = form1, data = newDat, mix_name = mix_name, 
                            q = nq, b = boots, rs = rs, validation = 0, plan_strategy = plan_strategy,
                            b1_pos = b1_pos))
    }, error = function(e) NULL, 
      warning = function(e) ifelse(rs == TRUE, message("WQSRS failed"), message("WQS failed")))
    
    if (is.null(gwqs1))
      lm1 <- NULL else lm1 <- gwqs1$fit
    if (is.null(lm1)) {
      retvec <- NA
    } else {
      retvec <- lm1$coef[2]
    }
    return(retvec)
  }
  
  pbapply::pboptions(type = "timer")
  betas <- pbapply::pbapply(reorgmat, 2, getbetas)
  
  if (any(is.na(betas))) {
    print(paste0(length(which(is.na(betas))), " failed model attempts"))
  }
  
  calculate_pval <- function(x, true, posb1 = b1_pos) {
    if (posb1) {
      length(which(x > true))/length(betas)
    } else {
      length(which(x < true))/length(betas)
    }
  }
  
  pval <- calculate_pval(betas, ref_beta1, b1_pos)
  
  perm_retlist <- list(pval = pval, testbeta1 = ref_beta1, betas = betas, call = cl)
    
  ret_ref_wqs <- ifelse(boots == length(model$bindex), NULL, perm_ref_wqs)
  
  results <- list(gwqs_main = model, 
                  gwqs_perm = ret_ref_wqs, 
                  perm_test = perm_retlist)
  
  class(results) <- "wqsperm"
  
  results
}

#' @rawNamespace S3method(print, wqsperm)
#' @rdname methods

print.wqsperm <- function(x, ...){
  
  cat("Permutation test WQS coefficient p-value: \n", 
      x$perm_test$pval,
      "\n")

  main_sum <- summary(x$gwqs_main)
  final_weights <- x$gwqs_main$final_weights
  
  print(main_sum)

}

#' @rawNamespace S3method(summary, wqsperm)
#' @rdname methods

summary.wqsperm <- function(x, ...){
  
  cat("Permutation test WQS coefficient p-value: \n", 
      x$perm_test$pval,
      "\n")
  
  main_sum <- summary(x$gwqs_main)
  final_weights <- x$gwqs_main$final_weights
  
  print(main_sum)

}

