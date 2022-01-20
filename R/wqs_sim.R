#' WQS simulated dataset generator
#' 
#' \code{wqs_sim} generates a simulated dataset of mixture components, covariates, and outcomes 
#' based on an initial set of specifications. 
#'
#' @param nmix Number of mixture components in simulated dataset.
#' @param ncovrt Number of covariates in simulated dataset.
#' @param nobs Number of observations in simulated dataset.
#' @param ntruewts Number of mixture components that have a non-zero association with
#' the outcome (i.e. are not noise). 
#' @param ntruecovrt Number of covariates that have a non-zero association with the 
#' outcome (i.e. are not noise).
#' @param corrstruct Correlation matrix.
#' @param eps Error term.
#' @param truewqsbeta Simulated WQS beta_1 value. If NULL, then this value will be
#' randomly sampled. 
#' @param truebeta0 Simulated beta_0 value. If NULL, then this value will be
#' randomly sampled. 
#' @param truewts Simulated vector of mixture weights. If NULL, then this value will be
#' randomly sampled. 
#' @param truegamma Simulated gamma vector. If NULL, then this value will be
#' randomly sampled. 
#' @param rnd_wqsbeta_dir Direction of randomly sampled truewqsbeta (if truewqsbeta = NULL). 
#' You can choose between "positive", "negative", or NULL. If "positive" or "negative", 
#' the truewqsbeta will be sampled from a half normal distribution in either of those respective
#' directions. If NULL, then truewqsbeta will be sampled from a normal distribution. 
#' @param seed Random seed.
#' @param q Number of quantiles. 
#'
#' @return \code{wqs_perm} returns a list of:
#' \item{weights}{Simulated weights.}
#' \item{coef}{Simulated beta coefficients.}
#' \item{Data}{Simulated dataset.}
#' \item{yhat}{simulated predicted y values from the data generating model.}
#' \item{wqs}{Weighted quantile sum vector (quantile-transformed mixture components multiplied by weights).}
#' \item{modmat}{Model matrix.}
#' \item{Xq}{Quantile-transformed mixture components.}
#' 
#' @import mvtnorm extraDistr
#' @export wqs_sim
#'
#' @examples
wqs_sim <- function(nmix = 10, ncovrt = 10, nobs = 500, ntruewts = 10, ntruecovrt = 5, 
                    corrstruct = 0, eps = 1, truewqsbeta = NULL, truebeta0 = NULL, 
                    truewts = NULL, truegamma = NULL, rnd_wqsbeta_dir = "none", seed = 101,
                    q = 10) {
  
  if (length(corrstruct) == 1) {
    Rho <- diag(nmix + ncovrt)
    Rho[upper.tri(Rho)] <- Rho[lower.tri(Rho)] <- corrstruct
  } else {
    Rho <- corrstruct
  }
  
  weights <- rep(0, nmix)
  if (is.null(truewts)) {
    set.seed(seed)
    truewts <- extraDistr::rdirichlet(1, rep(1, ntruewts))
    weights[1:ntruewts] <- truewts
  } else {
    if (length(truewts) == nmix & sum(abs(truewts)) != 1) {
      truewts <- truewts / sum(truewts)
    }
    if (length(truewts) < nmix) {
      weights[1:length(truewts)] <- truewts
      weights[(length(truewts) + 1):nmix] <-
        (1 - sum(truewts)) / (nmix - length(truewts))
    } else {
      weights[1:length(truewts)] <- truewts
    }
  }
  
  if (round(sum(weights), 3) != 1.0) {
    warning(print(paste0("weights add up to ", sum(weights))))
  }
  
  set.seed(seed)
  Xmat <- rmvnorm(nobs, mean = rep(0, nmix + ncovrt), sigma = Rho)
  if (is.null(q)) {
    Xmatquant <- Xmat
  } else {
    Xmatquant <- Xmat
    Xmatquant[, 1:nmix] <- apply(
      Xmatquant[, 1:nmix],
      2,
      FUN = function(x) {
        as.numeric(as.character(cut(
          x,
          breaks = quantile(x, probs = seq(0, 1, by = (1 / q))),
          include.lowest = T,
          labels = 0:(q - 1)
        )))
      }
    )
  }
  
  if (ncovrt < ntruecovrt) {
    ntruecovrt <- ncovrt
  }
  if (!is.null(truegamma)) {
    if (length(truegamma) == 1) {
      covrtbetas <- rep(truegamma, ncovrt)
    } else {
      covrtbetas <- truegamma
    }
  } else {
    set.seed(seed)
    covrtbetas <- c(rnorm(ntruecovrt), rep(0, length = ncovrt - ntruecovrt))
  }
  
  set.seed(seed)
  if (!is.null(truebeta0)) {
    beta0 <- truebeta0
  } else {
    beta0 <- rnorm(1)
  }
  
  if (!is.null(truewqsbeta)) {
    wqsbeta <- truewqsbeta
  } else {
    set.seed(seed)
    if (rnd_wqsbeta_dir == "positive") {
      wqsbeta <- extraDistr::rhnorm(1)
    } else if (rnd_wqsbeta_dir == "negative") {
      wqsbeta <- extraDistr::rhnorm(1) * -1
    } else {
      wqsbeta <- rnorm(1)
    }
  }
  
  wqs <- Xmatquant[, 1:nmix] %*% weights
  if (ncovrt > 0) {
    modmat <- cbind(1, wqs, Xmat[, c((nmix + 1):(nmix + ncovrt))])
    dimnames(modmat)[[2]] <-
      c("Intercept", "wqs", paste0("C", 1:ncovrt))
    betas <- c(beta0, wqsbeta, covrtbetas)
    names(betas) <- c("beta0", "beta1", paste0("gamma", 1:ncovrt))
  } else {
    modmat <- cbind(1, wqs)
    dimnames(modmat)[[2]] <-
      dimnames(modmatq)[[2]] <- c("Intercept", "wqs")
    betas <- c(beta0, wqsbeta)
    names(betas) <- c("beta0", "beta1")
  }
  
  yhat <- modmat %*% betas
  set.seed(seed)
  epsilon <- rnorm(nobs, sd = eps)
  y <- yhat + epsilon
  Data <- data.frame(cbind(y, Xmat))
  
  if (ncovrt > 0) {
    names(Data) <- c("y", paste0("T", 1:nmix), paste0("C", 1:ncovrt))
    colnames(Xmatquant) <-
      c(paste0("T", 1:nmix), paste0("C", 1:ncovrt))
    colnames(Xmat) <- c(paste0("T", 1:nmix), paste0("C", 1:ncovrt))
  } else {
    names(Data) <- c("y", paste0("T", 1:nmix))
    colnames(Xmatquant) <- c(paste0("T", 1:nmix))
    colnames(Xmat) <- c(paste0("T", 1:nmix))
  }
  
  wtmat <- data.frame(mix_name = paste0("T", 1:nmix), true_weight = weights)
  
  retlist <-
    list(
      weights = wtmat,
      coef = betas,
      Data = Data,
      yhat = yhat,
      wqs = wqs,
      modmat = modmat,
      Xq = Xmatquant
    )
  return(retlist)
}