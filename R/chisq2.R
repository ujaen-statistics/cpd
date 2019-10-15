#' @name chisq2.test
#' @title Pearson's Chi-squared Test for Count Data
#' @description bla bla bla
#' @param obs blaba
#' @param p.esp bla.
#' @param npar blabla
#' @param grouping blabla
#' @return test balb alb al a.
#' @importFrom stats pchisq
#' @export
chisq2.test <- function (obs, p.esp, npar = NULL, grouping = TRUE) {
  
  if (length(obs) != length(p.esp)) 
    stop("'obs' and 'p.esp' must have the same length")
  if (any(obs < 0) || anyNA(obs)) 
    stop("all entries of 'obs' must be nonnegative and finite")
  if ((n <- sum(obs)) == 0) 
    stop("at least one entry of 'obs' must be positive")
  if (sum(p.esp) < 0.9999999999999999) 
    stop("expected probabilities must sum 1")
  
  n <- sum(obs)
  m <- length(obs)
  
  esp <- n * p.esp
  
  nesp <- c()
  nobs <- c()
  i <- 1
  j <- 1
  suma.e <- 0
  suma.o <- 0
  
  if(grouping == TRUE){
    while(i <= length(esp)){
      suma.e <- suma.e + esp[i]
      suma.o <- suma.o + obs[i]
      if(suma.e >= 5){
        nesp[j] <- suma.e
        nobs[j] <- suma.o
        j <- j+1
        suma.e <- 0
        suma.o <- 0
      }
      i <- i+1
    }
    if (suma.e > 0){
      nesp[j-1] <- nesp[j-1] + suma.e
      nobs[j-1] <- nobs[j-1] + suma.o
    }
  }
  else{
    nesp <- esp
    nobs <- obs
  }
  
  if(is.null(npar)){
    DF <- length(nesp) - 1
  }
  else DF <- length(nesp) - npar - 1
  if(DF <= 0) DF <- NaN
  
  STATISTIC <- sum((nobs - nesp) ^ 2 / nesp)
  
  PVALUE <- pchisq(STATISTIC, DF, lower.tail = FALSE)
  
  ans<-list(statistic = STATISTIC, df = DF, p.value = PVALUE, observed = obs, observed.grouped = nobs,
                 expected = esp, expected.grouped = nesp, residuals = sqrt(STATISTIC))
  class(ans)<-"htest"
  ans
}