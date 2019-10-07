#' Maximum-likelihood fitting of the Complex Triparametric Pearson (CTP) distribution
#'
#' @description
#' Maximum-likelihood fitting of the Complex Triparametric Pearson (CTP) distribution with parameters \eqn{a}, \eqn{b} and \eqn{\gamma}. 
#'
#' @usage
#' fitctp(x, astart = 0, bstart = 1, gammastart = 1.1, method = "L-BFGS-B", 
#'        moments = FALSE, hessian = TRUE, control = list(), ...)
#'
#' @param x A numeric vector of length at least one containing only finite values.
#' @param astart An starting value for the parameter \eqn{a}; by default 0.
#' @param bstart An starting value for the parameter \eqn{b}; by default 1.
#' @param gammastart An starting value for the parameter \eqn{\gamma}; by default 1.1.
#' @param method The method to be used in fitting the model. The default method is "L-BFGS-B" (optim).
#' @param moments If \code{TRUE} the estimates of \eqn{a}, \eqn{b} and \eqn{\gamma} by the method of moments are used as starting values (if it is posible). By default this argument is \code{FALSE}.
#' @param hessian If \code{TRUE} the hessian of the objective function at the minimum is returned.
#' @param control A list of parameters for controlling the fitting process.
#' @param ...  Additional parameters.
#'
#' @return An object of class "fitctp" is a list containing the following components:
#'
#' \itemize{
#' \item \code{n}, the number of observations,
#' \item \code{initialValues}, a vector with the starting values used,
#' \item \code{coefficients}, the parameter ML estimates of the CTP distribution,
#' \item \code{se}, a vector of the standard error estimates,
#' \item \code{hessian}, a symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function,
#' \item \code{cov}, an estimate of the covariance matrix of the model coefficients,
#' \item \code{corr}, an estimate of the correlation matrix of the model estimates,
#' \item \code{loglik}, the maximized log-likelihood,
#' \item \code{aic}, Akaike Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{bic}, Bayesian Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{code}, a code that indicates successful convergence of the fitter function used (see nlm and optim helps),
#' \item \code{converged},  logical value that indicates if the optimization algorithms succesfull,
#' \item \code{method}, the name of the fitter function used.
#' }
#' 
#' @importFrom stats nlm optim coef runif
#' @export
#' @references 
#' 
#' \insertRef{RCSO2004}{cpd}
#' 
#' \insertRef{ROC2018}{cpd}
#' 
#' @seealso
#' Maximum-likelihood fitting for the CBP distribution: \code{\link{fitcbp}}.
#'
#' @examples
#' set.seed(123)
#' x <- rctp(500, -0.5, 1, 2)
#' fitctp(x)
#' fitctp(x, astart = 1, bstart = 1.1, gammastart = 3)
#' fitctp(x, moments = TRUE)

fitctp <- function(x, astart = 0, bstart = 1, gammastart = 1.1, method = "L-BFGS-B", moments=FALSE, hessian = TRUE, control = list(), ...)
{

  #Checking
  if ( mode(x) != "numeric")
    stop( paste ("non-numeric argument to mathematical function", sep = ""))

  if( is.null(control$maxit) )
    control$maxit = 10000

  if( is.null(control$trace) )
    control$trace = 0

  if ( !is.numeric(control$maxit) || control$maxit <= 0 )
    stop( "maximum number of iterations must be > 0" )

  if (moments){
    #Try to generate initial values using moments method
    n <- length(x)
    m1 <- mean(x)
    m2 <- sum(x ^ 2) / n
    m3 <- sum(x ^ 3) / n

    #System equations
    A <- matrix(c(1, 1 + m1, 1 + 2 * m1 + m2, m1, m1 + m2, m1 + 2 * m2 + m3, -m1, -m2, -m3), ncol = 3, nrow = 3)
    b <- matrix(c(0, -m2, -2 * m3 - m2), ncol = 1, nrow = 3)

    #Solve the system
    sol <- try(solve(A)%*%b,silent=TRUE)
    if ('try-error' %in% class(sol)){
      warning(paste("The moment method hasn't solution continue with the default initial values [",paste(astart,bstart,gammastart,sep=", "),"]"))
    }else{
      #New initial values
      astart <- sol[2] / 2
      bstart <- sqrt(sol[1] - astart ^ 2)
      gammastart <- sol[3] + 1
    }

  }
  #Imaginary unit
  i<-sqrt(as.complex(-1))

  #Log-likelihood
  if (method != "L-BFGS-B") {
    logL <- function(p){
            a <- p[1]
            b <- p[2]
            gama <- exp(p[3])
            respuesta <- -sum(2*Re(cgamma(gama-a+b*i,log=TRUE)+log(cgamma(a+b*i+x))-log(cgamma(a+b*i)))-lgamma(gama-2*a)-lgamma(gama+x))
            return(respuesta)
      }

      pstart <- c(astart,bstart,log(gammastart))

    }else {
      logL <- function(p){
        a <- p[1]
        b <- p[2]
        gama <- p[3]
        respuesta <- -sum(2*Re(cgamma(gama-a+b*i,log=TRUE)+log(cgamma(a+b*i+x))-log(cgamma(a+b*i)))-lgamma(gama-2*a)-lgamma(gama+x))
        return(respuesta)
      }

      pstart <- c(astart,bstart,gammastart)
    }

  #Optimization process

  if (method == "nlm") {
    fit <- nlm(logL, p = pstart, hessian = hessian, iterlim = control$maxit, print.level = control$trace)
    fit$value <- fit$minimum
    fit$par <- fit$estimate
    fit$convergence <- fit$code
    if (fit$convergence < 3)
      fit$converged = TRUE
    else fit$converged = FALSE
    methodText = "nlm"
  }
  else if (any(method == c("Nelder-Mead", "BFGS", "CG", "SANN"))) {
    fit <- optim(pstart, logL, method = method, hessian = hessian,
                 control = list(maxit = control$maxit, trace = control$trace))
    methodText <- method
    if (fit$convergence == 0)
      fit$converged = TRUE
    else fit$converged = FALSE
  }
  else if (any(method == c("L-BFGS-B"))) {
    fit<-optim(pstart, logL, method = method, lower = c(-Inf,0,0.0000001), upper = c(Inf,Inf,Inf), hessian = hessian, control = list(maxit = control$maxit, trace = control$trace))
    methodText <- method
    if (fit$convergence == 0)
      fit$converged = TRUE
    else fit$converged = FALSE
  }else{
    stop("Incorrect method")
  }

  if (method != "L-BFGS-B") {
    coef.table<-rbind(fit$par,deparse.level=0)
    dimnames(coef.table)<-list("",c("a","b","log(gamma)"))
  }
  else {
    coef.table<-rbind(fit$par,deparse.level=0)
    dimnames(coef.table)<-list("",c("a","b","gamma"))
  }

  #Results
  results<-list(
    x = x,
    n=length(x),
    loglik = -(fit$value+sum(lfactorial(x))),
    aic = 2*(fit$value+sum(lfactorial(x)))+6,
    bic = 2*(fit$value+sum(lfactorial(x)))+2*log(length(x)),
    coefficients =  coef.table,
    code = fit$convergence,
    hessian = fit$hessian,
    cov = solve(fit$hessian),
    se = sqrt(diag(solve(fit$hessian))),
    corr = solve(fit$hessian)/(sqrt(diag(solve(fit$hessian))) %o%
                                 sqrt(diag(solve(fit$hessian)))),
    code = fit$convergence,
    converged = fit$converged,
    initialValues=c(astart,bstart,gammastart),
    method = methodText
  )
  class(results) <- "fitCTP"
  return(results)
}


#' Maximum-likelihood fitting of the Complex Biparametric Pearson (CBP) distribution
#'
#' @description
#' Maximum-likelihood fitting of the Complex Biparametric Pearson (CBP) distribution with parameters \eqn{b} and \eqn{\gamma}. 
#'
#' @usage
#' fitcbp(x, bstart = 1, gammastart = 1.1, method = "L-BFGS-B", 
#'        moments = FALSE, hessian = TRUE, control = list(), ...)
#' 
#' @param x A numeric vector of length at least one containing only finite values.
#' @param bstart An starting value for the parameter \eqn{b}; by default 1.
#' @param gammastart An starting value for the parameter \eqn{\gamma}; by default 1.1.
#' @param method The method to be used in fitting the model. The default method is "L-BFGS-B" (optim).
#' @param moments If \code{TRUE} the estimates of \eqn{b} and \eqn{\gamma} by the method of moments are used as starting values (if it is posible). By default this argument is \code{FALSE}.
#' @param hessian If \code{TRUE} the hessian of the objective function at the minimum is returned.
#' @param control A list of parameters for controlling the fitting process.
#' @param ...  Additional parameters.
#'
#' @return An object of class "fitcbp" is a list containing the following components:
#'
#' \itemize{
#' \item \code{n}, the number of observations,
#' \item \code{initialValues}, a vector with the starting values used,
#' \item \code{coefficients}, the parameter ML estimates of the CTP distribution,
#' \item \code{se}, a vector of the standard error estimates,
#' \item \code{hessian}, a symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function,
#' \item \code{cov}, an estimate of the covariance matrix of the model coefficients,
#' \item \code{corr}, an estimate of the correlation matrix of the model estimates,
#' \item \code{loglik}, the maximized log-likelihood,
#' \item \code{aic}, Akaike Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{bic}, Bayesian Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters,
#' \item \code{code}, a code that indicates successful convergence of the fitter function used (see nlm and optim helps),
#' \item \code{converged},  logical value that indicates if the optimization algorithms succesfull,
#' \item \code{method}, the name of the fitter function used.
#' }
#' 
#' @importFrom stats nlm optim coef runif
#' @export
#' 
#' @references 
#' 
#' \insertRef{RCS2003}{cpd}
#' 
#' @seealso
#' Maximum-likelihood fitting for the CTP distribution: \code{\link{fitctp}}.
#'
#' @examples
#' set.seed(123)
#' x <- rcbp(500, 1.75, 3.5)
#' fitcbp(x)
#' fitcbp(x, bstart = 1.1, gammastart = 3)
#' fitcbp(x, moments = TRUE)

fitcbp <- function(x, bstart = 1, gammastart = 1.1, method = "L-BFGS-B", moments = FALSE, hessian = TRUE, control = list(), ...)
{

  #Checking
  if ( mode(x) != "numeric")
    stop( paste ("non-numeric argument to mathematical function", sep = ""))

  if( is.null(control$maxit) )
    control$maxit = 10000

  if( is.null(control$trace) )
    control$trace = 0

  if ( !is.numeric(control$maxit) || control$maxit <= 0 )
    stop( "maximum number of iterations must be > 0" )

  if (moments){
    #Try to generate initial values using moments method
    n <- length(x)
    m1 <- mean(x)
    m2 <- sum(x ^ 2) / n
    if (m2 - m1 ^ 2 -m1 <= 0){
      warning(paste("Data must be overdispersed. The moment method hasn't solution continue with the default initial values [",
                    paste(bstart,gammastart,sep=", "),"]"))
    }else{
      #New initial values
      bstart <- sqrt(m1 * m2 / (m2 - m1  ^ 2 - m1))
      gammastart <- bstart ^ 2 / m1 +1
    }
  }
  
  #Imaginary unit
    i<-sqrt(as.complex(-1))

  #Log-likelihood
    if (method != "L-BFGS-B") {
    logL <- function(p){
      b <- p[1]
      gama <- exp(p[2])
      respuesta <- -sum(2*Re(cgamma(gama+b*i,log=TRUE)+log(cgamma(b*i+x))-log(cgamma(b*i)))-lgamma(gama)-lgamma(gama+x))
      return(respuesta)
    }

    pstart <- c(bstart,log(gammastart))

  }else {
    logL <- function(p){
      b <- p[1]
      gama <- p[2]
      respuesta <- -sum(2*Re(cgamma(gama+b*i,log=TRUE)+log(cgamma(b*i+x))-log(cgamma(b*i)))-lgamma(gama)-lgamma(gama+x))
      return(respuesta)
    }

    pstart <- c(bstart,gammastart)
  }

  #Optimization process
  if (method == "nlm") {
    fit <- nlm(logL, p = pstart, hessian = hessian, iterlim = control$maxit, print.level = control$trace)
    fit$value <- fit$minimum
    fit$par <- fit$estimate
    fit$convergence <- fit$code
    if (fit$convergence < 3)
      fit$converged = TRUE
    else fit$converged = FALSE
    methodText = "nlm"
  }
  else if (any(method == c("Nelder-Mead", "BFGS", "CG", "SANN"))) {
    fit <- optim(pstart, logL, method = method, hessian = hessian,
                 control = list(maxit = control$maxit, trace = control$trace))
    methodText <- method
    if (fit$convergence == 0)
      fit$converged = TRUE
    else fit$converged = FALSE
  }
  else if (any(method == c("L-BFGS-B"))) {
    fit<-optim(pstart, logL, method = method, lower = c(0,0.0000001), upper = c(Inf,Inf), hessian = hessian, control = list(maxit = control$maxit, trace = control$trace))
    methodText <- method
    if (fit$convergence == 0)
      fit$converged = TRUE
    else fit$converged = FALSE
  }else{
    stop("Incorrect method")
  }

  if (method != "L-BFGS-B") {
    coef.table<-rbind(fit$par,deparse.level=0)
    dimnames(coef.table)<-list("  ",c("b","log(gamma)"))
  }
  else {
    coef.table<-rbind(fit$par,deparse.level=0)
    dimnames(coef.table)<-list("  ",c("b","gamma"))
    }

  #Results
  results<-list(
    x = x,
    n=length(x),
    loglik = -(fit$value+sum(lfactorial(x))),
    aic = 2*(fit$value+sum(lfactorial(x)))+4,
    bic = 2*(fit$value+sum(lfactorial(x)))+2*log(length(x)),
    coefficients =  coef.table,
    code = fit$convergence,
    hessian = fit$hessian,
    cov = solve(fit$hessian),
    se = sqrt(diag(solve(fit$hessian))),
    corr = solve(fit$hessian)/(sqrt(diag(solve(fit$hessian))) %o%
                                 sqrt(diag(solve(fit$hessian)))),
    code = fit$convergence,
    converged = fit$converged,
    initialValues=c(bstart,gammastart),
    method = methodText
  )
  class(results) <- "fitCBP"
  return(results)
}

#' @export
logLik.fitCTP <- function (object, ...){
  val <- object$loglik
  attr(val, "nobs") <- object$n
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  val
}

#' @method print fitCTP
#' @export
print.fitCTP<-function (x, digits = getOption("digits"), ...) {
  if (length(coef(x))) {
    ans=format(rbind(dimnames(x$coefficients)[[2L]],
                     format(x$coefficients,digits = digits),
              sapply(format(x$se,digits = digits),function(v){paste("(",v,")",sep="")})),justify = "centre")
    colnames(ans)<-ans[1,]
    ans<-ans[2:3,]
    print.default(ans, print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")

  if(!x$converged){
    cat("Error of convergence")
  }
  invisible(x)
}

#' @export
logLik.fitCBP <- function (object, ...){
  val <- object$loglik
  attr(val, "nobs") <- object$n
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  val
}

#' @method print fitCBP
#' @export
print.fitCBP<-function (x, digits = getOption("digits"), ...) {
  if (length(coef(x))) {
    ans=format(rbind(dimnames(x$coefficients)[[2L]],
                     format(x$coefficients,digits = digits),
                     sapply(format(x$se,digits = digits),function(v){paste("(",v,")",sep="")})),justify = "centre")
    colnames(ans)<-ans[1,]
    ans<-ans[2:3,]
    print.default(ans, print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  
  if(!x$converged){
    cat("Error of convergence")
  }
  invisible(x)
}


#' @method summary fitCBP
#' @importFrom stats pchisq pnorm stepfun
#' @importFrom dgof ks.test
#' @export
summary.fitCBP<-function (x, ...) {
  digits = getOption("digits")
  if ("fitCBP" %in% class(x)){
    myfun <- stepfun(0:max(x$x), c(0, pcbp(0:max(x$x), b = x$coefficients[1], gamma = x$coefficients[2])))
    x$kstest <- dgof::ks.test(x$x,myfun , simulate.p.value=TRUE, B=1000)
    
    xmax<-max(x$x)
    p.esp<-dcbp(0:(xmax-1),x$coefficients[1],x$coefficients[2])
    p.esp[xmax]<-1-sum(p.esp[0:(xmax-1)])
    
    coefName<-rbind("",cbind(dimnames(x$coefficients)[[2L]]))
    coefValu<-rbind("Estimate",cbind(as.vector(format(x$coefficients,digits = digits))))
    se<-rbind("Std. Error",cbind(as.vector(format(x$se,digits = digits))))
    x$zvalue<-x$coefficients/x$se
    zvalue<-rbind("z-value",cbind(as.vector(format(x$zvalue,digits = digits))))
    x$pvalue<- 2 * pnorm(abs(x$zvalue),lower.tail = FALSE)
    prz<-rbind("Pr(>|z|)",cbind(as.vector(format(x$pvalue,digits = digits))))
    ans<-cbind(format(cbind(coefValu,se,zvalue,prz),justify = "centre"))
      
    cat("Parameters:")
    prmatrix(ans, rowlab=coefName, collab=rep("",4),quote=FALSE)
    cat(paste("\nLoglikelihood: ",round(x$loglik,2),"   AIC: ",round(x$aic,2),"   BIC: ",round(x$bic,2),"\n",sep=""))
    cat("\nGoodness-of-fit test:\n")
    cat(paste("Chi-2: ",x$chi2," (",x$chi2p,")   Kolmogorov-Smirnov: ",format(x$kstest$statistic,digits=digits)," (p-value: ",x$kstest$p.value,")\n",sep=""))
    cat("\nCorrelation Matrix:\n")
    prmatrix(x$corr, rowlab=dimnames(x$coefficients)[[2L]], collab=dimnames(x$coefficients)[[2L]],quote=FALSE)
  }
  invisible(x)
}





#Pearson's Chi-squared Test for Count Data 
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
  
  structure(list(statistic = STATISTIC, df = DF, p.value = PVALUE, observed = obs, observed.grouped = nobs,
                 expected = esp, expected.grouped = nesp, residuals = sqrt(STATISTIC)))
}
