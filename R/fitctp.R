#Ajuste por maxima-verosimilitud de la distribucion CTP(a,b,gamma)
#a partir de un conjunto de datos SIN AGRUPAR
#a y b son reales, gamma > 0
#' Maximum likelihood estimation fit for CTP(a,b,gamma)
#'
#' @param x A param x
#' @param astart A param astart
#' @param bstart A param bstart
#' @param gammastart A param gammastart
#' @param method A param method ......
#' @param moments A param TRUE FALSE that say if we want to use initial values generating for moments method if posible.
#' @param hessian A param hessian TRUE FALSE
#' @param control A param control
#' @param ...  other params
#'
#' @return la salida del ...
#' @export
#'
#' @examples
#'x<-c(3,5,7,2,6,9,7,5)
#'fitctp(x,1,1.1,3)
#'
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
  #Unidad imaginaria
  i<-sqrt(as.complex(-1))

  #Definimos la log-verosimilitud
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


  #Optimizamos la log-verosimilitud

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


  #Result
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


#' Maximum likelihood estimation fit for CBP(b,gamma)
#'
#' @param x A param x
#' @param bstart A param bstart
#' @param gammastart A param gammastart
#' @param method A param method ......
#' @param moments A param TRUE FALSE that say if we want to use initial values generating for moments method if posible.
#' @param hessian A param hessian TRUE FALSE
#' @param control A param control
#' @param ...  other params
#'
#' @return la salida del ...
#' @export
#'
#' @examples
#'x<-c(3,5,7,2,6,9,7,5)
#'fitcbp(x,1.1,3)
#'
fitcbp <- function(x, bstart = 1, gammastart = 1.1, method = "L-BFGS-B", moments=FALSE, hessian = TRUE, control = list(), ...)
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
  #Unidad imaginaria
  i<-sqrt(as.complex(-1))

  #Definimos la log-verosimilitud
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


  #Optimizamos la log-verosimilitud

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
    dimnames(coef.table)<-list("",c("b","log(gamma)"))
  }
  else {
    coef.table<-rbind(fit$par,deparse.level=0)
    dimnames(coef.table)<-list("",c("b","gamma"))
  }


  #Result
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
  class(results) <- "fitCTP"
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
    cat("Coefficients")
    cat(":\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")

  if(!x$converged){
    cat("Error of convergence")
  }
  invisible(x)
}
