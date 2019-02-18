#' The Complex Triparametric Pearson (CTP) Distribution
#'
#' @description
#' Probability mass function for the Complex Triparametric Pearson (CTP) distribution with parameters a, b and gamma.
#'
#' @details bla blala lakj akj j fkjlh jakh bjksh jkh ajbh akjh kajsh jb ajbfd
#'   alskjfdaj sdfhlaskjhf lasdkjhf kjasdhfkjadshfkjahds fkjh asdkjfh asdjkhf
#'   kajf sdhlkasjdhf lkasdjhf lkasdjhf kasdjh flkasjh flksajdh flkasjhd jhasdf l
#'
#' @param x vector of (non-negative integer) quantiles
#' @param p vector of (non-negative integer) probabilities
#' @param n number of observations. If length(n) > 1, the length is
#'  taken to be the number required.
#' @param a parameter a (real)
#' @param b parameter b (real)
#' @param g parameter gamma (positive)
#' @param lower.tail A param g
#' @name ctp
#'
#' @return The CTP distribution with parameters a, b and gamma has pmf
#' \eqn{f(x|a,b,\gamma)=f_0\frac{(a+ib)_x(a-ib)_x}{(\gamma)_x}\frac{1}{x!},\quad x=0,1,...}
#' with \eqn{f_0=\frac{\Gamma(\gamma-a-ib)\Gamma(\gamma-a+ib)}{\Gamma(\gamma)\Gamma(\gamma-2a)}} the normalizing constant.
#' If a=0 the CTP is a Complex Biparametric Pearson (CBP) distribution, so the pmf of the CBP distribution is obtained.

#'
#' @source
#' Rodr\'iguez-Avi, J., S\'aez-Castillo, A. J., Conde-Sánchez, A. and Olmo-Jim\'enez, M. J. (2004). A triparametric discrete distribution with complex parameters. Statistical Papers, 45 (1), 81-95.
#' Olmo-Jim\'enez, M. J., Rodr\'iguez-Avi, J. and Cueva-López, V. (2018). A review of the CTP distribution: a comparison with other over- and underdispersed count data models. Journal of Statistical Computation and Simulation, 88 (14), 2684-2706.
#'
NULL
#> NULL



#' @rdname ctp
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @examples
#' #Ejemplo uso dctp
#' dctp(3,1,2,5)
#' dctp(c(3,4),1,2,5)
#'
#'
dctp <- function(x, a, b, g) {

  if ( mode(c(x,a,b,g)) != "numeric")
    stop( "non-numeric argument to mathematical function" )

  if( g <= 2 * a )
    stop("gamma must be greater than 2a")

  x <- as.vector(x)
  errors<-c()


  for ( i in 1:length(x) ){
      if ( round(x[[i]]) != x[[i]] || x[[i]]<0){
        warning( paste ("non-integer or negative x[[", i, "]]", sep = ""))
        errors <- c(errors,i)
      }

  }

  if (b==0){
     lpmf<-2*(lgamma(g-a)+lgamma(x+a)-lgamma(a))-lgamma(g-2*a)-lgamma(g+x)-lgamma(x+1)
  } else {
     i <- sqrt(as.complex(-1))
     lpmf<-Re(cgamma(g-a-b*i,log=TRUE)+cgamma(g-a+b*i,log=TRUE)+cgamma(x+a+b*i,log=TRUE)+cgamma(x+a-b*i,log=TRUE)) -lgamma(g-2*a)-Re(cgamma(a+b*i,log=TRUE)+cgamma(a-b*i,log=TRUE))-lgamma(g+x)-lgamma(x+1)
  }
  result <- exp(lpmf)
  for (i in 1:length(errors)){
    result[ errors[[i]] ] = 0
  }
  return(result)
}


#' The Complex Triparametric Pearson (CTP) Distribution
#'
#' @description
#' Probability mass function for the Complex Triparametric Pearson (CTP) distribution with parameters a, b and gamma.
#'
#' @details bla blala lakj akj j fkjlh jakh bjksh jkh ajbh akjh kajsh jb ajbfd
#'   alskjfdaj sdfhlaskjhf lasdkjhf kjasdhfkjadshfkjahds fkjh asdkjfh asdjkhf
#'   kajf sdhlkasjdhf lkasdjhf lkasdjhf kasdjh flkasjh flksajdh flkasjhd jhasdf l
#'
#' @param x vector of (non-negative integer) quantiles
#' @param p vector of (non-negative integer) probabilities
#' @param n number of observations. If length(n) > 1, the length is
#'  taken to be the number required.
#' @param b parameter b (real)
#' @param g parameter gamma (positive)
#' @param lower.tail A param g
#' @name cbp
#'
#' @return The CTP distribution with parameters a, b and gamma has pmf
#' \eqn{f(x|a,b,\gamma)=f_0\frac{(a+ib)_x(a-ib)_x}{(\gamma)_x}\frac{1}{x!},\quad x=0,1,...}
#' with \eqn{f_0=\frac{\Gamma(\gamma-a-ib)\Gamma(\gamma-a+ib)}{\Gamma(\gamma)\Gamma(\gamma-2a)}} the normalizing constant.
#' If a=0 the CTP is a Complex Biparametric Pearson (CBP) distribution, so the pmf of the CBP distribution is obtained.

#'
#' @source
#' Rodr\'iguez-Avi, J., S\'aez-Castillo, A. J., Conde-Sánchez, A. and Olmo-Jim\'enez, M. J. (2004). A triparametric discrete distribution with complex parameters. Statistical Papers, 45 (1), 81-95.
#' Olmo-Jim\'enez, M. J., Rodr\'iguez-Avi, J. and Cueva-López, V. (2018). A review of the CTP distribution: a comparison with other over- and underdispersed count data models. Journal of Statistical Computation and Simulation, 88 (14), 2684-2706.
#'
NULL
#> NULL

#' @rdname cbp
#' @importFrom fAsianOptions cgamma
#' @export
#' @examples
#' #Ejemplo uso dctp
#' dcbp(3,2,3)
#' dcbp(c(3,4),2,3)
#'
#'
dcbp <- function(x, b, g) {

  if ( mode(c(x,b,g)) != "numeric")
    stop( "non-numeric argument to mathematical function" )

  if( g <= 0)
    stop("gamma must be greater than 0")

  x <- as.vector(x)
  errors<-c()


  for ( i in 1:length(x) ){
    if ( round(x[[i]]) != x[[i]] || x[[i]]<0){
      warning( paste ("non-integer or negative x[[", i, "]]", sep = ""))
      errors <- c(errors,i)
    }

  }

  if (b==0){
    lpmf<-2*(lgamma(g)+lgamma(x))-lgamma(g)-lgamma(g+x)-lgamma(x+1)
  } else {
    i <- sqrt(as.complex(-1))
    lpmf<-Re(cgamma(g-b*i,log=TRUE)+cgamma(g+b*i,log=TRUE)+cgamma(x+b*i,log=TRUE)+cgamma(x-b*i,log=TRUE)) -lgamma(g)-Re(cgamma(b*i,log=TRUE)+cgamma(b*i,log=TRUE))-lgamma(g+x)-lgamma(x+1)
  }
  result <- exp(lpmf)
  for (i in 1:length(errors)){
    result[ errors[[i]] ] = 0
  }
  return(result)
}
