#' @rdname ctp
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @examples
#' pctp(3,1,2,3)
#' pctp(c(3,4),1,2,3)

pctp <- function(x, a, b, g, lower.tail = TRUE) {
  if ( mode(c(x,a,b,g)) != "numeric")
    stop( "non-numeric argument to mathematical function" )

  if( g <= 2 * a )
    stop( "gamma must be greater than 2a" )

  x<-as.vector(x)
  maxX=x[1]
  n<-length(x)
  for ( i in 1:n ){
    x[i] = floor( x[i] )
    if (x[i] > maxX)
      maxX=x[i]
  }
  icomplex <- sqrt(as.complex(-1))
  lf0 <- 2 * Re(cgamma(g - a + b * icomplex, log = TRUE)) - lgamma(g) - lgamma(g - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating distribution function
  digits=options()$digit
  while( i <= maxX+1 && Fd[i]<(1-10^-digits)){
    pmfAux <- exp(log(pmfAux)+log(((a+i-1)^2+b^2))-log((g+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    print(paste(i,Fd[i],pmfAux,sep=":"))
    i <- i + 1
  }

  result<-vector(mode="numeric",length=n)
  for ( i in 1:n ){
    if ( x[i] < 0 )
      result[i]=0
    else
      result[i]=Fd[x[i]+1]

    if (! lower.tail){
      result[i]<-1-result[i]
    }
  }

  return (result)
}

#' @rdname cbp
#' @importFrom fAsianOptions cgamma
#' @export
#'
#'
#' @examples
#' pcbp(3,2,3)
#' pcbp(c(3,4),2,3)

pcbp <- function(x, b, g, lower.tail = TRUE ) {
  if ( mode(c(x,b,g)) != "numeric")
    stop( "non-numeric argument to mathematical function" )

  if( g <= 0 )
    stop( "gamma must be greater than 0" )

  x<-as.vector(x)
  maxX=x[1]
  n<-length(x)
  for ( i in 1:n ){
    x[i] = floor( x[i] )
    if (x[i] > maxX)
      maxX=x[i]
  }
  icomplex <- sqrt(as.complex(-1))
  lf0 <- 2 * Re(cgamma(g  + b * icomplex, log = TRUE)) - lgamma(g) - lgamma(g )
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating distribution function
  while( i <= maxX+1 ){
    pmfAux <- exp(log(pmfAux)+log(((i-1)^2+b^2))-log((g+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }

  result<-vector(mode="numeric",length=n)
  for ( i in 1:n ){
    if ( x[i] < 0 )
      result[i]=0
    else
      result[i]=Fd[x[i]+1]

    if (! lower.tail){
      result[i]<-1-result[i]
    }
  }

  return (result)
}
