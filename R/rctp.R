#' @rdname ctp
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @export
#'
#' @examples
#' rctp(4,1,1,3)

rctp<-function(n, a, b, g){
  if ( mode(c(n,a,b,g)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( g <= 2 * a )
    stop("gamma must be greater than 2a")

  icomplex <- sqrt(as.complex(-1))

  randoms <- runif(n, 0, 1)
  maxRandoms <- max(randoms)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(cgamma(g - a + b * icomplex, log = TRUE)) - lgamma(g) - lgamma(g - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating Density Distribution
  while( Fd[[i]] < maxRandoms ){
    pmfAux <- exp(log(pmfAux)+log(((a+i-1)^2+b^2))-log((g+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:length(randoms)){
    pMin=1
    pMax=length(Fd)
    while (Fd[pMin] < randoms[i]){
      mitad = floor ((pMin + pMax)/2)
      if (Fd[mitad] <= randoms[i] && Fd[mitad+1] >= randoms[i])
        pMin = mitad + 1
      else if (Fd[mitad] <= randoms[i])
          pMin = mitad
      else
          pMax = mitad
    }

    result[[i]]=pMin-1
  }
  return (result)
}


#' @rdname cbp
#' @importFrom fAsianOptions cgamma
#' @export
#'
#' @export
#'
#' @examples
#' rcbp(4,1,3)

rcbp<-function(n, b, g){
  if ( mode(c(n,b,g)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( g <= 0)
    stop("gamma must be greater than 0")

  icomplex <- sqrt(as.complex(-1))

  randoms <- runif(n, 0, 1)
  maxRandoms <- max(randoms)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(cgamma(g + b * icomplex, log = TRUE)) - lgamma(g) - lgamma(g)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating Density Distribution
  while( Fd[[i]] < maxRandoms ){
    pmfAux <- exp(log(pmfAux)+log(((i-1)^2+b^2))-log((g+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
  }
  #Searching values
  for (i in 1:length(randoms)){
    pMin=1
    pMax=length(Fd)

    while (Fd[pMin] < randoms[i]){
      mitad = floor ((pMin + pMax)/2)
      if (Fd[mitad] <= randoms[i] && Fd[mitad+1] >= randoms[i])
        pMin = mitad + 1
      else if (Fd[mitad] <= randoms[i])
        pMin = mitad
      else
        pMax = mitad
    }

    result[[i]]=pMin-1
  }
  return (result)
}




rctpPepe<-function(n, a, b, g,limit=10^(-10)){
  if ( mode(c(n,a,b,g)) != "numeric")
    stop( "non-numeric argument to mathematical function")

  if( g <= 2 * a )
    stop("gamma must be greater than 2a")

  icomplex <- sqrt(as.complex(-1))

  randoms <- runif(n, 0, 1)
  maxRandoms <- max(randoms)
  result<-vector(mode="numeric",length=n)

  lf0 <- 2 * Re(cgamma(g - a + b * icomplex, log = TRUE)) - lgamma(g) - lgamma(g - 2 * a)
  pmfAux<-exp(lf0)
  i=1
  Fd <-c(pmfAux)
  #Generating Density Distribution
  while( Fd[[i]] < maxRandoms ){
    pmfAux <- exp(log(pmfAux)+log(((a+i-1)^2+b^2))-log((g+i-1))-log(i))
    Fd <- c( Fd, Fd[[i]] + pmfAux )
    i <- i + 1
    if (Fd[[i]]>1-limit) Fd[[i]] <- 1
  }
  #Searching values
  for (i in 1:length(randoms)){
    pMin=1
    pMax=length(Fd)
    while (Fd[pMin] < randoms[i]){
      mitad = floor ((pMin + pMax)/2)
      if (Fd[mitad] <= randoms[i] && Fd[mitad+1] >= randoms[i])
        pMin = mitad + 1
      else if (Fd[mitad] <= randoms[i])
        pMin = mitad
      else
        pMax = mitad
    }

    result[[i]]=pMin-1
  }
  return (result)
}

