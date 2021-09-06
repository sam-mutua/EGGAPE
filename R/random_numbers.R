#' Title EGGAPE: Exponentiated Generalized Gull Alpha Power Exponential Distribution
#'
#' @param n Number of observations
#' @param a shape parameter of the distribution
#' @param b shape parameter of the distribution
#' @param alpha shape parameter of the distribution
#' @param lambda scale parameter of the distribution
#'
#' @return reggape generates random variates
#' @export
#'
#' @examples
#' reggape(100,1.5,1.8,0.4,0.5)
reggape <- function(n,a, b, alpha,lambda)
{
  if(n >= 1 & a > 0 & b > 0 & alpha > 0 & lambda > 0)
  {
    u <- stats::runif(n,0,1)
    fn <- ((1-u^(1/b))^(1/a))
    fg <- -1+fn
    jki <- (fg*log(alpha))/alpha
    gh <- gsl::lambert_Wm1(jki)
    oi <- log(alpha) + gh
    po <- log(alpha)/oi
    ft <- (1/lambda)*log(po)
    ft
  }
  else if (n <= 0 & a > 0 & b > 0 & alpha > 0 & lambda > 0){ stop("Error \"n\" must be greater than or equal to 1")
  }
  else stop("Ensure all the arguments \"n\",\"a\", \"b\", \"alpha\", \"lambda\" meet the required criterion ")
}

#' Title EGGAPE: Exponentiated Generalized Gull Alpha Power Exponential Distribution
#'
#' @param x vector of quantities
#' @param a shape parameter of the distribution
#' @param b shape parameter of the distribution
#' @param alpha shape parameter of the distribution
#' @param lambda scale parameter of the distribution
#'
#' @return deggape gives the density function
#' @export
#'
#' @examples
#' deggape(c(seq(0,1, by = 0.1)), 1.5,1.5,2.5,0.2)
deggape <- function(x,a, b, alpha, lambda)
{
  num <- ((a*b*alpha^(exp(-lambda*x)))*(lambda*exp(-lambda*x))) - (a*b*alpha^(exp(-lambda*x))) * log(alpha)*lambda*(1-exp(-lambda*x))*exp(-lambda*x)
  denom <- ((1-(1-(alpha*(1-exp(-lambda*x)))/(alpha^(1-exp(-lambda*x))))^a)^(b-1))*(1-(alpha*(1-exp(-lambda*x)))/(alpha^(1-exp(-lambda*x))))^(a-1)
  pdf <- num*denom
  pdf
}

#' Title EGGAPE: Exponentiated Generalized Gull Alpha Power Exponential Distribution
#'
#' @param p Vector of probabilities
#' @param a Shape parameter
#' @param b Shape parameter
#' @param alpha Shape parameter
#' @param lambda Scale parameter
#'
#' @return numeric
#' @export
#'
#' @examples
#' qeggape(0.1,0.3,0.4,0.4,0.5)
qeggape <- function(p,a, b, alpha,lambda)
{
  fn <- ((1-p^(1/b))^(1/a))
  fg <- -1+fn
  jki <- (fg*log(alpha))/alpha
  gh <- gsl::lambert_Wm1(jki)
  oi <- log(alpha) + gh
  po <- log(alpha)/oi
  ft <- (1/lambda)*log(po)
  ft
}
