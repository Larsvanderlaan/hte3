
library(tmvtnorm)



sim.CATEHighDim <- function(n, hard = TRUE, positivity = TRUE, ...) {
  d <- 20
  Sigma <- outer(1:d, 1:d, Vectorize(function(i,j) {
    if(i==j) {return(1)}
    else {
      return(0.4)
    }
  }))

  W <- rtmvnorm(n, mean = rep(0, d),
                sigma = Sigma,
                lower=rep(-2, d),
                upper=rep( 2, d), algorithm = "gibbs")/2


  colnames(W) <- paste0("W", 1:d)
  W1 <- W[,1]
  W5 <- W[,5]
  W9 <- W[,9]
  W11 <- W[,11]
  W19 <- W[,19]
  W10 <- W[,10]
  W15 <- W[,15]
  if(positivity){
    pA1W <- plogis((W1 + W5 + W9 + W11 + W19)/5  )
  } else {
    pA1W <- plogis((W1 + W5 + W9 + W11 + W19) /1.3  )

  }
  quantile(pA1W)
  A <- rbinom(n, 1 ,  pA1W)
  if(!hard) {
    CATE <- 1 + (W1 + W5 + W9 + W15 + W10)/5

  } else {
    CATE <- 1 + (sin(4*W1) + sin(4*W5) + cos(4*W9) + 1.5*(W15^2 - W10^2))/5

  }

  EY0W <- (cos(4*W1) + cos(4*W5) + sin(4*W9) + 1/(1.5+W15) + 1/(1.5+ W10))/5 +
    (sin(5*W1) + sin(5*W5) + sin(5*W9))/2

  EY1W <- EY0W + CATE
  Y <- rnorm(n, EY0W + A*CATE, 1)

  return(data.table(W, A, Y, EY0W, EY1W, pA1W))
}








