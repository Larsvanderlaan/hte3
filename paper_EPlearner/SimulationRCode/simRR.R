
# Sims CATE
## no positivity
## easy RR



sim.RR <- function(n, hard = TRUE, positivity = TRUE) {

  if(!positivity & !hard) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3)/3)
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    LRR <- 0.1*(-1 + W1 + W2 + W3)
    EY0W <- plogis(-1 + 0.3*(W1 + W2 + W3 + sin(4*W1) +  sin(4*W2) +  sin(4*W3)))
    #EY0W <- plogis(-0.8 + -1+W1 + sin(5*W2) + 1/(W3 + 2))
    EY1W <- pmin(EY0W * exp(LRR),1)
    Y <- rbinom(n, size = 1, EY0W * (1-A) + A* EY1W)
  }

  ## hardCATE
  if(!positivity & hard) {


    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3)/3)
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    LRR <- 0.1*(-1 + W1 + W2 + W3 + sin(4*W1) +  sin(4*W2) +  sin(4*W3))
    EY0W <- plogis(-1 + 0.3*(W1 + W2 + W3 + sin(4*W1) +  sin(4*W2) +  sin(4*W3)))
    EY1W <- pmin(EY0W * exp(LRR),1)
    Y <- rbinom(n, size = 1, EY0W * (1-A) + A* EY1W)


  }
  if(positivity & !hard) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3))
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    LRR <- 0.1*(-1 + W1 + W2 + W3)
    EY0W <- plogis(-1 + 0.3*(W1 + W2 + W3 + sin(4*W1) +  sin(4*W2) +  sin(4*W3)))
    EY1W <- pmin(EY0W * exp(LRR),1)
    Y <- rbinom(n, size = 1, EY0W * (1-A) + A* EY1W)
  }

  ## hardCATE
  if(positivity & hard) {


    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3))
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    LRR <- 0.1*(-1 + W1 + W2 + W3 + sin(4*W1) +  sin(4*W2) +  sin(4*W3))
    EY0W <- plogis(-1 + 0.3*(W1 + W2 + W3 + sin(4*W1) +  sin(4*W2) +  sin(4*W3)))
    EY1W <- pmin(EY0W * exp(LRR),1)
    Y <- rbinom(n, size = 1, EY0W * (1-A) + A* EY1W)


  }
  return(data.table(W, A, Y, pA1W, EY1W, EY0W))
}




