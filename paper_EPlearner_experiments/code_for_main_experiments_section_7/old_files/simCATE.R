

sim.CATE <- function(n, hard = TRUE, positivity = TRUE, randomized = F,  ...) {
  if(randomized) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- 0.5
    A <- rbinom(n, 1 ,  pA1W)
    CATE <- 1 + W1
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 1)
    return(data.table(W, A, Y, EY0W, EY1W, pA1W))

  }
  if(!positivity & !hard) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3)/3)
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    CATE <- 1 + W1
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 1)
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
    CATE <- 1 + rowSums(as.matrix(W)) + rowSums(sin(5*as.matrix(W)))
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 1)


  }
  ##  positivity
  ## easy CATE
  if(positivity & !hard) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3))
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    CATE <- 1 + rowSums(as.matrix(W))
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 1)
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
    CATE <- 1 + rowSums(as.matrix(W)) + rowSums(sin(5*as.matrix(W)))
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 1)

  }
  return(data.table(W, A, Y, EY0W, EY1W, pA1W))
}
