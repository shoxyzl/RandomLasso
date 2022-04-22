randomlasso <- function(x,y,bootstrap){
  x <- as.matrix(x)
  y <- as.factor(y)
  n <- nrow(x)
  p <- ncol(x)
  if(missing(bootstrap)){
    B <- 500
  }
  else{
    B <- bootstrap
  }
  .part1 <- function(x,y){
    beta.hat <- replicate(p,0)
    id <- sample(1:n,n,replace = T)
    n.q1 <- runif(1,min = 1,max = p-1)
    q1 <- sample(p,n.q1,replace = F)
    x_hat <- as.matrix(x[id,q1])
    y_hat <- y[id]
    beta.hat[q1] <- lasso(x_hat,y_hat,continuous = F)
    return(list(beta.hat))
  }
  for(j in 1:B){
    beta1[j] <- .part1(x,y)
  }
  importance.measure <- abs(Reduce('+', beta1))
  return()
}
