library(parallel)
library(glmnet)
# sample size n

X = matrix(rnorm(100 * 50), 100, 50)
y = rnorm(100)

# generic regularized regression function
regularized_glm <- function(X, y, alpha=0, continuous=T, exclude=NULL) {
  X <- as.matrix(X)
  if (continuous) {
    cv <- cv.glmnet(X, y, alpha=alpha, exclude=exclude)
    model <- glmnet(X, y, alpha=alpha, exclude=exclude, lambda=cv$lambda.min)
  } 
  else {
    y <- as.factor(y)
    cv <- cv.glmnet(X, y, alpha=alpha, exclude=exclude, family='binomial')
    model <- glmnet(X, y, alpha=alpha, exclude=exclude, family='binomial', lambda=cv$lambda.min)
  }
  return (model)
}


lasso <- function(X, y, continuous=T, exclude=NULL) {
  return(regularized_glm(X, y, 1, continuous, exclude))
}

ridge <- function(X, y, continuous=T, exclude=NULL) {
  return(regularized_glm(X, y, 0, continuous, exclude))
}


bootstrap <- function(X,y,n) {
  bootstrap_ind <- sample(1:n, n, replace=TRUE)
  return(list(X = X[bootstrap_ind,], y = y[bootstrap_ind]))
}


random_lasso <- function(X, y, B, q1, q2) {
  generate_step <- random_lasso.generate(X, y, B, q1)
  
}

random_lasso.generate <- function(X, y, B, q1) {
  n <- dim(X)[1]
  q <- dim(X)[2]
  beta_hat <- matrix(NA, nrow=B, ncol=q)
  B <- 100
  for (i in 1:B) {
    beta_hat[i,] <- t(as.matrix(random_lasso.generate_i(X, y, q1, T)))
  }
  imp_j <- function (beta_j_hat) abs(1/B*sum(beta_j_hat))
  importance <- apply(t(beta_hat), 1, imp_j)
  return(list(beta_hat = beta_hat, importance = importance))
}

random_lasso.generate_i <- function(X, y, q1, continuous) {
  # 1. Draw Bootstrap Sample
  # 2. select q1 candidate variables
  n <- dim(X)[1]
  q <- dim(X)[2]
  
  # Bootstrap
  bootstrap_ind <- sample(1:n, n, replace=TRUE)
  X.b <- X[bootstrap_ind,]
  y.b <- y[bootstrap_ind]
  exclude <- sample(1:q, q-q1, replace=FALSE)
  beta_j_hat <- coef(lasso(X.b, y.b, T, exclude))
  return(as.vector(beta_j_hat[2:(q+1),]))
}

random_lasso.select <- function(X, y, B, q2, importance, continuous) {
  
}

random_lasso.select_i <- function(X, y, B, q2, importance, continuous) {
  n <- dim(X)[1]
  q <- dim(X)[2]
  X.b <- X[bootstrap_ind,]
  y.b <- y[bootstrap_ind]
  include <- sample(1:q, q, replace=FALSE, prob = importance)
  beta_j_hat <- coef(lasso(X.b, y.b, T, -include))
  return(as.vector(beta_j_hat[2:(q+1),]))
}



q1 = 20
imp <- random_lasso.generate_i(X,y,10,q1,T)
imp
imp