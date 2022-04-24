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


bootstrap <- function(X,n) {
  X.bootstrap_ind <- sample(1:n, n, replace=TRUE)
  return(X[X.bootstrap_ind,])
}


random_lasso <- function(X, B, q1, q2) {
  
}

random_lasso.generate <- function(X, B, q1) {
  
}

random_lasso.generate_i <- function(X, y, B, q1, continuous) {
  # 1. Draw Bootstrap Sample
  # 2. select q1 candidate variables
  n <- dim(X)[1]
  q <- dim(X)[2]
  X.bootstrap <- bootstrap(X, n)
  exclude <- sample(1:q, q-q1, replace=FALSE)
  beta_j_hat <- coef(lasso(X.bootstrap,y,T,exclude))
  importance <- abs(1/B*sum(beta_j_hat))
  return(list(beta_j_hat = beta_j_hat, importance = importance))
}



q1 = 20
random_lasso.generate_i(X,y,10,q1,T)
