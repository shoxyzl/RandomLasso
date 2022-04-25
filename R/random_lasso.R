library(parallel)
library(glmnet)
# sample size n



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

random_lasso.generate_i <- function(X, y, q1, continuous) {
  
  n <- dim(X)[1]
  q <- dim(X)[2]
  
  bootstrap_ind <- sample(1:n, n, replace=TRUE)
  X.b <- X[bootstrap_ind,]
  y.b <- y[bootstrap_ind]
  
  exclude <- sample(1:q, q-q1, replace=FALSE)
  
  beta_j_hat <- coef(lasso(X.b, y.b, continuous, exclude))
  return(as.vector(beta_j_hat[2:(q+1),]))
}


random_lasso.generate <- function(X, y, B, q1, continuous) {
  n <- dim(X)[1]
  q <- dim(X)[2]
  beta_hat <- matrix(NA, nrow=B, ncol=q)
  for (i in 1:B) {
    beta_hat[i,] <- t(as.matrix(random_lasso.generate_i(X, y, q1, continuous)))
  }
  importance <- apply(t(beta_hat), 1, function(x) abs(1/B*sum(x)))
  return(importance)
}

random_lasso.select_i <- function(X, y, B, q2, importance, continuous) {
  n <- dim(X)[1]
  q <- dim(X)[2]
  
  bootstrap_ind <- sample(1:n, n, replace=TRUE)
  X.b <- X[bootstrap_ind,]
  y.b <- y[bootstrap_ind]
  
  include <- sample(1:q, q2, replace=FALSE, prob = importance)
  exclude <- (1:q)[-include]
  
  beta_j_hat <- coef(lasso(X.b, y.b, continuous, exclude))
  return(as.vector(beta_j_hat[2:(q+1),]))
}

random_lasso.select <- function(X, y, B, q2, importance, continuous) {
  
  n <- dim(X)[1]
  q <- dim(X)[2]
  
  beta_hat <- matrix(NA, nrow=B, ncol=q)
  for (i in 1:B) {
    beta_hat[i,] <- t(as.matrix(random_lasso.select_i(X, y, B, q2, importance, continuous)))
  }
  
  beta_j_final <- apply(t(beta_hat), 1, function(x) 1/B * sum(x))
  return(beta_j_final)
}

random_lasso <- function(X, y, B, q1, q2, continuous) {
  importance <- random_lasso.generate(X, y, B, q1, continuous)
  beta_j_hat <- random_lasso.select(X, y, B, q2, importance, continuous)
  return(beta_j_hat)
}



X = matrix(rnorm(100 * 50), 100, 50)
y = rnorm(100)
B = 100
q1 = 20
q2 = 20

beta_j_hat <- random_lasso(X,y,B,q1,q2,T)
