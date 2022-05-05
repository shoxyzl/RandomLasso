library(parallel)
library(glmnet)
library(tune)


################################################################################
############################# Utility Functions ################################
################################################################################


bootstrap <- function(X, y) {
  n <- dim(X)[1]
  bootstrap_ind <- sample(1:n, n, replace=TRUE)
  X.b <- X[bootstrap_ind,]
  y.b <- y[bootstrap_ind]
  return(list(X = X.b, y=y.b))
}


demean <- function(X, y, continuous=NULL) {
  # tell X and y that they are no-good wastes of space.
  n <- dim(X)[1]
  demeaned_X <- X
  feature_means <- apply(X, 2, mean)
  for (i in 1:n){
    row <- X[i,]
    demeaned_row <- X[i,] - feature_means
    demeaned_X[i,] <- demeaned_row
  }
  if (continuous & !missing(y)) {
    demeaned_y <- y - mean(y)
    return(list(X = demeaned_X, X_mean = feature_means,
                y = demeaned_y, y_mean = mean(y)))
  }
  else {
    return(list(X = demeaned_X, X_mean = feature_means))
  }

}


################################################################################
########################### Generic Regularized GLM ############################
################################################################################

regularized_glm <- function(X, y, alpha=0, continuous=T, exclude=NULL,intercept=T) {
  X <- as.matrix(X)

  if (continuous) {
    if(intercept){
      cv <- cv.glmnet(X, y, alpha=alpha, exclude=exclude)
      model <- glmnet(X, y, alpha=alpha, exclude=exclude, lambda=cv$lambda.min)
    }
    else{
      cv <- cv.glmnet(X, y, alpha=alpha, exclude=exclude)
      model <- glmnet(X, y, alpha=alpha, exclude=exclude,
                      lambda=cv$lambda.min, intercept=FALSE)
    }
  }
  else {
    if(intercept){
      y <- as.factor(y)
      cv <- cv.glmnet(X, y, alpha=alpha, exclude=exclude, family='binomial')
      model <- glmnet(X, y, alpha=alpha, exclude=exclude,
                      family='binomial', lambda=cv$lambda.min)
    }
    else{
      y <- as.factor(y)
      cv <- cv.glmnet(X, y, alpha=alpha, exclude=exclude, family='binomial')
      model <- glmnet(X, y, alpha=alpha, exclude=exclude,
                      family='binomial', lambda=cv$lambda.min,intercept=FALSE)
    }
  }
  return (model)
}

################################################################################
############################### Lasso & Ridge ##################################
################################################################################


#' lasoo
#'
#' Fit a generalized linear model by lasso regression.
#' Can deal with continuous and binary dependent variable Y.
#'
#' @param X independent variable
#' @param y dependent variable
#' @param continuous whether the dependent variable is continuous. Default is true.
#' @param exclude Indices of variables to be excluded from the model. Default is none.
#' @param intercept Should intercept(s) be fitted (default=TRUE) or set to zero (FALSE)
#'
#' @return
#' @export
#'
#' @examples lasso(X = matrix(rnorm(100 * 20), 100, 20), y = rnorm(100))
lasso <- function(X, y, continuous=TRUE, exclude=NULL,intercept=TRUE) {
  return(regularized_glm(X, y, 1, continuous, exclude,intercept))
}

ridge <- function(X, y, continuous=TRUE, exclude=NULL,intercept=TRUE) {
  return(regularized_glm(X, y, 0, continuous, exclude,intercept))
}

################################################################################
############################### Generate Step ##################################
################################################################################

random_lasso.generate_i <- function(X, y, q1, continuous) {
  bs <- bootstrap(X,y)

  q <- dim(X)[2]
  exclude <- sample(1:q, q-q1, replace=FALSE)

  beta_j_hat <- coef(lasso(bs$X, bs$y, continuous, exclude, intercept = F))
  return(as.vector(beta_j_hat)[-1])
}


random_lasso.generate <- function(X, y, B, q1, continuous) {
  q <- dim(X)[2]

  beta_hat <- matrix(NA, nrow=B, ncol=q)
  for (i in 1:B) {
    beta_hat[i,] <- t(as.matrix(random_lasso.generate_i(X, y, q1, continuous)))
  }

  importance <- apply(t(beta_hat), 1, function(x) abs(1/B*sum(x)))
  return(importance)
}

################################################################################
############################### Select Step ####################################
################################################################################


random_lasso.select_i <- function(X, y, B, q2, importance, continuous) {
  bs <- bootstrap(X,y)

  q <- dim(X)[2]
  include <- sample(1:q, q2, replace=FALSE, prob = importance)
  exclude <- (1:q)[-include]

  beta_j_hat <- coef(lasso(bs$X, bs$y, continuous, exclude,intercept = F))
  return(as.vector(beta_j_hat)[-1])
}

random_lasso.select <- function(X, y, B, q2, importance, continuous) {
  q <- dim(X)[2]

  beta_hat <- matrix(NA, nrow=B, ncol=q)
  for (i in 1:B) {
    beta_hat[i,] <- t(as.matrix(random_lasso.select_i(X, y, B, q2,
                                                      importance, continuous)))
  }

  beta_j_final <- apply(t(beta_hat), 1, function(x) 1/B * sum(x))
  return(beta_j_final)
}

################################################################################
############################### Random Lasso ###################################
################################################################################


random_lasso <- function(X, y, B, q1, q2, continuous) {
  # Step 1:
  #   - Create some number of bootstrap samples
  #   - Select a subset of covariates for each sample
  #     - randomly choosing some number of features "q1"
  #   - Apply linear regression with lasso regularization
  #     to each bootstrap sample with the selected q1 features.
  #   - Calculate an importance measure from the fitted coefficients.
  #
  # Step 2:
  #
  original_X <- X
  original_y <- y
  demeaned <- demean(X, y, continuous)
  X <- demeaned$X
  if (continuous) {
    y <- demeaned$y
  }
  importance <- random_lasso.generate(X, y, B, q1, continuous)
  beta_hat <- random_lasso.select(X, y, B, q2, importance, continuous)
  if (continuous) {
    return(list(X.original = original_X,
                X.demeaned = X,
                X.mean = demeaned$X_mean,
                y.original = original_y,
                y.demeaned = y,
                y.mean = demeaned$y_mean,
                beta_hat = beta_hat,
                model_params = list(continuous = continuous,
                                     B = B,
                                     q1 = q1,
                                     q2 = q2)
                ))
  }
  else {
    return(list(X.original = original_X,
                X.demeaned = X,
                X.mean = demeaned$X_mean,
                y.original = y,
                beta_hat = beta_hat,
                model_params = list(continuous = continuous,
                                    B = B,
                                    q1 = q1,
                                    q2 = q2)
    ))
  }
}


random_lasso.predict <- function(random_lasso_model, X) {
  # if (random_lasso_model$model_params$continuous) {}
  return (X %*% random_lasso_model$beta_hat)
}

inv_logit <- function(beta) {
  return( exp(sum(beta)) / (1 + exp(sum(beta))) > 0.5)
}




cv.random_lasso <- function(X, y, B, Q1, Q2, continuous) {
  for (q1 in Q1) {
    for (q2 in Q2) {

    }
  }
}

dat <- read.delim('http://www.ams.sunysb.edu/~pfkuan/Teaching/AMS597/Data/leukemiaDataSet.txt',
                  header=T,sep='\t')
X <- as.matrix(dat[,-1])
y<- dat[,1]

B = 100
q1 = 2
q2 = 3

beta_j_hat <- random_lasso(dat[,-1],dat[,1],B,q1,q2,F)

pred <- random_lasso.predict(beta_j_hat, X)
