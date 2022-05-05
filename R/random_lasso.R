library(parallel)
library(glmnet)


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


inv_logit <- function(y_hat) {
  return( as.numeric(exp(y_hat) / (1 + exp(y_hat)) > 0.5))
}

random_lasso.predict <- function(random_lasso_model, X) {
  y_hat <- X %*% random_lasso_model$beta_hat
  if (random_lasso_model$model_params$continuous) {
    return (y_hat)
  }
  else {
    return (inv_logit(y_hat))
  }
}

classification_accuracy <- function(y, y_hat) {
  return(sum(y==y_hat)/length(y))
}

regression_mse <- function(y, y_hat) {
  return(sum((y-y_hat)^2)/length(y))
}


cv.random_lasso <- function(X, y, B, Q1, Q2, continuous, training_size=0.7) {
  I = length(Q1)
  J = length(Q2)
  n = dim(X)[1]
  
  if (continuous) {
    perf <- matrix(Inf, nrow=I, ncol=J)
  }
  else {
    perf <- matrix(0, nrow=I, ncol=J)
  }
  for (i in 1:I) {
    q1 = Q1[i]
    for (j in 1:J) {
      q2 = Q2[j]
        trainID <- sample(1:n,round(training_size*n))
        X_train <- X[trainID,]
        y_train <- y[trainID]
        X_test <- X[-trainID,]
        y_test <- y[-trainID]
        
        rlm <- random_lasso(X_train, y_train, B, q1, q2, continuous)
        y_hat <- random_lasso.predict(rlm, X_test)
        
        
        if (continuous==F) {
          performance <- classification_accuracy(y_test, y_hat)
        }
        else {
          performance <- regression_mse(y_test, y_hat)
        }
        
        perf[i,j] <- performance
    }
  }
  return (perf)
}

data('BinomialExample')
data('QuickStartExample')

X <- BinomialExample$x
y <- BinomialExample$y

X <- QuickStartExample$x
y <- QuickStartExample$y



# dat <- read.delim('http://www.ams.sunysb.edu/~pfkuan/Teaching/AMS597/Data/leukemiaDataSet.txt',
#                   header=T,sep='\t')
#  <- as.matrix(dat[,-1])
# y <- dat[,1]
# as.numeric(y)
# y[which(y==y[1])] = 1
# y[which(y!=y[1])] = 0
# y <- as.numeric(y)

# y

B = 100

perf_mat <- cv.random_lasso(X, y, B, Q1=c(10, 20), Q2=c(5, 10, 15), continuous=T, training_size=0.7)
