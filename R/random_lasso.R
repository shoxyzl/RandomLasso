# Authors:
# Zhilin Yang & Branden Ciranni

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

inv_logit <- function(y_hat) {
  return( as.numeric(exp(y_hat) / (1 + exp(y_hat)) > 0.5))
}

classification_accuracy <- function(y, y_hat) {
  return(sum(y==y_hat)/length(y))
}

regression_mse <- function(y, y_hat) {
  return(sum((y-y_hat)^2)/length(y))
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


#' lasso
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


#' Random Lasso
#'
#' @param X Matrix of covariates
#' @param y array-like of Variable (Continuous or Binary)
#' @param B Number of Bootstrap Samples
#' @param q1 Number of random features to choose in generate step.
#' @param q2 Number of features to choose based on probabilities in select step.
#' @param continuous True if y is continuous, else False and Binary
#'
#' @return Beta Coefficient Values for fitted model
#' @export
#'
#' @examples 
#' 
#' data('QuickStartExample')
#'
#' B = 100
#' q1 = 10
#' q2 = 5
#' X <- QuickStartExample$x
#' y <- QuickStartExample$y
#' continuous = T
#' 
#' rlm <- random_lasso(X, y, B, q1, q2, continuous)
random_lasso <- function(X, y, B, q1, q2, continuous) {
  # original values before mean correction
  original_X <- X
  original_y <- y
  
  # peform mean correction as suggested in the paper
  demeaned <- demean(X, y, continuous)
  X <- demeaned$X
  if (continuous) {
    y <- demeaned$y
  }
  
  # calculate importance values from step 1: generate method
  importance <- random_lasso.generate(X, y, B, q1, continuous)
  
  # use importance values to select q2 features and fit model with these.
  beta_hat <- random_lasso.select(X, y, B, q2, importance, continuous)
  
  # return list with information relevant to the model for prediction and
  # diagnostics
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


################################################################################
################################# Prediction ###################################
################################################################################


#' Random Lasso Prediction
#'
#' @param random_lasso_model Fitted Random Lasso Model from `random_lasso` above
#' @param X Matrix of Covariates
#'
#' @return Either the predicted y_hat value for continuous targets, or a binary
#' logit output for binary targets.
#' @export
#'
#' @examples
#' data('BinomialExample')
#'
#' X <- BinomialExample$x
#' y <- BinomialExample$y
#' 
#' B = 100
#' q1 = 15
#' q2 = 10
#' continuous = F
#' 
#' 
#' rlm <- random_lasso(X, y, B, q1, q2, continuous)
#' random_lasso.predict(rlm, X)
random_lasso.predict <- function(random_lasso_model, X) {
  y_hat <- X %*% random_lasso_model$beta_hat
  if (random_lasso_model$model_params$continuous) {
    return (y_hat)
  }
  else {
    return (inv_logit(y_hat))
  }
}


################################################################################
########################### q1 & q2 Cross Validation ###########################
################################################################################


#' Cross Validation to find optimal q1 & q2 values.
#'
#' @param X Matrix of Covariates
#' @param y array-like of target variable
#' @param B Number of Bootstrap Samples
#' @param Q1 Vector of q1 values for grid search
#' @param Q2 Vector of q2 values for grid search
#' @param continuous True if y is continuous, else False and Binary.
#' @param training_size Ratio of Train/Test Split.
#'
#' @return Performance Matrix containing either MSE for Continuous Y, or
#' Prediction Accuracy for Binary Y.
#' @export
#'
#' @examples
#' 
#' data('QuickStartExample')
#'
#' B = 100
#' X <- QuickStartExample$x
#' y <- QuickStartExample$y
#' 
#' 
#' perf_mat <- cv.random_lasso(X, y, B, 
#'                              Q1=c(10, 20), 
#'                              Q2=c(5, 10, 15), 
#'                              continuous=T, 
#'                              training_size=0.7)
#'                              
#' perf_mat
#'          [,1]     [,2]     [,3]
#' [1,] 3.189807 1.229554 1.382035
#' [2,] 2.252646 1.223567 1.331062
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
      
        # Split data into train/test sets
        trainID <- sample(1:n,round(training_size*n))
        X_train <- X[trainID,]
        y_train <- y[trainID]
        X_test <- X[-trainID,]
        y_test <- y[-trainID]
        
        # fit rlm model and predict
        rlm <- random_lasso(X_train, y_train, B, q1, q2, continuous)
        y_hat <- random_lasso.predict(rlm, X_test)
        
        
        # calculate performance metric
        if (continuous==F) {
          performance <- classification_accuracy(y_test, y_hat)
        }
        else {
          performance <- regression_mse(y_test, y_hat)
        }
        
        # populate performance matrix
        perf[i,j] <- performance
    }
  }
  return (perf)
}


