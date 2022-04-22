## lasso regression
lasso <- function(x,y,continuous=T){
  x <- as.matrix(x)
  if(continuous){
    cv.lambda <- glmnet::cv.glmnet(x,y,alpha=1,nfolds = 5)
    fit <- glmnet::glmnet(X,y,alpha = 1, lambda = cv.lambda$lambda.min)
  }
  else{
    y.factor <- as.factor(y)
    cv.lambda <- glmnet::cv.glmnet(x,y.factor,alpha=1,family='binomial',nfolds = 5)
    fit <- glmnet::glmnet(x,y.factor, alpha = 1,family="binomial",lambda = cv.lambda$lambda.min,intercept = F,standardize = FALSE)
  }
  c <- coef(fit)[-1]
  return(c)
}
