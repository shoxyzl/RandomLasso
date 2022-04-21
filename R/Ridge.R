## ridge regression
ridge <- function(x,y,lambda,continuous=T){
  x <- as.matrix(x)
  if(continuous){
    if(missing(lambda)){
      cv.lambda <- glmnet::cv.glmnet(x,y,alpha=0)
      fit <- glmnet::glmnet(x,y,alpha = 0, lambda = cv.lambda$lambda.min)
    }
    else{
      fit <- glmnet::glmnet(x,y,alpha = 0, lambda = lambda)
    }
  }
  else{
    y<- as.factor(y)
    if(missing(lambda)){
      cv.lambda <- glmnet::cv.glmnet(x,y,alpha=0,family='binomial')
      fit <- glmnet::glmnet(x,y, alpha = 0,family="binomial",lambda = cv.lambda$lambda.min)
    }
    else{
      fit <- glmnet::glmnet(x,y, alpha = 0,family="binomial",lambda = lambda)
    }
  }
  beta.hat <- coef(fit)
  return(beta.hat)
}
