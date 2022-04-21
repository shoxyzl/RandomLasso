## ridge regression
rr <- function(x,y,lambda=c("min","1st"),continuous=T){
  if(continuous){
    cv.lambda <- glmnet::cv.glmnet(x,y,alpha=0)
    if(lambda=='min'){
      fit <- glmnet::glmnet(x,y,alpha = 0, lambda = cv.lambda$lambda.min)
    }
    else{
      fit <- glmnet::glmnet(x,y,alpha = 0, lambda = cv.lambda$lambda.1st)
    }
  }
  else{
    y.factor <- as.factor(y)
    cv.lambda <- glmnet::cv.glmnet(data.matrix(x),y.factor,alpha=0,family='binomial')
    if(lambda=='min'){
      fit <- glmnet::glmnet(data.matrix(x),y.factor, alpha = 0,family="binomial",lambda = cv.lambda$lambda.min)
    }
    else{
      fit <- glmnet::glmnet(data.matrix(x),y.factor, alpha = 0,family="binomial",lambda = cv.lambda$lambda.1st)
    }
  }
  return(coef(fit))
}
