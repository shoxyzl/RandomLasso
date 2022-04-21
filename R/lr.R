## lasso regression
lr <- function(x,y,lambda=c("min","1st"),continuous=T){
  X <- data.matrix(x)
  if(continuous){
    cv.lambda <- glmnet::cv.glmnet(x,y,alpha=1)
    if(lambda=='min'){
      fit <- glmnet::glmnet(X,y,alpha = 1, lambda = cv.lambda$lambda.min)
    }
    else{
      fit <- glmnet::glmnet(X,y,alpha = 1, lambda = cv.lambda$lambda.1st)
    }
  }
  else{
    y.factor <- as.factor(y)
    cv.lambda <- glmnet::cv.glmnet(data.matrix(x),y.factor,alpha=1,family='binomial')
    if(lambda=='min'){
      fit <- glmnet::glmnet(X,y.factor, alpha = 1,family="binomial",lambda = cv.lambda$lambda.min)
    }
    else{
      fit <- glmnet::glmnet(X,y.factor, alpha = 1,family="binomial",lambda = cv.lambda$lambda.1st)
    }
  }
  c <- coef(fit)
  variables <- row.names(c)[c[,1]!=0]
  variables <- variables[!(variables %in% '(Intercept)')]
  return(variables)
}

