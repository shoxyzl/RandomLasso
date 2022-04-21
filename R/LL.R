## linear and logistic
ll <- function(y,x,logistic = F){
  if(logistic){
    y.factor <- as.factor(y)
    data <- data.frame(y.factor,x)
    fit <- lm(y.factor~.,data = data)
  }
  else{
    data <- data.frame(y,x)
    fit <- lm(y~., data=data)
  }
  return(summary(fit))
}
