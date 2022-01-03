library(mrfDepth) # para AO
library(e1071) # para svm

dirOut <- function(test, train, outly.method="SDO") {
  pTrain <- dim(train)[1]
  pTest <- dim(test)[1]
  n <- dim(train)[2]
  dirOut <- matrix(0, pTest, n)
  trainTestBind <- rbind(train,test)
  
  if (outly.method=="normal") {
    diroutly.method <- function(col) {(col[(pTrain+1):length(col)]-mean(col[1:pTrain]))/sd(col[1:pTrain])}
  }
  else if (outly.method=="SDO") {
    diroutly.method <- function(col) {(col[(pTrain+1):length(col)]-median(col[1:pTrain]))/mad(col[1:pTrain])}
  }
  else if (outly.method=="AO"){
    diroutly.method <- function(col) {adjOutl(col[1:pTrain],col[(pTrain+1):length(col)])$outlyingnessZ}
  }

  dirOut <- apply(trainTestBind, 2, diroutly.method)
  out_avr <- apply(dirOut, 1, mean)
  out_dev <- sqrt(apply(apply(dirOut, 2, (function(col) (col-out_avr)**2)),1,mean))
  M <- cbind(out_avr,out_dev)
  colnames(M) <- c("out_avr","out_dev")
  return(M)
}

clasificador <- function(test, train, y, outly.method="SDO", classifier="svm") {
  pTrain <- dim(train)[1]
  pTest <- dim(test)[1]
  levels <- levels(y)
  out_avrTrain <- matrix(0, pTrain, length(levels))
  out_avrTest <- matrix(0, pTest, length(levels))
  out_devTrain <- matrix(0, pTrain, length(levels))
  out_devTest <- matrix(0, pTest, length(levels))
  for (i in 1:length(levels)) {
    out_avrTrain[,i] <- dirOut(train,train[y==levels[i],], outly.method=outly.method)[,1]
    out_avrTest[,i] <- dirOut(test,train[y==levels[i],], outly.method=outly.method)[,1]
    out_devTrain[,i] <- dirOut(train,train[y==levels[i],], outly.method=outly.method)[,2]
    out_devTest[,i] <- dirOut(test,train[y==levels[i],], outly.method=outly.method)[,2]
  }
  if (classifier=="svm") {
    model <- svm(cbind(out_avrTrain,out_devTrain),y)
    predictions <- predict(model, cbind(out_avrTest,out_devTest))
  }
  else if (classifier=="qda") {
    model <- qda(cbind(out_avrTrain,out_devTrain),y)
    predictions <- predict(model, cbind(out_avrTest,out_devTest))$class
  }
  return(list(predictions = predictions,
              out_avrTrain = out_avrTrain, out_avrTest = out_avrTest,
              out_devTrain = out_devTrain, out_devTest = out_devTest))
}

# Ejemplo de uso
library(fda) # para growth
data<-growth
data <- as.data.frame(rbind(t(data$hgtm),t(data$hgtf)))
target <- as.factor(c(rep(0,39),rep(1,54)))
p <- dim(data)[1]
indicesTrain <- sample(1:p,floor(p*0.7))
train <- data[indicesTrain,]
test <- data[-indicesTrain,]
targetTrain <- target[indicesTrain]
targetTest <- target[-indicesTrain]
predichos<-clasificador(test,train,targetTrain,outly.method = "SDO")$predictions
mean(predichos==targetTest)


