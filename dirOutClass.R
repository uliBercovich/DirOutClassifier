library(mrfDepth) # AO
library(e1071) # svm
library(MASS) # qda
library(class) # knn

dirOut <- function(test, train, outly.method, L.norm, depth) {

  pTrain <- dim(train)[1]
  pTest <- dim(test)[1]
  n <- dim(train)[2]
  
  ### Univariate method for outlyingness calculation
  if (outly.method=="Mahalanobis") {
    Zcenter <- apply(train,2,mean)
    Zsd <- apply(train,2,sd)
    dirOutMatrix <- (test-rep(Zcenter,each=pTest))/rep(Zsd,each=pTest)
  }
  else if (outly.method=="SDO") {
    Zcenter <- apply(train,2,median)
    Zmad <- apply(train,2,mad)
    dirOutMatrix <- (test-rep(Zcenter,each=pTest))/rep(Zmad,each=pTest) 
    }
  else if (outly.method=="AO"){
    Zcenter <- apply(train,2,median)
    adjOutlCol <- apply(test,2,function(columna) {adjOutl(columna[1:pTrain],columna[(pTrain+1):pTest])})
    dirOutMatrix <- sapply(1:n,function(ncol) {c(adjOutlCol[[ncol]]$outlyingnessX,adjOutlCol[[ncol]]$outlyingnessZ)})
    dirOutMatrix <- dirOutMatrix*sign(test-rep(Zcenter,each=pTest))
    }
  
  
  ### False for outlyingness, True for depth
  if (depth==TRUE) {
    dirOutMatrix <- 1/(1+abs(dirOutMatrix))
  }
  out_avr <- apply(dirOutMatrix, 1, mean)
  
  ### L norm to calculate variation
  if (L.norm=="Inf") {
    out_var <- apply(matrix(apply(dirOutMatrix, 2, (function(columna) abs(columna-out_avr))),pTest),1,max)
  }
  else {
    out_var <- (apply(matrix(apply(dirOutMatrix, 2, (function(columna) (abs(columna-out_avr))**L.norm)),pTest),1,mean))**(1/L.norm)
  }
  M <- cbind(out_avr,out_var)
  colnames(M) <- c("out_avr","out_var")
  return(M)
}


dirOutClass <- function(test, train, y, outly.method="SDO", classifier="knn", L.norm=2, depth=FALSE, out.dev=TRUE) {
 
  pTrain <- dim(train)[1]
  pTest <- dim(test)[1]
  levels <- levels(y)
  trainTestBind <- rbind(train,test)
  
  ### Dimensionality reduction
  out_avrTrain <- matrix(0, pTrain, length(levels))
  out_avrTest <- matrix(0, pTest, length(levels))
  out_varTrain <- matrix(0, pTrain, length(levels))
  out_varTest <- matrix(0, pTest, length(levels))
  for (i in 1:length(levels)) {
    M <- dirOut(trainTestBind,train[y==levels[i],], outly.method=outly.method, L.norm=L.norm, depth=depth)
    out_avrTrain[1:pTrain,i] <- M[1:pTrain,1]
    out_avrTest[1:pTest,i] <- M[(pTrain+1):(pTrain+pTest),1]
    out_varTrain[1:pTrain,i] <- M[1:pTrain,2]
    out_varTest[1:pTest,i] <- M[(pTrain+1):(pTrain+pTest),2]
  }
  
  ### Classify 
  if (classifier=="svm") {
    if (out.dev==FALSE) {
      model <- svm(out_avrTrain,y)
      predictions <- predict(model, out_avrTest)
    }
    else {
      model <- svm(cbind(out_avrTrain,out_varTrain),y)
      predictions <- predict(model, cbind(out_avrTest,out_varTest))
    }
  }
  else if (classifier=="qda") {
    if (out.dev==FALSE) {
      model <- qda(out_avrTrain,y)
      predictions <- predict(model, out_avrTest)$class
    }
    else {
      model <- qda(cbind(out_avrTrain,out_varTrain),y)
      predictions <- predict(model, cbind(out_avrTest,out_varTest))$class
    }
  }
  else if (classifier=="knn") {
    if (out.dev==FALSE) {
      predictions <- knn(out_avrTrain,out_avrTest,y,k=floor(sqrt(pTrain)))
    }
    else {
      predictions <- knn(cbind(out_avrTrain,out_varTrain), cbind(out_avrTest,out_varTest),y,k=floor(sqrt(pTrain)))
    }
  }
  
  return(list(predictions = predictions,
              out_avrTrain = out_avrTrain, out_avrTest = out_avrTest,
              out_varTrain = out_varTrain, out_varTest = out_varTest))
}

