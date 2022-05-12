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


dirOutClass <- function(test, train, y, outly.method="SDO", classifier="knn", L.norm=2, depth=FALSE, out.var=TRUE) {
 
library(mrfDepth) # AO
library(e1071) # svm
library(MASS) # qda
library(class) # knn

dirOut <- function(test, train, outly.method, L.norm, depth, abs) {

  nTrain <- dim(train)[1]
  nTest <- dim(test)[1]
  p <- dim(train)[2]
  
  ### Univariate method for outlyingness calculation
  if (outly.method=="Mahalanobis") {
    Zcenter <- apply(train,2,mean)
    Zsd <- apply(train,2,sd)
    dirOutMatrix <- (test-rep(Zcenter,each=nTest))/rep(Zsd,each=nTest)
  }
  else if (outly.method=="SDO") {
    Zcenter <- apply(train,2,median)
    Zmad <- apply(train,2,mad)
    dirOutMatrix <- (test-rep(Zcenter,each=nTest))/rep(Zmad,each=nTest) 
    }
  else if (outly.method=="AO"){
    Zcenter <- apply(train,2,median)
    adjOutlCol <- apply(test,2,function(columna) {adjOutl(columna[1:nTrain],columna[(nTrain+1):nTest])})
    dirOutMatrix <- sapply(1:p,function(ncol) {c(adjOutlCol[[ncol]]$outlyingnessX,adjOutlCol[[ncol]]$outlyingnessZ)})
    dirOutMatrix <- dirOutMatrix*sign(test-rep(Zcenter,each=nTest))
    }
  
  
  ### False for outlyingness, True for depth
  if (depth==TRUE) {
    dirOutMatrix <- 1/(1+abs(dirOutMatrix))
  }

  out_avr <- apply(dirOutMatrix, 1, mean)
  
  ### L norm to calculate variation
  if (L.norm=="Inf") {
    out_var <- apply(matrix(apply(dirOutMatrix, 2, (function(columna) abs(columna-out_avr))),nTest),1,max)
  }
  else {
    out_var <- (apply(matrix(apply(dirOutMatrix, 2, (function(columna) (abs(columna-out_avr))**L.norm)),nTest),1,mean))**(1/L.norm)
  }
  
  ### True for the absolute value of out_avr
  if (abs==TRUE) {
    dirOutMatrix <- abs(dirOutMatrix)
  }
  out_avr <- apply(dirOutMatrix, 1, mean)
  
  M <- cbind(out_avr,out_var)
  colnames(M) <- c("out_avr","out_var")
  return(M)
}


dirOutClass <- function(test, train, target, outly.method="Mahalanobis", L.norm=2, abs=FALSE, depth=FALSE, classifier="qda", classify=TRUE, out.avr=TRUE, out.var=TRUE) {
  
  train <- as.matrix(train)
  test <- as.matrix(test)
  nTrain <- dim(train)[1]
  nTest <- dim(test)[1]
  trainTestBind <- rbind(train,test)
  
  if (classify==FALSE) {
    
    ### If classify==FALSE, only calculate the Directional Outlyingness
    predictions <- NA
    M <- dirOut(trainTestBind,train, outly.method=outly.method, L.norm=L.norm, depth=depth, abs=abs)
    out_avrTrain <- M[1:nTrain,1]
    out_avrTest <- M[(nTrain+1):(nTrain+nTest),1]
    out_varTrain <- M[1:nTrain,2]
    out_varTest <- M[(nTrain+1):(nTrain+nTest),2]
  }
  
  else {
    
    ### Dimensionality reduction using Directional Outlyingness
    levels <- levels(target)
    out_avrTrain <- matrix(NA, nTrain, length(levels))
    out_avrTest <- matrix(NA, nTest, length(levels))
    out_varTrain <- matrix(NA, nTrain, length(levels))
    out_varTest <- matrix(NA, nTest, length(levels))
    
    for (i in 1:length(levels)) {
      M <- dirOut(trainTestBind,train[target==levels[i],], outly.method=outly.method, L.norm=L.norm, depth=depth, abs=abs)
      out_avrTrain[1:nTrain,i] <- M[1:nTrain,1]
      out_avrTest[1:nTest,i] <- M[(nTrain+1):(nTrain+nTest),1]
      out_varTrain[1:nTrain,i] <- M[1:nTrain,2]
      out_varTest[1:nTest,i] <- M[(nTrain+1):(nTrain+nTest),2]
    }
    
    ### Classification
    ### If out.avr==FALSE or out.avr==FALSE, the method won't use those values
    if (classifier=="svm") {
      if (out.avr==FALSE) {
        model <- svm(out_varTrain,target)
        predictions <- predict(model, out_varTest)
      }
      else if (out.var==FALSE) {
        model <- svm(out_avrTrain,target)
        predictions <- predict(model, out_avrTest)
      }
      else {
        model <- svm(cbind(out_avrTrain,out_varTrain),target)
        predictions <- predict(model, cbind(out_avrTest,out_varTest))
      }
    }
    
    else if (classifier=="qda") {
      if (out.avr==FALSE) {
        model <- qda(out_varTrain,target)
        predictions <- predict(model, out_varTest)$class
      }
      else if (out.var==FALSE) {
        model <- qda(out_avrTrain,target)
        predictions <- predict(model, out_avrTest)$class
      }
      else {
        model <- qda(cbind(out_avrTrain,out_varTrain),target)
        predictions <- predict(model, cbind(out_avrTest,out_varTest))$class
      }
    }
    
    else if (classifier=="knn") {
      if (out.avr==FALSE) {
        predictions <- knn(out_varTrain,out_varTest,target,k=floor(sqrt(nTrain)))
      }
      if (out.var==FALSE) {
        predictions <- knn(out_avrTrain,out_avrTest,target,k=floor(sqrt(nTrain)))
      }
      else {
        predictions <- knn(cbind(out_avrTrain,out_varTrain), cbind(out_avrTest,out_varTest),target,k=floor(sqrt(nTrain)))
      }
    }
  }

  return(list(predictions = predictions,
              out_avrTrain = out_avrTrain, out_avrTest = out_avrTest,
              out_varTrain = out_varTrain, out_varTest = out_varTest))
}
