MoransI <- function(x, phy=my_tree, trait=NULL, sig = TRUE, sig.t = 0.05, MI.t = 0.05){
  pvr <- x@Eigen$vectors
  model_pre_eigenvectors <- NULL
  tmpTrait <- trait
  model_eigenvectors <- NULL
  model_EV <- NULL
  W <- 1/cophenetic(my_tree)
  diag(W) <- 0
  
  
  for (k in 1:(ncol(pvr)-1)){
    all_EV <- colnames(pvr)
    tmpRes <- matrix(nrow = length(all_EV), ncol = 2) 
    
    for (i in 1:length(all_EV)){
      new_eigenvector <- pvr[,all_EV[i]]
      
      if (is.null(model_pre_eigenvectors)) {
        model_eigenvectors <- new_eigenvector
      } else {
        model_eigenvectors <- cbind(model_pre_eigenvectors, new_eigenvector)
      }
      
      tmpLM <- lm(tmpTrait ~ model_eigenvectors)
      tmpMI <- Moran.I(tmpLM$residuals, weight = W, alternative = "two.sided")
      tmpRes[i, 1] <- round(tmpMI$p.value, 5)
      tmpRes[i, 2] <- tmpMI$observed
    }
    
    idx <- which(tmpRes[,2] == min(tmpRes[,2]))
    
    selected_EV <- all_EV[idx]
    selected_eigenvector <- pvr[,selected_EV]
    all_EV <- all_EV[-idx]
    pvr <- pvr[,!colnames(pvr) %in% selected_EV]
    
    if (is.null(model_pre_eigenvectors)) {
      
      model_pre_eigenvectors <- selected_eigenvector
      model_EV <- selected_EV
    } else {
      model_pre_eigenvectors <- cbind(model_pre_eigenvectors, selected_eigenvector)
      model_EV <- c(model_EV, selected_EV)
    } 
    
    if (tmpRes[idx,1] > 0.05 | (tmpRes[idx,1] < 0.05 & tmpRes[idx,2] < 0)){
      break
    }
    
    k <- k + 1
    
  }
  
  #model_pre_eigenvectors
  return(model_EV)
  
  
}