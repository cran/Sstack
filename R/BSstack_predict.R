#' Predict using a set of Stacked Random Forest models.
#'
#' Gives predictions for a linear bootstrapped stacked Random Forest predictors. Gives the predictions
#' of each individual model as well as the linearlly combined predictions.
#'
#' Required Packages:
#'   randomForest
#'
#' @param BSmodel List containing the individual Random Forest models, their weights, and feature names. Generated using BSstack function.
#' @param Xi NxM datatable containing input features to be predicted. Must contain all features used in the individual RF models.
#'
#' @return NxL+1 matrix where L is the number of individual RF models. Predictions for the ith RF model is found in the ith column of this matrix while predictions for the stacked model is in the final column.
#'
#' @export

BSstack_predict <- function(BSmodel, Xi){
  if(length(BSmodel) == 7){
    BSmodel = BSmodel[[1]]
  }

  N <- nrow(Xi)
  w <- BSmodel[[1]]
  RFn <- BSmodel[[2]]
  fNames <- BSmodel[[3]]

  Nstack <- length(RFn)

  P <- matrix(data = 0, nrow = N, ncol = (Nstack+1))

  for(l in 1:Nstack){
    Xp <- Xi[,fNames[[l]]]

    if(anyNA(Xp)){
      stop("ERROR: Input samples missing features.")
    }


    P[,l] <- predict(RFn[[l]], Xp)
  }

  if(N==1){
    P[,Nstack+1] <- c(1, P[1:Nstack]) %*% w
  }
  else{
    P[,Nstack+1] <- cbind(matrix(1,N,1),P[1:N,1:Nstack]) %*% w
  }

  return(P)
}
