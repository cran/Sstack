#' Vertical stacking Random Forest models.
#'
#' Generate the weights for a vertically stacked set of Random Forest (RF) models given a set of heterogeneous datasets.
#' For vertical stacking at least one dataset must contain full record (all features). Subfunction of BSstack but can be used stand-alone.
#'
#' Required Packages:
#'   dplyr, randomForest, foreach
#'
#' @param T Number of trees for the individual RF models. (int)
#' @param mtry Number of variables available for splitting at each tree node. If a scalar is given then all models use the given values. If a 1D array is given then each individual model uses the given value.
#' @param nodesize Minimum size of terminal nodes. If a scalar is given then all models use the given values. If a 1D array is given then each individual model uses the given value. By default all models use 5.
#' @param iter The number of time to bootstrap sample the data. (int)
#' @param ECHO Bool, enable to provide output to the user in terms of overlapping samples and runtime for CV.
#' @param Xn List containing each dataset to be stacked. If not supplied will be generated from X1, X2, ...
#' @param Xfull Data table containing samples with full record. Used for generating the weights. Will attempt to find if not given.
#' @param Y  Nsample x 1 data table of responses for ALL samples. Must have matching rownames with each individual dataset.
#' @param X1 Data table of first dataset to be stacked. Rownames should be contained within Y.
#' @param X2 Data table of second dataset to be stacked. Rownames should be contained within Y.
#' @param ... Further data tables, X3, X4, ..., Xl.
#'
#' @return Weights and offsets for each individual RF model.
#'
#' @export

BSVerticalStack <- function(T = 50, mtry=NULL, nodesize=5, iter = 25, ECHO = TRUE, Y, Xfull = NULL, Xn = NULL, X1, X2, ...){

  # If not given assign all datasets to Xn.
  if(is.null(Xn)){
      Xl <- list(...)
      Nstack <- 2 + length(Xl) # Possible Number of stacked groups
      Xn <- list()
      Xn[[1]] <- X1
      Xn[[2]] <- X2

      if(Nstack > 3){
        for(l in 3:Nstack){
          Xn[[l]] <- Xl[[l-2]]
        }
      }

    } else {Nstack <- length(Xn)}

    if(is.null(mtry)){
      mtryS <- rep(1,times = Nstack)
      mtryS[1] <- max(floor(ncol(Xn[[1]])/3),1)
      mtryS[2] <- max(floor(ncol(Xn[[2]])/3),1)
    } else if (length(mtry) == 1){
      mtryS <- rep(mtry,times=Nstack)
    }
    else if (!is.vector(mtry) | length(mtry) != Nstack){
      stop("mtry must be either a scalar or vector with size = number of databases.")
    }
    else{
      mtryS <- mtry
    }

    # Assign nodesize values for each individual RF
    if(is.null(nodesize)){
      nodesizeS <- rep(5,times = Nstack)
    } else if (length(nodesize) == 1){
      nodesizeS <- rep(nodesize,times=Nstack)
    }
    else if (!is.vector(nodesize) | length(nodesize) != Nstack){
      stop("nodesize must be either a scalar or vector with size = number of databases.")
    }
    else{
      nodesizeS <- nodesize
    }

    Nstack <- length(Xn)

    # All feature names
    fNall <- NULL

    # Number of features in each set
    Nl = matrix(0,Nstack,1)

    # Save the feature names of each dataset
    fNames = list()
    for (l in 1:Nstack){
      fNall <- union(fNall, colnames(Xn[[l]]))
      fNames[[l]] = colnames(Xn[[l]])
      Nl[l] = nrow(Xn[[l]])
      if(is.null(mtry)){
        mtryS[l] <- max(floor(ncol(Xn[[l]])/3),1)
      }
    }


  if(is.null(Xfull)){
  # Find dataset(s) with full record and save them for BS sampling
  nOverl <- 0
  Xfull <- data.frame()

  for (l in 1:Nstack){
    if(setequal(fNall, fNames[[l]])){
      rNames <- rownames(Xfull)
      Xfull <- bind_rows(Xfull,Xn[[l]])

      # Perserve rownames
      if(nOverl > 0){
        rownames(Xfull)[1:nOverl] <- rNames
      }
      rownames(Xfull)[(nOverl+1):(nOverl+Nl[[l]])] <- rownames(Xn[[l]]) # Preserve colnames
      nOverl <- nOverl + Nl[[l]]
    }
  }
  }

  Cf <- rownames(Xfull)
  nOverl <- length(Cf)

  # Warning and errors for low number of samples.
  if (nOverl > 0){
    if( ECHO ) {print(paste(nOverl, " Overlapping samples found."))}
  } else {
    stop("No overlapping samples found, can not stack these datasets.")
  }

  if (nOverl < Nstack*15){
    warning("Small number of overlapping samples found, stacking may provide little to no benefit.")
  }


  # Parallel foreach loop to get a set of weights.
  witer<-foreach(n1=1:iter,.combine=cbind, .packages='randomForest') %dopar%{
    # Get a set of bootstrap samples.
    BS <- sample(length(Cf), length(Cf), replace = TRUE)
    BSCF <- Cf[BS]
    NotBSCF <- setdiff(Cf, BSCF) # Non-picked samples used as validation set
    Nv <- length(NotBSCF)
    Yv <- Y[NotBSCF,1]

    M <- list()

    for (l in 1:Nstack){
      # Take Bootstrap sampled data and non-common samples.
      Xil <- rbind(Xn[[l]][BSCF,], Xn[[l]][!rownames(Xn[[l]]) %in% Cf, ])
      Yil <- rbind(Y[BSCF,1, drop=FALSE], Y[setdiff(rownames(Xn[[l]]),Cf),1, drop=FALSE])
      M[[l]] <- randomForest(Xil[complete.cases(Xil),],Yil[complete.cases(Xil),1], xtest=NULL, ytest=NULL,ntree=T, mtry=mtryS[l], nodesize = nodesizeS[l])
    }

    # Predictions for each individual model (first column for constant offset)
    P <- matrix(1,Nv,Nstack+1)

    for (l in 1:Nstack){
      P[,l+1] <- predict(M[[l]],Xfull[NotBSCF,fNames[[l]]])
    }

    # Fit predictions to a glm model
    LM <- glm.fit(P, Yv, family = gaussian())

    # Return Coefficients
    LM$coefficients
  }

  # Final value is the average of all the weights.
  W <- rowMeans(witer)
  return(W)
}
