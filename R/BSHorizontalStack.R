BSHorizontalStack <- function(T = 100, mtry = NA, nodesize = 5, iter = 100, Xn = NULL, ECHO = TRUE, Cf = NULL, Y, X1, X2, ...){
  if(is.null(Xn)){
    Xl <- list(...)
    Nstack <- 2 + length(Xl)

    Xn <- list()
    Xn[[1]] <- X1
    Xn[[2]] <- X2

    if (Nstack > 2){
      for(n1 in 3:Nstack){
        Xn[[n1]] <- Xl[[n1-2]]
      }
    }
  }
  Nstack <- length(Xn)

  # Set of samples common in all datasets (has full feature set)
  sNover <- intersect(rownames(Xn[[1]]),rownames(Xn[[2]]))

  # Calculate the length of each stacked group and update intersection
  Nl = matrix(0,Nstack,1)

  # Assigne parameters for individual RFs
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

  for (l in 1:Nstack){
    sNover <- intersect(sNover, rownames(Xn[[l]]))
    Nl[l] = nrow(Xn[[l]])
    if(is.null(mtry)){
      mtryS[l] <- max(floor(ncol(Xn[[l]])/3),1)
    }
  }

  if(is.null(Cf)){
    # all feature names
    sNover <- NULL

    # Save the feature names of each dataset
    fNames = list()
    for (l in 1:Nstack){
      sNover <- intersect(sNover, colnames(Xn[[l]]))
      Nl[l] = nrow(Xn[[l]])
    }
  }


  # If validation set specified use that set, otherwise use set found.
  if(is.null(Cf)){
    Cf <- sNover
  }

  # Warning and errors for low number of samples.
  if (length(Cf) > 0){
    if(ECHO) {print(paste(length(Cf), " Overlapping samples found."))}
  } else {
    stop("No overlapping samples found, can not stack these datasets.")
  }

  if (length(Cf) < Nstack*15){
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

    # Build an individual model using sampled data.
    M <- list()

    for(l in 1:Nstack){
      # Take Bootstrap sampled data and non-common samples.
      Xil <- rbind(Xn[[l]][BSCF,], Xn[[l]][!rownames(Xn[[l]]) %in% Cf, ])
      Yil <- rbind(Y[BSCF,1, drop=FALSE], Y[setdiff(rownames(Xn[[l]]),Cf),1, drop=FALSE])
      M[[l]] <- randomForest(Xil,Yil[,1], xtest=NULL, ytest=NULL,ntree=T, mtry=mtryS[l], nodesize = nodesizeS[l])
    }

    # Predict with non-sampled data.
    P <- matrix(1,Nv,Nstack+1)


    for (l in 1:Nstack){
      P[,l+1] <- predict(M[[l]], Xn[[l]][NotBSCF,])
    }


    LM <- glm.fit(P, Yv, family = gaussian())

    # Return Coefficients
    LM$coefficients
  }

  # Final value is the average of all the weights.
  W <- rowMeans(witer)
  return(W)
}
